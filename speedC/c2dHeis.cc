//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

// calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[])
{
  	std::clock_t tStart = std::clock();

  	if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto Ly = input.getInt("Ly", 3);
	auto Lx = input.getInt("Lx", 8);
    auto h = input.getReal("h", 1.);
    auto truncE = input.getReal("truncE", 1E-8);
    auto maxDim = input.getInt("maxDim", 512);
	auto GSETDVP = input.getYesNo("GSETDVP",true);
	println();

  	// We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.2f_maxDim_%d_c2dHeis.dat", Ly, Lx, h, maxDim);

    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "maxBondDim" << " " << "enPhi" << " " << "svn(x)" << " " << std::endl;

  	auto N = Ly * Lx;
    auto sites = SpinHalf(N);
  	auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", false});

	auto ampo = AutoMPO(sites);
    for(auto j : lattice){
        ampo += 0.5, "S+", j.s1, "S-", j.s2;
        ampo += 0.5, "S-", j.s1, "S+", j.s2;
        ampo += 1.0, "Sz", j.s1, "Sz", j.s2;
    }
    /////// Checkerboard Geometry //////
    int col = 1;
    for(auto j : range1(N)){
        if(Ly%2!=0){
            if(j%2==0)
                ampo += +h, "Sz", j;
            else
                ampo += -h, "Sz", j;
        }
        else{
            if(j%2==0)
                ampo += +h * pow(-1., col), "Sz", j;
            else
                ampo += -h * pow(-1., col), "Sz", j;
            if(j%Ly == 0)
                col++;    
        }
    }
    auto H = toMPO(ampo);

	// initial state
	auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i, i%2==0 ? "Up" : "Dn");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

	// 2d heisenberg model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, 100, 200, 200, 400, 400, 800, maxDim;
    sweeps.cutoff() = truncE;
    sweeps.noise() = 1E-7,1E-8,0.0;

	// calculate ground state
	auto [enPsi, psi] = dmrg(H, initState, sweeps, {"Silent=", true});

	// make |phi> = Sz|psi>
	int loc = (Lx / 2 - 1) * Ly + 1; 
	psi.position(loc);
	auto newA = 2.0 * sites.op("Sz", loc) * psi(loc);
	newA.noPrime();
	psi.set(loc, newA);
    psi.orthogonalize({"Cutoff=",truncE,"MaxDim=",maxDim});
    auto enPhi = inner(psi,H,psi); //energy after disturbing the ground state

	// calculate von neumann entanglement entropy
	std::vector<double> svn0(Lx-1), svn(Lx-1);
	for(int j = 1; j < Lx; j++){
		auto b = j*Ly;
		svn0[j-1] = vonNeumannS(psi, b);
	}

	// store to file
	dataFile << 0 << " " << enPsi << " " << maxLinkDim(psi) << " " << enPhi << " ";
	for (int j = 0; j < Lx-1; j++){
		dataFile << svn0[j] << " ";
	}
	dataFile << std::endl;

	printfln("time = %0.1f; phi energy = %0.3f, max link dim is %d", 0, enPhi, maxLinkDim(psi));

	// time evolution parameters.
    double tval = 0., dt = 0.1;
    double delta1 =  0.414490771794376*dt; // for 4th order TDVP
    double delta2 = -0.657963087177503*dt;
    double finalTime = double(Lx)/3.14159;
    int nt=int(finalTime/dt);

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxDim;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 10;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxDim;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 10;

    println("\nStarting 4th GSE-TDVP\n");

	////////////////////////////////////////////////////////////////////////////
	///////// time evolve //////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	for (int n = 1; n <= nt; n++){
    	tval += dt; // update time

		std::clock_t tStartTDVP = std::clock(); // check time for performance
		if(GSETDVP){
            // time evolve with GSE-TDVP
            std::vector<int> dimK = {maxLinkDim(psi), maxLinkDim(psi)};
            addBasis(psi, H, dimK, {"Cutoff",truncE,
                                            "Method", "DensityMatrix",
                                            "KrylovOrd",3,
                                            "Quiet",true});
            // check if bond dimension has grown enough
            if(maxLinkDim(psi)>=maxDim){
                GSETDVP = false;
                printfln("\n --- Starting 2-TDVP --- ");
            }
            // one-site TDVP
            tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1});
            tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"Truncate",true,"NumCenter",1});
            enPhi = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1});
        }
        else{
			if(n==1)
				printfln("\n --- Starting 2-TDVP --- ");
            // two-site TDVP
            tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",2});
            tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"Truncate",true,"NumCenter",2});
            enPhi = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",2});
        }

        auto tdvpTime = (double)(std::clock() - tStartTDVP)/CLOCKS_PER_SEC;    
        
        // calculate svn
        for(int j = 1; j < Lx; j++){
			auto b = j*Ly;
			svn[j-1] = vonNeumannS(psi, b);
		}

		// store to file
		dataFile << tval << " " << enPsi << " " << maxLinkDim(psi) << " " << enPhi << " " ;
		for(int j = 0; j<Lx-1; j++){ //save svn
			dataFile << svn[j]-svn0[j] << " ";
		}
		dataFile << std::endl;

      	printfln("\n----\n Iteration %d, time = %0.2f; phi energy = %0.3f, maxDim = %d, tdvp time = %0.3fs\n----\n", n, tval, enPhi, maxLinkDim(psi),tdvpTime);


  	} // for n

	dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

} // main

//calculate entanglement
Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //choose orthogonality center and perform svd
    psi.position(b);
    auto l = leftLinkIndex(psi,b);
    auto s = siteIndex(psi,b);
    auto [U,S,V] = svd(psi(b),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    for(auto n : range1(dim(u))){
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    return SvN;

}//vonNeumannS
