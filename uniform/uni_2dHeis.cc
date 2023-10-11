//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>);
// calculate energy for periodic boundary at x = i
double calculatePBCenergy(int, int, MPS, SiteSet, ITensor, ITensor, ITensor);
// calculate energy of horizontal bonds for one column
std::vector<double> calculateLRenergy(int, int, MPS, SiteSet, 
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){
        printfln("Usage: %s input_file",argv[0]);
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto Lx = input.getInt("Lx", 12);
    auto Ly = input.getInt("Ly", 3);
    auto h = input.getReal("h", 5.0);
    auto tau = input.getReal("tau",1.0);
    auto truncE = input.getReal("truncE", 1E-8);
    auto maxDim = input.getInt("maxDim", 128);
    auto tanhshift = input.getReal("tanhshift",4.0);
    auto dt = input.getReal("dt",0.1);
    auto GSETDVP = input.getYesNo("GSETDVP",true);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_tau_%0.1f_maxDim_%d_gse_%d_2dHeis_uni.dat",Ly,Lx,h,tau,maxDim,GSETDVP);

    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    datafile << "time" << " " << "en(t)" << " " << "enf(t)" << " " << "enf(t)-en0" << " " << "svN(t)" << " "
             << "localEn(t)" << " " << "localEn(t)-localEn0" << " " << "SzSz(t)" << " " << "Sperp(t)" << std::endl;

    auto N = Ly * Lx;
    auto sites = SpinHalf(N);

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += 0.5, "S+", j.s1, "S-", j.s2;
        ampo += 0.5, "S-", j.s1, "S+", j.s2;
        ampo += 1.0, "Sz", j.s1, "Sz", j.s2;
    }
    auto Hfinal = toMPO(ampo);
    // H at t=0 (gapped ground state)
    auto hval = h*(0.5 + 0.5*tanh( tanhshift ));
    /////// Checkerboard Geometry //////
    int col = 1;
    for(auto j : range1(N)){
        if(Ly%2!=0){
            if(j%2==0)
                ampo += +hval, "Sz", j;
            else
                ampo += -hval, "Sz", j;
        }
        else{
            if(j%2==0)
                ampo += +hval * pow(-1., col), "Sz", j;
            else
                ampo += -hval * pow(-1., col), "Sz", j;
            if(j%Ly == 0)
                col++;
        }
    }
    auto H = toMPO(ampo);
    
    //initial state
    auto state = InitState(sites); 
    for(int i = 1; i <= N; i++){
        state.set(i, i%2==0 ? "Up" : "Dn");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    // 2d heisenberg model parameters
    auto sweeps = Sweeps(15);
    sweeps.maxdim() = 20, 50, 100, 100, 200, 200, 400, 400, 800, maxDim;
    sweeps.cutoff() = truncE;
    sweeps.noise() = 1E-6, 1E-7, 1E-8, 0;

    // calculate initial local energy density
    std::vector<double> localEn0(Lx-1,0.0), localEn(Lx-1,0.0); // local energy density vector
    std::vector<double> szsz(Lx,0.0), sperpsperp(Lx,0.0);

    //make 2D vector of ITensor for local energy operators
    //long-range interactions have the same structure as nearest-neighbour when we use swap gates
    std::vector<std::vector<ITensor>> PM(Lx, std::vector<ITensor>(Ly-1));
    std::vector<std::vector<ITensor>> MP(Lx, std::vector<ITensor>(Ly-1));
    std::vector<std::vector<ITensor>> ZZ(Lx, std::vector<ITensor>(Ly-1));

    std::vector<std::vector<ITensor>> PM_LR(Lx-1, std::vector<ITensor>(Ly));
    std::vector<std::vector<ITensor>> MP_LR(Lx-1, std::vector<ITensor>(Ly));
    std::vector<std::vector<ITensor>> ZZ_LR(Lx-1, std::vector<ITensor>(Ly));

    // make local energy tensors
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<=Ly; j++){
            int index = (i-1)*Ly+j;
            //MPS long-range
            if(i<Lx && j==1){
                for(int m = 0; m<Ly; m++){
                    PM_LR[i-1][m] = 0.5*sites.op("S+",index+2*m)*sites.op("S-",index+2*m+1);
                    MP_LR[i-1][m] = 0.5*sites.op("S-",index+2*m)*sites.op("S+",index+2*m+1);
                    ZZ_LR[i-1][m] =     sites.op("Sz",index+2*m)*sites.op("Sz",index+2*m+1);                    
                }
            }
            // MPS nearest-neighbour
            if(j<Ly){
                PM[i-1][j-1] = 0.5*sites.op("S+",index)*sites.op("S-",index+1);
                MP[i-1][j-1] = 0.5*sites.op("S-",index)*sites.op("S+",index+1);
                ZZ[i-1][j-1] =     sites.op("Sz",index)*sites.op("Sz",index+1);
            }
        }
    }

    //DMRG to find critical ground state
    auto [en0,psi] = dmrg(Hfinal,initState,sweeps,{"Silent=",true});
    auto dim0 = maxLinkDim(psi); //defines max dimension of critical state
    localEn0 = calculateLocalEnergy(Lx, Ly, sites, psi, PM, MP, ZZ, PM_LR, MP_LR, ZZ_LR);

    //DMRG to find ground state at t=0
    auto en = dmrg(psi,H,sweeps,{"Silent=",true}); // doing this saves on memory
    auto enf = inner(psi, Hfinal, psi);

    // calculate von Neumann S
    std::vector<double> svN(Lx-1,0.0);
    for(auto j : range1(Lx-1)){
        svN[j-1] = vonNeumannS(psi, j*Ly);
    }
    // calculate local energy density
    localEn = calculateLocalEnergy(Lx, Ly, sites, psi, PM, MP, ZZ, PM_LR, MP_LR, ZZ_LR);
    // calculate spin auto-correlators
    for(int b = 1; b<=Lx; b++){
        auto [zz, perp] = spinspin( (Lx-1)/2*Ly+1, (b-1)*Ly + 1, psi, sites);
        szsz[b-1] = zz;
        sperpsperp[b-1] = perp;
    }

    // store data to file
    datafile << 0.0 << " " << en << " " << enf << " " << enf-en0 << " ";
    for (int j=0; j < Lx-1; j++){
        datafile << svN[j] << " ";
    }
    for (int j=0; j < Lx-1; j++){
        datafile << localEn[j] << " ";
    }
    for (int j=0; j < Lx-1; j++){
        datafile << localEn[j]-localEn0[j] << " ";
    }
    for (int j = 0; j < Lx; j++){
        datafile << szsz[j] << " ";
    }
    for (int j = 0; j < Lx; j++){
        datafile << sperpsperp[j] << " ";
    }
    datafile << std::endl;

    // time evolution parameters
    double tval = 0.0; //time
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 12*tau;
    int nt = int(finalTime/dt);

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxDim;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 10;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxDim;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 10;

    printfln("\ncritical GS: energy = %0.3f, maxDim = %d\n", en0, dim0);
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d\n", tval, en, svN[(Lx+1)/2], maxLinkDim(psi));

    ////////////////////////////////////////////////////////////////////////////
    ///////// time evolve //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    printfln("\n --- Starting GSE-TDVP --- ");
    auto stopCondition = false;

    for(int n=1; n<=nt && !stopCondition; n++){

        tval += dt; //update time vector

	    // update magnetic field
        hval = h*(0.5 + 0.5*tanh( -tval/tau + tanhshift ));

        ampo = AutoMPO(sites);
        // autompo hamiltonian
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
                    ampo += +hval, "Sz", j;
                else
                    ampo += -hval, "Sz", j;
            }
            else{
                if(j%2==0)
                    ampo += +hval * pow(-1., col), "Sz", j;
                else
                    ampo += -hval * pow(-1., col), "Sz", j;
                if(j%Ly == 0)
                    col++;
            }
        }
        H = toMPO(ampo);

        std::clock_t tStartTDVP = std::clock(); // check time for performance
        
        if(GSETDVP){
            // time evolve with GSE-TDVP
            //std::vector<int> dimK = {maxLinkDim(psi), maxLinkDim(psi)};
            std::vector<Real> truncK = {0.1*truncE, 0.1*truncE};
            //addBasis(psi, H, dimK, {"Cutoff",truncE,
            addBasis(psi, H, truncK, {"Cutoff",truncE,
                                            "Method", "DensityMatrix",
                                            "KrylovOrd",3,
                                            "Quiet",true});
            // check if bond dimension has grown enough
            if(maxLinkDim(psi)>=maxDim){
                GSETDVP = false;
                printfln("\n --- Starting 2-TDVP at t = %0.1f --- ", tval+dt);
            }
            // one-site TDVP
            tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1,"ErrGoal",truncE});
            tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"Truncate",true,"NumCenter",1,"ErrGoal",truncE});
            en = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1,"ErrGoal",truncE});
            if(maxLinkDim(psi)>dim0){
                GSETDVP = false;
                printfln("\n --- Starting 2-TDVP at t = %0.2f --- ", tval+dt);
            }
        }
        else{
            if(n==1)
				printfln("\n --- Starting 2-TDVP --- ");
            // two-site TDVP
            tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",2,"ErrGoal",truncE});
            tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"Truncate",true,"NumCenter",2,"ErrGoal",truncE});
            en = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",2,"ErrGoal",truncE});
        }

        auto tdvpTime = (double)(std::clock() - tStartTDVP)/CLOCKS_PER_SEC;

        //calculate energy wrt final Hamiltonian
        enf = innerC(psi, Hfinal, psi).real();
        //calculate entanglement entropy
        for(auto j : range1(Lx-1)){
            svN[j-1] = vonNeumannS(psi, j*Ly);
        }
        // calculate local energy density <psi(t)|H(x,y)|psi(t)>
        localEn = calculateLocalEnergy(Lx, Ly, sites, psi, PM, MP, ZZ, PM_LR, MP_LR, ZZ_LR);
        // calculate spin auto-correlators
        for(int b = 1; b<=Lx; b++){
            auto [zz, perp] = spinspin( (Lx-1)/2*Ly+1, (b-1)*Ly + 1, psi, sites);
            szsz[b-1] = zz;
            sperpsperp[b-1] = perp;
        }

        // store data to file
        datafile << 0.0 << " " << en << " " << enf << " " << enf-en0 << " ";
        for (int j=0; j < Lx-1; j++){
            datafile << svN[j] << " ";
        }
        for (int j=0; j < Lx-1; j++){
            datafile << localEn[j] << " ";
        }
        for (int j=0; j < Lx-1; j++){
            datafile << localEn[j]-localEn0[j] << " ";
        }
        for (int j = 0; j < Lx; j++){
            datafile << szsz[j] << " ";
        }
        for (int j = 0; j < Lx; j++){
            datafile << sperpsperp[j] << " ";
        }
        datafile << std::endl;

        printfln("\nt = %0.2f, enf-en0 = %0.3g, SvN = %0.3f, maxDim = %d, wall time = %0.3fs", tval, enf-en0, svN[(Lx+1)/2], maxLinkDim(psi), tdvpTime);

        if( abs(en - enf) < 1E-5){
            stopCondition = true;
            printfln("stop condition met, |en-enf| = %0.10f", abs(en-enf));
        }

    }

    datafile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

// calculates local energy for 2D MPS using gates
std::vector<double> calculateLocalEnergy(int Lx, int Ly, SiteSet sites, MPS psi,
                                std::vector<std::vector<ITensor>> PM,
                                std::vector<std::vector<ITensor>> MP,
                                std::vector<std::vector<ITensor>> ZZ,
                                std::vector<std::vector<ITensor>> PM_LR,
                                std::vector<std::vector<ITensor>> MP_LR,
                                std::vector<std::vector<ITensor>> ZZ_LR){

    std::vector<double> tempEn(Lx,0.0);

    std::vector<double> localEnergy(Lx-1, 0.0); // interpolated energy density

    // MPS nearest-neighbour interaction
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<Ly; j++){ 

            int index = (i-1)*Ly + j;
            psi.position(index);
            auto ket = psi(index)*psi(index+1);

            tempEn[i-1] += eltC(dag(prime(ket,"Site")) * PM[i-1][j-1] * ket).real();
            tempEn[i-1] += eltC(dag(prime(ket,"Site")) * MP[i-1][j-1] * ket).real();
            tempEn[i-1] += eltC(dag(prime(ket,"Site")) * ZZ[i-1][j-1] * ket).real();

        }// for j
    }// for i

    // y-periodic interactions
    for(int i=1; i<=Lx; i++){

        tempEn[i-1] += calculatePBCenergy(i, Ly, psi, sites, PM[i-1][Ly-2], MP[i-1][Ly-2], ZZ[i-1][Ly-2]);

    }// for i

    // MPS long-range interactions
    for(int i=1; i<Lx; i++){

        int index = (i-1)*Ly + 1;
        auto lrEnergy = calculateLRenergy(i, Ly, psi, sites, PM_LR, MP_LR, ZZ_LR);

        for(int m = 0; m<Ly; m++){ // long-range

            localEnergy[i-1] += lrEnergy[m];

        }// for m
    }// for i

    // interpolate tempEn into localEnergy
    for(int i=1; i<Lx; i++){

        // interpolation and average of the nearest-neighbour interactions
        if( i==1 ){ // left boundary
            localEnergy[i-1] += tempEn[i-1] + 0.5*tempEn[i];
        }
        else if( i==Lx-1 ){ // right boundary
            localEnergy[i-1] += 0.5 * tempEn[i-1] + tempEn[i];
        }
        else{
            localEnergy[i-1] += 0.5 * (tempEn[i-1] + tempEn[i]);
        }
    } // for i

    return localEnergy;

}//localEnergy

double calculatePBCenergy(int i, int Ly, MPS psi, SiteSet sites, ITensor PM, ITensor MP, ITensor ZZ){

    double energy;
    int index = (i-1)*Ly + 1;

    psi.position(index);

    for(int n=0; n<Ly-2; n++){

        int b = index+n;//define gate bond
        auto g = BondGate(sites,b,b+1);
        auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
        AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
        psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
        psi.position(g.i2()); //orthogonality center moves to the right

    } // for n

    auto ket = psi(index+Ly-2)*psi(index+Ly-1);
    energy =  eltC( dag(prime(ket,"Site")) * PM * ket).real();
    energy += eltC( dag(prime(ket,"Site")) * MP * ket).real();
    energy += eltC( dag(prime(ket,"Site")) * ZZ * ket).real();

    return energy;
} // calculatePBCenergy

std::vector<double> calculateLRenergy(int i, int Ly, MPS psi, SiteSet sites,
                                    std::vector<std::vector<ITensor>> PM_LR,
                                    std::vector<std::vector<ITensor>> MP_LR,
                                    std::vector<std::vector<ITensor>> ZZ_LR){

    std::vector<double> energy(Ly,0.);

    int index = (i-1)*Ly + 1;

    // long range interactions with smart ordering of gates 
    for(int m=0; m<=Ly-2; m++){

        psi.position(index+Ly+m);

        for(int n=Ly; n>m+1; n--){

            int b = index + n + m;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi.position(g.i1()); //orthogonality center moves to the left

        } // for n
    } // for m

    for(int m = 0; m<Ly; m++){

        psi.position(index+2*m);
        auto ket = psi(index+2*m)*psi(index+2*m+1);
        energy[m] = eltC( dag(prime(ket,"Site")) * PM_LR[i-1][m] * ket).real();
        energy[m] += eltC( dag(prime(ket,"Site")) * MP_LR[i-1][m] * ket).real();
        energy[m] += eltC( dag(prime(ket,"Site")) * ZZ_LR[i-1][m] * ket).real();

    } // for m

    return energy;
}

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

//calculate spin-spin correlator
std::tuple<double, double> spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrZ, corrPerp;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left 
        }
        auto ket = psi(center)*psi(center+1);
        auto SzSz = sites.op("Sz",center)*sites.op("Sz",center+1);
        auto PM = sites.op("S+",center)*sites.op("S-",center+1);
        auto MP = sites.op("S-",center)*sites.op("S+",center+1);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPerp = eltC( dag(prime(ket,"Site")) * PM * ket).real();
        corrPerp += eltC( dag(prime(ket,"Site")) * MP * ket).real();
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SzSz = sites.op("Sz",center-1)*sites.op("Sz",center);
        auto PM = sites.op("S+",center-1)*sites.op("S-",center);
        auto MP = sites.op("S-",center-1)*sites.op("S+",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPerp = eltC( dag(prime(ket,"Site")) * PM * ket).real();
        corrPerp += eltC( dag(prime(ket,"Site")) * MP * ket).real();
    }
    else{
        corrZ = 1.;
        corrPerp = 1.;
    }

    return {corrZ, corrPerp};

}//SxSx
