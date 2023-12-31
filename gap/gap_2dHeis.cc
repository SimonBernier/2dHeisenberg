//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int Lx=24, Ly=3;

    if(argc > 1)
        Lx = std::stoi(argv[1]);
    if(argc > 2)
        Ly = std::stoi(argv[2]);
    
    printfln("Ly = %d, Lx = %d", Ly, Lx);

    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_2dHeis_gap.dat",Ly,Lx);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "h" << " " << "en0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "gap" << " " << "var1" << " " << "maxBondDim1" << " " << std::endl;

    auto N = Lx * Ly;
    auto sites = SpinHalf(N);

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    //
    // h = 0
    //

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += 0.5, "S+", j.s1, "S-", j.s2;
        ampo += 0.5, "S-", j.s1, "S+", j.s2;
        ampo += 1.0, "Sz", j.s1, "Sz", j.s2;
    }
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i, i%2==0 ? "Up" : "Dn");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    // 2d heisenberg model parameters
    auto sweeps1 = Sweeps(15);
    sweeps1.maxdim() = 20, 50, 100, 100, 200, 200, 400, 400, 800, 1600;
    sweeps1.cutoff() = 1E-8;
    sweeps1.noise() = 1E-7,1E-8,0.0;

    auto sweeps2 = Sweeps(10);
    sweeps2.maxdim() = 20, 50, 100, 100, 200, 200, 400, 400, 800, 1600;
    sweeps2.cutoff() = 1E-8;
    sweeps2.noise() = 0.0;

    //
    //solve for ground state
    //
    auto [en0,psi0] = dmrg(H,initState,sweeps1,{"Quiet=",true});
    auto var = sqrt( abs( inner(H,psi0,H,psi0) - en0*en0 ) );
    auto maxBondDim = maxLinkDim(psi0);
    printfln("\nh = %0.2f", 0.);
    println("first state");
    printfln("Energy = %0.10f, var = %0.6f, maxLinkDim = %d", en0, var, maxBondDim);

    dataFile << 0.0 << " " << en0 << " " << var << " " << maxBondDim << " ";

    //
    // excited state
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(H,wfs,initState,sweeps2,{"Quiet=",true,"Weight=",20.0});
    var = sqrt( abs( inner(H,psi1,H,psi1) - en1*en1 ) );
    maxBondDim = maxLinkDim(psi1);

    println("\nsecond state");
    printfln("Gap = %0.10f, var = %0.6f, maxLinkDim = %d", en1-en0, var, maxBondDim);

    dataFile << en1-en0 << " " << var << " " << maxBondDim << " " << std::endl;

    //
    // h > 0
    //
    double h = 0.;
    double dh = 0.5;
    int Nh = int(9./dh);
    for( int i = 1; i <= Nh; i++){
        h += dh;

        ampo = AutoMPO(sites);
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
        H = toMPO(ampo);

        //
        //solve for ground state
        //
        en0 = dmrg(psi0,H,sweeps2,{"Silent=",true});
        var = sqrt( abs( inner(H,psi0,H,psi0) - en0*en0) );
        maxBondDim = maxLinkDim(psi0);
        printfln("\nh = %0.2f", h);
        println("\nfirst state");
        printfln("Energy = %0.10f, var = %0.10f, maxLinkDim = %d", en0, var, maxBondDim);

        dataFile << h << " " << en0 << " " << var << " " << maxBondDim << " ";

        //
        // excited state
        //
        wfs = std::vector<MPS>(1);
        wfs.at(0) = psi0;
        en1 = dmrg(psi1,H,wfs,sweeps2,{"Silent=",true,"Weight=",20.0});
        var = inner(H,psi1,H,psi1) - en1*en1;
        maxBondDim = maxLinkDim(psi1);

        println("\nsecond state");
        printfln("Gap = %0.10f, var = %0.10f, maxLinkDim = %d", en1-en0, var, maxBondDim);

        dataFile << en1-en0 << " " << var << " " << maxBondDim << " " << std::endl;
    
    }

    dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}// runPT

