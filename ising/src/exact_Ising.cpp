#include <iostream>
#include <fstream>
#include "include/all_exact_Ising.hpp"

using namespace itensor;

int main(int argc, char *argv[])
    {
    
    char output[80];
    std::sprintf(output, "./data/data.dat");
    std::ofstream out_file(output, std::ios::out | std::ios::app);
    out_file.precision(9);
    
    // The system will be a 2x20 ladder
    int N = 8, dN = 0;
    int numg = 1, numh = 1, numN = 1;
    auto J = 4., dg = 0., dh = 0.;
    auto g = 1.*J/2.;
    auto h = 0.5*J/2.;
    
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
			if(argv[1][1] == 'd')		dN = atoi( &argv[1][2] );
			else if(argv[1][1] == 'f')	numN = atoi( &argv[1][2] );
			else				N = atoi( &argv[1][1] );

			break;
		case 'g':
			if(argv[1][1] == 'd')		dg = atof( &argv[1][2] )*J/2.;
			else if(argv[1][1] == 'f')	numg = atoi( &argv[1][2] );
			else				g = atof( &argv[1][1] )*J/2.;

			break;
		case 'h':
			if(argv[1][1] == 'd')		dh = atof( &argv[1][2] )*J/2.;
			else if(argv[1][1] == 'f')	numh = atoi( &argv[1][2] );
			else						h = atof( &argv[1][1] )*J/2.;
			break;
		case '-':
			//	std::sprintf(folder, "%s", &argv[1][1]);
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}

out_file << "# The Magnetization along the x-axis(longitudinal) in the site = Lsize / 2 \n# N		g		h		Mx		Energy\n";

for(int indN=0;indN<numN;++indN)
for(int indg=0;indg<numg;++indg)
for(int indh=0;indh<numh;++indh)
{
    N += dN;
    g += dg;
    h += dh;

    // Make N spin 1/2's
    auto sites = SpinHalf(N, {"ConserveQNs=", false});

    // Make the Hamiltonian for rung-decoupled Heisenberg ladder
    auto ampo = AutoMPO(sites);    
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += -J,"Sx",i,"Sx",i+1;
        ampo += -g, "Sz",i;
        ampo += -h, "Sx",i;
        }
        ampo += -g, "Sz",N;
        ampo += -h, "Sx",N;
        
    auto H = toMPO(ampo);
    //printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    
    
    auto sweeps0 = Sweeps(5); //number of sweeps is 5
    sweeps0.maxdim() = 20,50,200,200,400;
    sweeps0.cutoff() = 1E-12;

    auto psi0 = randomMPS(sites);

    auto [energy0,psi] = dmrg(H,psi0,sweeps0,{"Silent=",true});

    auto psi2 = MPS(psi);
    auto obsite = int(N/2);

    out_file << N  << "		" << g*2./J << "		" << h*2./J << "		" << real(measure(obsite, psi2, sites)) << "		" << energy0 << "\n";
}


    //std::cout << energy0 << "\n";
/*
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;

    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps0,{"Silent=",true,"Weight=",20.0});

    auto DeltaL = en1 - energy0;
    //auto obsite = 3;
    
    auto psi2 = MPS(psi);
    auto obsite = int(N/2);
*/    
    //    std::cout  << "	" << real(measure(obsite, psi2, sites)) << "	" << obsite << "\n";


    out_file.close();
    return 0;
    }
