#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>

 


int T=200;
namespace bq=boost::math::quadrature;

long double theta(long double x) {
    return (x >= 0) ? 1.0 : 0.0;
}

long double gaussiana(long double x, long double xbar, long double sigma) {
  long double E = exp(-pow(x-xbar,2)/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma);
  long double Z= (1+erf(xbar/(sqrt(2)*sigma)))/2;
  return E;
}


long double x(long double E){
  return 1-2*exp(-E);
}

long double rho(long double omega, int N, const std::vector<long double>& omegabar, long double sigma, const std::vector<long double>& cn_train, long double omegaZero) {
  long double rho_sum = 0;
  for (int i=0; i<N; i++){
    rho_sum += cn_train[i]*(boost::math::chebyshev_t(i, x(omega)) - boost::math::chebyshev_t(i, x(omegaZero)));
    //rho_sum += gaussiana(omega,omegabar[i],sigma);
  }
  //std::cout << "AAA " << theta(omega-omegaZero)*rho_sum << std::endl;
  return rho_sum;
}

long double integral(int tau, const std::vector<long double>& omegabar, long double s, int N, const std::vector<long double>& cn_train, long double omegaZero){
  
  const long double infLimit=0;
  const long double supLimit=std::numeric_limits<long double>::infinity();
  const auto f_cs=
    [=](const long double& omega) -> long double
    {
      return rho(omega,N,omegabar,s, cn_train, omegaZero)*(exp(-omega*tau)+exp(-(T-tau)*omega));
    };
  return bq::gauss_kronrod<long double,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);

}

long double rho(const std::vector<long double>& cn, const std::vector<long double>& En, long double S, long double E){
  
  if (cn.size() != En.size()) {
    std::cerr << "Errore: cn e En devono avere lo stesso numero di elementi." << std::endl;
    exit(0); 
  }
  
  long double Sum = 0;
  for (size_t i = 0; i < cn.size(); i++){
    //std::cout << cn[i] << std::endl;
    Sum += cn[i]*gaussiana(E,En[i],S);
  }
  return Sum;
}

std::vector<long double> cn_fact(int N, std::vector<long double>& cn_pre){

  for(int i=0; i<N; i++){

    if(i==0){}
    else{
      cn_pre[i]=cn_pre[i]/(pow(i,1+1e-7));
    }
    //std::cout << "cn:  " <<  cn_pre[i] << std::endl; 
  }

  return cn_pre;

}

std::vector<long double> riempiVettoreCasuale(int dimensione, long double min, long double max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<long double> distribuzione(min, max);
  
  std::vector<long double> vettore(dimensione);
  for (int i = 0; i < dimensione; i++) {
    vettore[i] = distribuzione(gen);
  }
  
  return vettore;
}


int main(int argc, char *argv[]){
  
  int Npeaks = 10001;
  std::vector<long double> cn(Npeaks);
  std::vector<long double> En(Npeaks);
  
  cn[0]=1;
  En[0]=0.8;
  
  long double S=0.44;
  /*
  std::vector<long double> valoriCasuali_cn = riempiVettoreCasuale(Npeaks - 1, 0, 0.01);
  std::vector<long double> valoriCasuali_En = riempiVettoreCasuale(Npeaks - 1, 1.2, 50);

  for (int i = 1; i < Npeaks; i++) {
    cn[i] = valoriCasuali_cn[i - 1];
    En[i] = valoriCasuali_En[i - 1];
  }

  
  
  std::ofstream outputFile("spectrum_real1.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Errore nell'apertura del file di output." << std::endl;
    exit(0);
  }
  
  outputFile << "# cn     En" << std::endl;
  for(int i=0; i<Npeaks; i++)
    outputFile << cn[i] << "     " << En[i] << std::endl;
  
  outputFile.close();
  
  std::ofstream outputFile1("spectrum_smeared1.txt");
  if (!outputFile1.is_open()) {
    std::cerr << "Errore nell'apertura del file di output." << std::endl;
    exit(0);
  }
  
  
  int Nbins=10000;
  std::vector<long double> Ebin = riempiVettoreCasuale(Nbins, 0, 50);
  std::sort(Ebin.begin(), Ebin.end());
  
  outputFile1 << "# Ebin    rho(Ebin)" << std::endl;
  for(int i=0; i<Nbins; i++)
    outputFile1 << Ebin[i] << "  " << rho(cn,En,S,Ebin[i]) << std::endl;

  outputFile1.close();
  */
  int N=6;
  int Nfiles = 10000;
  int Nbins=10000;
  std::vector<double> OmegaBins(Nbins);
  double inizio = 0.0;
  double fine = 2.5;
  double passo = (fine - inizio) / (Nbins - 1);
  
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
				N = atoi( &argv[1][1] );
			break;
        case 'f':
				Nfiles = atoi( &argv[1][1] );
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}

  for (int i = 0; i < Nbins; i++) {
    OmegaBins[i] = inizio + i * passo;
  }
  
  for(int s=0; s<Nfiles; s++){
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<long double> distribuzione(0.2, 1.3);
    long double E0 = distribuzione(gen);
    
    std::vector<long double> cn_pre = riempiVettoreCasuale(N, -1, 1);
    std::vector<long double> cn_train = cn_fact(N, cn_pre);
    char Fake_path[1024];
    sprintf(Fake_path,"fakedata1/sp_corr_N%d_%d.txt",N,s+1000);
    std::ofstream outputFile2(Fake_path);
    if (!outputFile2.is_open()) {
      std::cerr << "Errore nell'apertura del file di output." << std::endl;
      exit(0);
    }
    
    //outputFile2 << "C(t)" << std::endl;
    
    int t[int(T/2)];
    std::vector<long double> omegabar = riempiVettoreCasuale(N, inizio, fine);
    for(int i=0; i<int(T/2); i++){
      long double Ct_tot=0; 
      /*for(int k=0; k<Npeaks; k++)
	Ct_tot += cn[k]*(exp(-En[k]*i) + exp(-En[k]*(T-i)));*/
      Ct_tot = integral(i,omegabar,S,N,cn_train,E0);
      outputFile2  << Ct_tot << " \n";
    }
    
    //outputFile2 << std::endl << std::endl;
    
    
    
    for(int i=0; i<Nbins; i++)
      outputFile2 << rho(OmegaBins[i],N,omegabar,S, cn_train, E0) << " \n";
    
    outputFile2.close();

    /*
    if(s==0){
      
      std::ofstream outputFile3("FakeData/rho_train_last.txt");
      if (!outputFile3.is_open()) {
	std::cerr << "Errore nell'apertura del file di output." << std::endl;
	exit(0);
      }
      
      outputFile3 << "# Omegabin   rho(omega) " << std::endl;
      for(int i=0; i<Nbins; i++)
	outputFile3 << OmegaBins[i] << "  " << rho(OmegaBins[i],N,omegabar,S,cn_train,E0) << std::endl;
      
      outputFile3.close();
      
      }
    */

    
  }//Nfiles
  
  
  
  //system("zip -xvf FakeData.zip FakeData/*");
  
  return 0;
  
  
    
}
