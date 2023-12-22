#include <iostream>
#include <fstream>
#include <armadillo>

#include <boost/math/special_functions/chebyshev.hpp>

using namespace arma;
using namespace boost::math;

double a(1), b(10); // Estremi dell'intervallo

double approximatedFunction(const dvec& coefficients, double xvar, int N) {
    double result = 0.0;

    for (int n = 0; n <= N; ++n) {
        double yvar = 2./(b-a)*(xvar - a/2. - b/2.);
        result += coefficients(n) * chebyshev_t(n, yvar);
    }

    return result;
}

int main() {
    int N = 50; // Grado massimo della scomposizione
    double x; // Punto in cui calcolare f(x)

    dvec coefficients(N + 1); // Coefficienti a_n

    std::ifstream input_file("temp.dat");

	double value;
	int iiv(0);
	while( input_file >> value){
		coefficients(iiv) = value;
		++ iiv;
	}

	for(int i=1; i < 11; ++i){
        x = double(i);
        double result = approximatedFunction(coefficients, x, N);

        std::cout << "f(" << x << ") = " << result << std::endl;
    }
	input_file.close();
    return 0;
}
