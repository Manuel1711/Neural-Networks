#include <iostream>
#include <armadillo>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/chebyshev.hpp>

using namespace arma;
using namespace boost::math;


int main() {
    int N = 30; // Grado massimo della scomposizione
    double a(0), b(1); // Estremi dell'intervallo
    double integral(0.), error(0.);
	const int max_iterations = 15;

    dvec coefficients(N + 1, fill::zeros); // Vettore dei coefficienti

    for (int n = 0; n <= N; ++n) {

		// Choose an appropriate quadrature rule (e.g., 15-point Gauss-Kronrod)

		quadrature::gauss_kronrod<double, max_iterations> integrator;

		auto function = [n,a,b](double xvar) {
			double yvar = 2./(b-a)*(xvar - a/2. - b/2.);
			return chebyshev_t(n, yvar) / sqrt(1.0 - std::pow(yvar,2.)) * 2./(b-a) *(xvar*xvar);
		};

		// Perform the numerical integration
		//integrator.integrate(function, a, b, integral, );
		integral = quadrature::gauss_kronrod<long double,61>::integrate(function,a,b,15,1e-16);
		//std::cout << integral << std::endl;
        coefficients(n) = (2.0 / M_PI) * integral;
		if(n==0) coefficients(n) *= 0.5;
    }

    // Ora hai calcolato i coefficienti a_n

    // Puoi utilizzare questi coefficienti per stimare f(x) utilizzando la scomposizione
	std::cout.precision(9);
	std::cout << coefficients << std::endl;

    return 0;
}
