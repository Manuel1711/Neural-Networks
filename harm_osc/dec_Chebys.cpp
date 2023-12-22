#include <iostream>
#include <armadillo>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/chebyshev.hpp>

using namespace arma;
using namespace boost::math;


double discr_integration(const int n_cheb, const dvec axis, const dvec func){

	int dimen = size(func).n_rows;
	int dim_axis = size(func).n_rows;

	if(dimen != dim_axis)
		std::cerr << "Mismatch of axis and function\n";

	// Multiplication between chebyshev and function
	double a = axis(0), b = axis(dimen-1);

	auto function = [n_cheb,a,b,dimen,axis,func](double xvar) {
		for(int ind=1; ind < dimen; ++ind)
			if(xvar >= axis(ind-1) && xvar < axis(ind)){
				double yvar = 2./(b-a)*(xvar - a/2. - b/2.);
				double yyfunc = (func(ind)-func(ind-1))/(axis(ind) - axis(ind-1))*(xvar-axis(ind)) + func(ind);
				return chebyshev_t(n_cheb, yvar) / sqrt(1.0 - std::pow(yvar,2.)) * 2./(b-a)*yyfunc;
			};
	};

	// Perform the numerical integration
	//integrator.integrate(function, a, b, integral, );
	double result = quadrature::gauss_kronrod<long double,61>::integrate(function,a,b,15,1e-16);

/*
	// Computation of the integral
	double result = 0.0;

    for (int i = 1; i < dimen-2; ++i) {
        double h = axis(i+1) - axis(i);
        result += 0.5 * (yfunc(i-1) + yfunc(i)) * h;
    }

    //std::cout << result << "\n";
*/
    return result;

}

int main(void) {
    int N = 50; // Grado massimo della scomposizione
    double a(0), b(1); // Estremi dell'intervallo
    double integral(0.), error(0.);
	const int max_iterations = 15;

    dvec coefficients(N + 1, fill::zeros); // Vettore dei coefficienti

	dvec xx(10), yy(10);

	for(int ii=0; ii<10; ++ii){
		xx(ii) = ii + 1.;
		yy(ii) = xx(ii)*xx(ii);
	}
	a = xx(0); b = xx(9);

    for (int n = 0; n <= N; ++n) {

		// Choose an appropriate quadrature rule (e.g., 15-point Gauss-Kronrod)

		quadrature::gauss_kronrod<double, max_iterations> integrator;

/*
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
*/


///*
		integral = discr_integration(n, xx, yy);
		coefficients(n) = (2.0 / M_PI) * integral;
		if(n==0) coefficients(n) *= 0.5;
//*/

	}

    // Ora hai calcolato i coefficienti a_n

    // Puoi utilizzare questi coefficienti per stimare f(x) utilizzando la scomposizione
	std::cout.precision(9);
	std::cout << coefficients << std::endl;

    return 0;
}
