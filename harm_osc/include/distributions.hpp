#ifndef __DISTRIB
#define __DISTRIB

#include "recursive_functions.hpp"

struct outputNTK {
  cx_dvec mean;
  complex variance;
};


outputNTK distribNNGP(const cx_dvec x, const cx_dmat input_train, const cx_dmat output_train, const double eta, const int n_layer, const double sigma_w, const double sigma_b, const double n0){

    outputNTK result;

    int tdim = int(input_train.n_rows);
    complex Kxx = Kappa(x, x, n_layer+1, sigma_w, sigma_b, n0);
    cx_dvec KxCalX(tdim);
    cx_dmat KCalXCalX(tdim, tdim);

    for(int index=0; index < tdim; ++index){
        KxCalX(index) = Kappa(x, (input_train.row(index)).t(), n_layer+1, sigma_w, sigma_b, n0);
        for(int jndex=index; jndex < tdim; ++jndex){
            KCalXCalX(index, jndex) = Kappa((input_train.row(index)).t(), (input_train.row(jndex)).t(), n_layer+1, sigma_w, sigma_b, n0);
            KCalXCalX(jndex, index) = KCalXCalX(index, jndex);
        }
    }

    cx_dmat K_inverse = inv(KCalXCalX);
    //Aux = expm(-eta*KCalXCalX*time)
    //Aux1 = np.eye(tdim) - np.matmul(Aux,Aux)
    //Aux = np.eye(tdim) - Aux
    //print(Aux, output_train)
    result.mean = (KxCalX.t() * (K_inverse * output_train)).t();

    cx_dvec temp_var = Kxx - KxCalX.t()*(K_inverse*KxCalX);
    result.variance = temp_var(0);

    return result;
}

outputNTK distribNTKgp(
    const cx_dvec x,
    const cx_dmat input_train,
    const cx_dmat output_train,
    const double eta, const int n_layer, const double sigma_w, const double sigma_b, const double n0
) {
    outputNTK result;
    int tdim = int(input_train.n_rows);

    complex Kxx = Kappa(x, x, n_layer+1, sigma_w, sigma_b, n0);

    cx_dmat KCalXCalX(tdim, tdim);
    cx_dmat TCalXCalX(tdim, tdim);
    cx_dvec KxCalX(tdim);
    cx_dvec TxCalX(tdim);

    for (int index = 0; index < tdim; ++index) {
        KxCalX(index) = Kappa(x, (input_train.row(index)).t(), n_layer + 1, sigma_w, sigma_b, n0);
        TxCalX(index) = Theta(x, (input_train.row(index)).t(), n_layer + 1, sigma_w, sigma_b, n0);

        for (int jndex = index; jndex < tdim; ++jndex) {
            KCalXCalX(index, jndex) = Kappa((input_train.row(index)).t(), (input_train.row(jndex)).t(), n_layer + 1, sigma_w, sigma_b, n0);
            KCalXCalX(jndex, index) = KCalXCalX(index, jndex);
            TCalXCalX(index, jndex) = Theta((input_train.row(index)).t(), (input_train.row(jndex)).t(), n_layer + 1, sigma_w, sigma_b, n0);
            TCalXCalX(jndex, index) = TCalXCalX(index, jndex);
        }
    }

    cx_dmat T_inverse = inv(TCalXCalX); //inv(TCalXCalX/2. + TCalXCalX.st()/2.);
    cx_dmat Aux = eye<cx_dmat>(tdim, tdim); // Identity matrix

    result.mean = (TxCalX.t() * (T_inverse * (Aux * output_train))).t();

    cx_dvec temp = (Kxx + TxCalX.t() * (T_inverse * (Aux * KCalXCalX * (Aux * (T_inverse * TxCalX))))) ;

    cx_dvec var_aux = TxCalX.t() * (T_inverse * (Aux * KxCalX));
    result.variance = temp(0) - var_aux(0) - conj(var_aux(0));

    return result;
}

#endif
