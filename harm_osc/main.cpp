#include "include/all.h"

int main() {

    TOT_MAX = 1000;
    BOOT_MAX = 50; // TOT_MAX

    int time_max, num_samples;

    cx_dmat input_train, out_train;
    cx_dvec input_test, out_test;

    train_and_test_data(time_max, num_samples, input_train, out_train, input_test, out_test);

    harm_oscill(time_max, num_samples, input_test, out_test);

    // Define your parameter values
    double sigma_b = 10.0;
    double sigma_w = 3.0;
    double n0 = time_max;
    double eta = 0.1;
    int n_layer = 3;

    char foutput[80];

    std::sprintf(foutput, "data1.dat");

    FILE * net;
    net = fopen(foutput , "wt");

    std::sprintf(foutput, "test1.dat");
    FILE * test;
    test = fopen(foutput, "wt");

    outputNTK result = distribNTKgp(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);
    //outputNTK result = distribNNGP(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);

    dvec int_x = linspace(0, 0.3, num_samples);
    for(int i=0; i<num_samples; ++i){
        fprintf(net, "%.9e		%.9e\n", int_x(i), real(result.mean(i)));
        fprintf(test, "%.9e		%.9e\n", int_x(i), real(out_test(i)));
    }

    return 0;


}

