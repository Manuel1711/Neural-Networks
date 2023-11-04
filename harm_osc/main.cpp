#define ACT_RELU
#include "include/all.h"

int main() {

    //TOT_MAX = 2843;
    //BOOT_MAX = TOT_MAX;

    fs::path folderPath = "fakedata/";
    fs::path stringInput = "bootstrap/bootstrap1k_secondofile.dat";//"data/corr.txt";

    TOT_MAX = count_filefolder(folderPath);
std::cout << TOT_MAX << "\n";
    BOOT_MAX = TOT_MAX;

    int num_test(0), time_max, num_samples;

    cx_dmat input_train, out_train;
    cx_dvec input_test, out_test;

    train_and_test_data(time_max, num_samples, input_train, out_train, input_test, out_test, folderPath);

    harm_oscill(num_test, time_max, num_samples, input_test, out_test, stringInput);

    // Define your parameter values
    double sigma_b = 10.0;
    double sigma_w = 3.0;
    double n0 = time_max;
    double eta = 0.1;
    int n_layer = 3;

    char foutput[80];

    std::sprintf(foutput, "data_secntk20.dat");

    FILE * net;
    net = fopen(foutput , "wt");

    std::sprintf(foutput, "test.dat");
    FILE * test;
    test = fopen(foutput, "wt");

    std::cout << "END DATA COLLECTION\n";


    int max_boot_harm = 20;
    cx_dvec mean_out_test = zeros<cx_dvec>(out_test.n_elem), err_out_test = zeros<cx_dvec>(out_test.n_elem);
//std::cout << "sdasda\n";
    #pragma omp parallel for
    for (num_test=0; num_test < max_boot_harm; ++num_test)
        {
        harm_oscill(num_test, time_max, num_samples, input_test, out_test, stringInput);
//std::cout << "sdas5555da\n";
        outputNTK result = distribNTKgp(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);
        //outputNTK result = distribNNGP(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);
        mean_out_test += result.mean;
        //std::cout << result.variance << "\n";
        err_out_test += square(result.mean);
        }

    mean_out_test = mean_out_test / double(max_boot_harm);
    err_out_test = sqrt((err_out_test- double(max_boot_harm)*square(mean_out_test))/ (double(max_boot_harm)-1));

    dvec int_x = linspace(0, 0.3, num_samples);
    for(int i=0; i<num_samples; ++i){
        fprintf(net, "%.9e		%.9e        %.9e\n", int_x(i), real(mean_out_test(i)), real(err_out_test(i)));
        fprintf(test, "%.9e		%.9e\n", int_x(i), real(out_test(i)));
    }


/*
    outputNTK result = distribNTKgp(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);

    //outputNTK result = distribNNGP(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);

    dvec int_x = linspace(0, 0.2, num_samples);
    //dvec int_x = linspace(0.0, 2.5, num_samples);
    //dvec int_x = linspace(0.1, 0.2, num_samples/3);

    for(int i=0; i<num_samples; ++i){
    //for(int i=0; i<num_samples/3; ++i){
        fprintf(net, "%.9e		%.9e\n", int_x(i), real(result.mean(i)));
        fprintf(test, "%.9e		%.9e\n", int_x(i), real(out_test(i)));
    }
*/

    return 0;


}

