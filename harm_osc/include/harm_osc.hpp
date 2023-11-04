#ifndef __HARM_OSC
#define __HARM_OSC

#include "global_var.h"

void harm_oscill(const int num_test, int& time_max, int& num_samples, cx_dvec& input_test, cx_dvec& out_test, const fs::path stringInput){
    dvec boot, time, corr;

    //std::ifstream inputFile("bootstrap/bootstrap1k_secondofile.dat");
    std::ifstream inputFile(stringInput);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the input file." << std::endl;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        // Check if the line starts with '#', indicating a comment
        if (line.empty() || line[0] == '#') {
            continue; // Skip comments
        }

        std::istringstream iss(line);

        double val1, val2, val3;

        while (iss >> val1 >> val2 >> val3) {
            boot.resize(boot.n_elem + 1);
            time.resize(time.n_elem + 1);
            corr.resize(corr.n_elem + 1);

            boot(boot.n_elem - 1) = val1;
            time(time.n_elem - 1) = val2;
            corr(corr.n_elem - 1) = val3;
        }
    }
    inputFile.close();

    time_max = int(max(time) + 1);
    int harm_boot_max = int(max(boot));

    //int num_train = harm_boot_max - 1;

/*******IMPORTANT:choice_of_boot_data************/
    //int num_test = 0;//harm_boot_max - 1;
/*******IMPORTANT:choice_of_boot_data************/

    cx_dmat input(harm_boot_max, time_max);
    for(int i=0; i<harm_boot_max; ++i){
        int auxind = i*time_max;
        for(int jj=0; jj<time_max; ++jj)
            input(i, jj) = corr(jj + auxind);
    }

    input_test = (input.row(num_test)).t();
//std::cout << input_test << '\n';
    //num_samples = 10000;
    double etaharm(0.05), smear_sigma(0.01);
    dvec int_x = linspace(0, 0.2, num_samples);
    //dvec int_x = linspace(0.0, 2.5, num_samples);
    //dvec int_x = linspace(0.1, 0.2, num_samples/3);

    double aux_zeta = 0.5 + 0.5*erf(etaharm/sqrt(2)/smear_sigma);

    //dvec pdf = std::pow(2., -1./2.)*(1.0 / (aux_zeta * smear_sigma * sqrt(2.0 * datum::pi))) * exp(-0.5 * square((int_x - etaharm) / smear_sigma));

    dvec pdf = std::pow(2., -3./2.)*(1.0 / (aux_zeta * smear_sigma * sqrt(2.0 * datum::pi))) * exp(-0.5 * square((int_x - etaharm) / smear_sigma));

    pdf = pdf + std::pow(2., -3./2.)/9.*(1.0 / (aux_zeta * smear_sigma * sqrt(2.0 * datum::pi))) * exp(-0.5 * square((int_x - 3*etaharm) / smear_sigma));

    for(int ii=0; ii<num_samples; ++ii)
    //for(int ii=0; ii<num_samples/3; ++ii)
        out_test(ii) = complex(pdf(ii), 0.);
}


#endif
