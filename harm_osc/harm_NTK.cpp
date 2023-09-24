#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <glob.h>

#include <armadillo>

#include <filesystem>
namespace fs = std::filesystem;

using namespace arma;
typedef cx_double complex;

int TOT_MAX = 2100;
int BOOT_MAX = 2100; // TOT_MAX



inline complex mathT(const complex a, const complex b, const complex c, const complex d) {
    complex a1 = abs(a);
    complex d1 = abs(d);
    complex fact_or = b / sqrt(a1 * d1);
    complex output = sin(acos(fact_or)) + (datum::pi - acos(fact_or)) * fact_or;
    return 1.0 / 2.0 / datum::pi * sqrt(a1 * d1) * output;
}

inline complex mathTp(const complex a, const complex b, const complex c, const complex d) {
    complex a1 = abs(a);
    complex d1 = abs(d);
    return 1.0 / 2.0 / datum::pi * (datum::pi - acos(b / sqrt(a1 * d1)));
}


// Function to calculate Kappa
complex Kappa(const cx_dvec& arr1, const cx_dvec& arr2, const int l_index, const double sigma_w, const double sigma_b, const double n0) {
    if (l_index == 1) {
        return dot(arr1, arr2) / n0 * sigma_w * sigma_w + sigma_b * sigma_b;
    } else {
        complex T11 = Kappa(arr1, arr1, l_index - 1, sigma_w, sigma_b, n0);
        complex T12 = Kappa(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
        complex T21 = T12;
        complex T22 = Kappa(arr2, arr2, l_index - 1, sigma_w, sigma_b, n0);
        return sigma_w * sigma_w * mathT(T11, T12, T21, T22) + sigma_b * sigma_b;
    }
}

// Function to calculate Theta
complex Theta(const cx_dvec& arr1, const cx_dvec& arr2, const int l_index, const double sigma_w, const double sigma_b, const double n0) {
    if (l_index == 1) {
        return dot(arr1, arr2) / n0 * sigma_w * sigma_w + sigma_b * sigma_b;
    } else {
        complex T11 = Kappa(arr1, arr1, l_index - 1, sigma_w, sigma_b, n0);
        complex T12 = Kappa(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
        complex T21 = T12;
        complex T22 = Kappa(arr2, arr2, l_index - 1, sigma_w, sigma_b, n0);
        return sigma_w * sigma_w * mathTp(T11, T12, T21, T22) * Theta(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
    }
}

struct outputNTK {
  cx_dvec mean;
  complex variance;
};


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


void train_and_test_data(int& time_max, int& num_samples, cx_dmat& input_train, cx_dmat& out_train, cx_dvec& input_test, cx_dvec& out_test){

    time_max = 200;
    num_samples = 100;
    const int Ntotdata = 300;

    cx_dmat input(TOT_MAX, time_max);
    cx_dmat output(TOT_MAX, num_samples);

    fs::path folderPath = "fakedata";

    // Check if the folder exists
    if (!fs::is_directory(folderPath)) {
        std::cerr << "The folder does not exist." << std::endl;
    }

    int index(0);
    // Iterate through all files in the folder
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        // Check if the file is a regular text file
        //if (entry.is_regular_file() && entry.path().extension() == ".txt") {
        if (entry.is_regular_file()) {

            std::ifstream file(entry.path());

            // Check if the file can be opened
            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << entry.path() << std::endl;
                continue;
            }

            double number;
            int ind_read(0.);
            cx_dvec xx(Ntotdata);
            while (file >> number) {
               xx(ind_read) = number;
               ++ind_read;
            }

            input.row(index) = (xx.subvec(0, time_max-1)).t();
            output.row(index) = (xx.subvec(Ntotdata-num_samples, Ntotdata-1)).t();

            // Close the file
            file.close();
            ++index;
        }
    }

    int num_train = BOOT_MAX - 2;
    int num_test = BOOT_MAX - 1; // BOOT_MAX - 1

    input_train.resize(num_train, time_max);
    input_test.resize(time_max);

    out_train.resize(num_train, num_samples);
    out_test.resize(num_samples);

    input_train = input.submat(0, 0, num_train, time_max-1);
    input_test = (input.row(num_test)).t();

    out_train = output.submat(0, 0, num_train, num_samples-1);
    out_test = (output.row(num_test)).t();

}

void harm_oscill(int& time_max, int& num_samples, cx_dvec& input_test, cx_dvec& out_test){
    dvec boot, time, corr;

    std::ifstream inputFile("bootstrap/bootstrap_primofile.dat");

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
    int num_test = 1;//harm_boot_max - 1;
/*******IMPORTANT:choice_of_boot_data************/

    cx_dmat input(harm_boot_max, time_max);
    for(int i=0; i<harm_boot_max; ++i){
        int auxind = i*time_max;
        for(int jj=0; jj<time_max; ++jj)
            input(i, jj) = corr(jj + auxind);
    }

    input_test = (input.row(num_test)).t();
//std::cout << input_test << '\n';
    num_samples = 100;
    double etaharm(0.05), smear_sigma(0.1);
    dvec int_x = linspace(0, 20, num_samples);

    dvec pdf = (1.0 / (smear_sigma * sqrt(2.0 * datum::pi))) * exp(-0.5 * square((int_x - etaharm) / smear_sigma));

    for(int ii=0; ii<num_samples; ++ii)
        out_test(ii) = complex(pdf(ii), 0.);
}


int main() {

    int time_max, num_samples;

    cx_dmat input_train, out_train;
    cx_dvec input_test, out_test;

    train_and_test_data(time_max, num_samples, input_train, out_train, input_test, out_test);

    harm_oscill(time_max, num_samples, input_test, out_test);

    // Define your parameter values
    double sigma_b = 1.0;
    double sigma_w = 1.0;
    double n0 = time_max;
    double eta = 0.1;
    int n_layer = 3;

    char foutput[80];

	std::sprintf(foutput, "data.dat");

	FILE * net;
	net = fopen(foutput , "wt");

    std::sprintf(foutput, "test.dat");
    FILE * test;
    test = fopen(foutput, "wt");

    outputNTK result = distribNTKgp(input_test, input_train, out_train, eta, n_layer, sigma_w, sigma_b, n0);

    dvec int_x = linspace(0, 20, num_samples);
    for(int i=0; i<num_samples; ++i){
        fprintf(net, "%.9e		%.9e\n", int_x(i), real(result.mean(i)));
        fprintf(test, "%.9e		%.9e\n", int_x(i), real(out_test(i)));
    }

    return 0;


}

