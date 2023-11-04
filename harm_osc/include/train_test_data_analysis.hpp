#ifndef __TRAIN_TEST_DATA
#define __TRAIN_TEST_DATA

#include "global_var.h"

void train_and_test_data(int& time_max, int& num_samples, cx_dmat& input_train, cx_dmat& out_train, cx_dvec& input_test, cx_dvec& out_test, const fs::path folderPath){

    time_max = 100;
    num_samples = 10000;
    const int Ntotdata = time_max + num_samples;

    cx_dmat input(TOT_MAX, time_max);
    cx_dmat output(TOT_MAX, num_samples);

    //fs::path folderPath = "fakedata";

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
//std::cout << entry.path() << '\n';
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

    int num_train = BOOT_MAX-1;
    int num_test = BOOT_MAX-1; // BOOT_MAX - 1

    input_train.resize(num_train, time_max);
    input_test.resize(time_max);

    out_train.resize(num_train, num_samples);
    //out_train.resize(num_train, num_samples/3);
    out_test.resize(num_samples);

    input_train = input.submat(0, 0, num_train, time_max-1);
    input_test = (input.row(num_test)).t();

    out_train = output.submat(0, 0, num_train, num_samples-1);
    //out_train = output.submat(0, num_samples/3, num_train, num_samples/3*2-1);
    out_test = (output.row(num_test)).t();

}

#endif
