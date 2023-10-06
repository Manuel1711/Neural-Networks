#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>

// Structure to represent a data point
struct DataPoint {
    int time;
    double value;
};

int main(int argc, char *argv[]) {
	char in_string[50], out_string[50];
	while( argc > 1 ) {

            switch(argv[1][0]) {
            case '-':
				if(argv[1][1] == 'i') std::sprintf(in_string, "%s", &argv[1][2]);
				else if(argv[1][1] == 'o') std::sprintf(out_string, "%s", &argv[1][2]);
                break;
        	}
        	++argv;
        	--argc;
	}

    // Define the number of bootstraps
    int numBootstraps = 1000;

    // Read the dataset from the file into memory
    // std::cout << in_string << std::endl;
    std::ifstream inputFile(in_string);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    std::vector<DataPoint> dataset;
    int time;
    double value;
    while (inputFile >> time >> value) {
        dataset.push_back({time, value});
    }
    inputFile.close();

    // Create a map to group values by time
    std::map<int, std::vector<double>> groupedData;
    for (const auto& dataPoint : dataset) {
        groupedData[dataPoint.time].push_back(dataPoint.value);
    }

    // Set up a random number generator for sampling with replacement
    std::random_device rd;
    std::mt19937 gen(rd());

    std::ofstream outputFile(out_string);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

    outputFile << "# BootstrapNum \t \t time \t\t Corr" << std::endl;
    // Perform bootstrapping
    for (const auto& pair : groupedData) {
            int currentTime = pair.first;
            const std::vector<double>& values = pair.second;

            for (int bootstrap = 0; bootstrap < numBootstraps; ++bootstrap) {
                std::vector<double> bootstrapSample;

                // Resample values for the current time
                for (int i = 0; i < int(values.size()); ++i) {
                    int randomIndex = std::uniform_int_distribution<int>(0, values.size() - 1)(gen);
                    bootstrapSample.push_back(values[randomIndex]);
                }

                // Calculate the mean of the bootstrap sample
                double mean = 0.0;
                for (double val : bootstrapSample) {
                    mean += val;
                }
                mean /= bootstrapSample.size();

                //Print or store the results for each bootstrap
                outputFile << bootstrap + 1 << " \t \t " << currentTime << " \t \t " << mean << std::endl;
            }
    }

    return 0;
}

