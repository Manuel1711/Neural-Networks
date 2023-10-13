#ifndef __COUNT_FILES
#define __COUNT_FILES

int count_filefolder(const fs::path folderPath){
    // Check if the folder exists
    if (!fs::is_directory(folderPath)) {
        std::cerr << "The folder does not exist." << std::endl;
    }

    int indd(0);
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

            ++indd;
            }
    }
    return indd;
}

#endif
