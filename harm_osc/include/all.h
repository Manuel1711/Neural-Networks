#ifndef __ALL_FUNC
#define __ALL_FUNC

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <glob.h>
#include <omp.h>
#include <armadillo>

#include <filesystem>



namespace fs = std::filesystem;

using namespace arma;
typedef cx_double complex;

#include "global_var.h"
#include "count_filefolder.hpp"

#include "distributions.hpp" 
#include "harm_osc.hpp"
#include "recursive_functions.hpp" 
#include "train_test_data_analysis.hpp"

#endif

