#ifndef UTILS_H
#define UTILS_H

#include <vector>
//#include <cstring>
// iostream is required for std::string why?
#include <iostream>

#ifdef USE_PRECONDITION
#define PRECONDITION(condition, message) assert(condition); // If you don't put the message, C++ will output the code.
#else
#define PRECONDITION(condition, message); // If USE_PRECONDITION is defined, compiler will replace calls with empty.
#endif

// #define USE_Rcpp //TODO: Why moving from simpleFcn.h to utils.h solves multiple main() definition problem when build in R?

double getQuantile(std::vector<double> vec, double quantile);

#endif // UTILS_H
