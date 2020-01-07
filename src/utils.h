#ifndef UTILS_H
#define UTILS_H

#include <vector>
//#include <cstring>
// iostream is required for std::string why?
#include <iostream>

#ifdef PURE_CPP
#else
#define USE_Rcpp
#endif

#ifdef USE_PRECONDITION
#define PRECONDITION(condition, message) assert(condition); // If you don't put the message, C++ will output the code.
#else
#define PRECONDITION(condition, message); // If USE_PRECONDITION is defined, compiler will replace calls with empty.
#endif

namespace DIAlign
{
namespace Utils
{
  double getQuantile(std::vector<double> vec, double quantile);
} // namespace Utils
} // namespace DIAlign

#endif // UTILS_H
