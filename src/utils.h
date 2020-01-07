#ifndef UTILS_H
#define UTILS_H

#include <vector>
//#include <cstring>
// iostream is required for std::string why?
#include <iostream>

#ifdef DIALIGN_PURE_CPP
#else
#define DIALIGN_USE_Rcpp
#endif

#ifdef DIALIGN_USE_DIALIGN_PRECONDITION
#define DIALIGN_PRECONDITION(condition, message) assert(condition); // If you don't put the message, C++ will output the code.
#else
#define DIALIGN_PRECONDITION(condition, message); // If DIALIGN_USE_DIALIGN_PRECONDITION is defined, compiler will replace calls with empty.
#endif

namespace DIAlign
{
namespace Utils
{
  double getQuantile(std::vector<double> vec, double quantile);
} // namespace Utils
} // namespace DIAlign

#endif // UTILS_H
