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

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace DIAlign
{
/** Utility namespace */
namespace Utils
{
  /**
   * @brief Returns the quantile of the vector
   *
   * @param vec The vector with the values
   * @param quantile The n-th quantile to compute
   *
  */
  double getQuantile(std::vector<double> vec, double quantile);
} // namespace Utils
} // namespace DIAlign

#endif // UTILS_H
