#ifndef UTILS_H
#define UTILS_H

#include <vector>

#ifdef USE_PRECONDITION
#define PRECONDITION(condition, message) assert(condition); // If you don't put the message, C++ will output the code.
#else
#define PRECONDITION(condition, message) ; // If you don't put the message, C++ will output the code.
#endif

namespace DIAlign 
{
  double getQuantile(std::vector<double> vec, double quantile);
} // namespace DIAlign

#endif // UTILS_H
