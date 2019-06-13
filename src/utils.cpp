#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "utils.h"

namespace DIAlign 
{

double getQuantile(std::vector<double> vec, double quantile){
  int n = vec.size();
  double p = quantile;
  double m = 1-p; // Type 7 definition as implemented in R.
  int j = floor(n*p + m);
  double g = n*p + m - j; // Replacing n*p with 24884.16 outputs correct result.
  double gamma = g;

#if 0
  sort(vec.begin(), vec.end());
  double sampleQuant = (1.0 - gamma)*vec[j-1] + gamma*vec[j];
#else

  // do a simple approximation (saves us from calling nth two times)
  double sampleQuant;
  if (p <= 0.5)
  {
    int idx = n*p;
    std::nth_element(vec.begin(), vec.begin()+idx, vec.end(), std::less<double>());
    sampleQuant = vec[idx];
  }
  else
  {
    int idx = n*(1-p);
    std::nth_element(vec.begin(), vec.begin()+idx, vec.end(), std::greater<double>());
    sampleQuant = vec[idx];
  }
#endif

  return sampleQuant;
}

} // namespace DIAlign
