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
  // TODO : n_th element
  // Rcpp::Rcout << "jh  = " << std::setprecision(8) << (double)n*p + m - j << std::endl;
  // std::nth_element(vec.begin(), vec.begin()+1, vec.end(), std::greater<double>());
  sort(vec.begin(), vec.end());
  double sampleQuant = (1.0 - gamma)*vec[j-1] + gamma*vec[j];
  return sampleQuant;
}

} // namespace DIAlign
