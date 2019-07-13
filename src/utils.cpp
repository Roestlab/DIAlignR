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
  double g = n*p + m - j;
  double gamma = g;

  // sort(vec.begin(), vec.end()); This algorithm is O(n^2).
  // nth_element does partial sorting and take O(n) in worst case. With having if condition
  // we are utilizing its O(log(n)) performance. With one nth_element, we get substantial error
  // in certain edge cases if results are compared to R output.
  double sampleQuant;
  if (p <= 0.5)
  {
    std::nth_element(vec.begin(), vec.begin()+j, vec.end(), std::less<double>());
    sampleQuant = gamma*vec[j];
    std::nth_element(vec.begin(), vec.begin()+j-1, vec.end(), std::less<double>());
    sampleQuant = sampleQuant + (1.0 - gamma)*vec[j-1];
  }
  else
  {
    int idx = n*(1-p);
    std::nth_element(vec.begin(), vec.begin()+n-j-1, vec.end(), std::greater<double>());
    sampleQuant = gamma*vec[n-j-1];
    std::nth_element(vec.begin(), vec.begin()+n-j, vec.end(), std::greater<double>());
    sampleQuant = sampleQuant + (1.0-gamma)*vec[n-j];
  }
  return sampleQuant;
}

} // namespace DIAlign
