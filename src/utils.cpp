#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <functional>

#include "utils.h"

namespace DIAlign
{
namespace Utils
{

#if 1
double getQuantile(std::vector<double> vec, double quantile){
  int n = vec.size();
  double p = quantile;
  double m = 1-p; // Type 7 definition as implemented in R.
  int j = floor(n*p + m);
  double g = n*p + m - j;
  double gamma = g;

  // sort(vec.begin(), vec.end()); This algorithm is O(n^2). Takes ca 10x longer in benchmarks.
  //
  // nth_element does partial sorting and take O(n) in worst case. With having if condition
  // we are utilizing its O(log(n)) performance. With one nth_element, we get substantial error
  // in certain edge cases if results are compared to R output.
  //
  // nth_element is a partial sorting algorithm that rearranges elements in [first, last) such that:
  //
  // The element pointed at by nth is changed to whatever element would occur in that position if [first, last) were sorted.
  // All of the elements before this new nth element are less than or equal to the elements after the new nth element.
  //
  // Note: we are exploiting the fact that after the first call of nth_element,
  // the vector is partially sorted and we only need to look at the elements
  // before the new nth element to find the second largest value.
  double sampleQuant;
  if (p <= 0.5)
  {
    std::nth_element(vec.begin(), vec.begin()+j, vec.end(), std::less<double>());
    // takes 100 ms
    sampleQuant = gamma*vec[j];
    std::nth_element(vec.begin(), vec.begin()+j-1, vec.begin()+j, std::less<double>()); // second largest value
    sampleQuant = sampleQuant + (1.0 - gamma)*vec[j-1];
  }
  else
  {
    int idx = n*(1-p);
    std::nth_element(vec.begin(), vec.begin()+n-j, vec.end(), std::greater<double>());
    sampleQuant = (1.0-gamma)*vec[n-j];

    std::nth_element(vec.begin(), vec.begin()+n-j-1, vec.begin()+n-j, std::greater<double>()); // second largest value
    sampleQuant = sampleQuant + gamma*vec[n-j-1];
  }
  return sampleQuant;
}

#endif
} // namespace Utils
} // namespace DIAlign
