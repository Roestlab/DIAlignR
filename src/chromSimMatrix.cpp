#include "chromSimMatrix.h"
#include <Rcpp.h>

double meanVecOfVec(std::vector<std::vector<double>> vec){
  double average = 0.0;
  // Sum-up mean of each vector using Range-based for loop.
  // const makes sure we do not accidentally chnage v. auto allows compiler to find type of v. & makes sure we are referening to v instead of making a copy that could cause performance loss.
  // for (auto&& v : vec) average += std::accumulate( v.begin(), v.end(), 0.0)/v.size();
  for (const auto& v : vec) average += std::accumulate( v.begin(), v.end(), 0.0)/v.size();
  return average / vec.size();
}

std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& d, double num){
  std::vector<std::vector<double>> result;
  result = d;
  // TODO: Need to understand how does transform work.
  for (auto& v : result) std::transform(v.begin(), v.end(), v.begin(), std::bind(std::divides<double>(), std::placeholders::_1, num));
  return result;
}

SimMatrix OuterProdMeanNormAllFunc(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  SimMatrix s;
  double mean_d1 = meanVecOfVec(d1);
  double mean_d2 = meanVecOfVec(d2);
  std::vector<std::vector<double>> d1_new = divideVecOfVec(d1, mean_d1);
  std::vector<std::vector<double>> d2_new = divideVecOfVec(d2, mean_d2);
  return s;
}

/***
vector<vector<double> > stuff;
stuff.push_back({1,3,2});
stuff.push_back({0,0,0});
stuff.push_back({4,4,4});
***/

