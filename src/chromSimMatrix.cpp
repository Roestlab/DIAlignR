#include "chromSimMatrix.h"
#include <Rcpp.h>

double meanVecOfVec(std::vector<std::vector<double>> vec){
  double average = 0.0;
  for (const auto& v : vec) average += std::accumulate( v.begin(), v.end(), 0.0)/v.size();
  return average / vec.size();
}

SimMatrix OuterProdMeanNormAllFunc(std::vector<std::vector<double>> d1, std::vector<std::vector<double>> d2){
  SimMatrix s;
  double mean_d1 = meanVecOfVec(d1);
  double mean_d2 = meanVecOfVec(d2);
  return s;
}

/***
vector<vector<double> > stuff;
stuff.push_back({1,3,2});
stuff.push_back({0,0,0});
stuff.push_back({4,4,4});
***/

