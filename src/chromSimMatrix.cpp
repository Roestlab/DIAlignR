#include "chromSimMatrix.h"
#include <Rcpp.h>

double meanVecOfVec(const std::vector<std::vector<double>>& vec){
  double average = 0.0;
  // Sum-up mean of each vector using Range-based for loop.
  // const makes sure we do not accidentally chnage v. auto allows compiler to find type of v. & makes sure we are referening to v instead of making a copy that could cause performance loss.
  // for (auto&& v : vec) average += std::accumulate( v.begin(), v.end(), 0.0)/v.size();
  for (const auto& v : vec) average += std::accumulate( v.begin(), v.end(), 0.0)/v.size();
  return average / vec.size();
}

double eucLenVecOfVec(const std::vector<std::vector<double>>& vec){
  double sos = 0.0; // sum of squares
  for (const auto& v : vec) sos += std::accumulate( v.begin(), v.end(), 0.0, square<double>());
  return std::sqrt(sos);
}

void distToSim(SimMatrix& s, double offset, double Numerator){
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::bind(std::plus<double>(), std::placeholders::_1, offset));
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::bind(std::divides<double>(), Numerator, std::placeholders::_1));
}

std::vector<std::vector<double>> meanNormalizeVecOfVec(const std::vector<std::vector<double>>& d){
  // Calculate overall mean and divide by it.
  double mean_d = meanVecOfVec(d);
  std::vector<std::vector<double>> d_new = divideVecOfVec(d, mean_d);
  return d_new;
}

std::vector<std::vector<double>> L2NormalizeVecOfVec(const std::vector<std::vector<double>>& d){
  // Calculate overall mean and divide by it.
  double eucLen_d = eucLenVecOfVec(d);
  std::vector<std::vector<double>> d_new = divideVecOfVec(d, eucLen_d);
  return d_new;
}

std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& d, double num){
  std::vector<std::vector<double>> result;
  result = d;
  // TODO: Need to understand how does transform work.
  for (auto& v : result) std::transform(v.begin(), v.end(), v.begin(), std::bind(std::divides<double>(), std::placeholders::_1, num));
  return result;
}

void ElemWiseSumOuterProd(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s){
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_col == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");
  int nrow = d1.size();
  int ncol = d2.size();
  for (int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      s.data[i*ncol + j] += d1[i]*d2[j]; // Summing outer product of vectors across fragment-ions.
    }
  }
}

void ElemWiseSumOuterEucl(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s){
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_col == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");
  int nrow = d1.size();
  int ncol = d2.size();
  for (int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      s.data[i*ncol + j] += (d1[i]-d2[j]) * (d1[i]-d2[j]); // Summing outer product of vectors across fragment-ions.
    }
  }
}

SimMatrix SumOuterProdMeanNormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  // Mean-normalize each vector of vector.
  std::vector<std::vector<double>> d1_new = meanNormalizeVecOfVec(d1);
  std::vector<std::vector<double>> d2_new = meanNormalizeVecOfVec(d2);
  // Calculate outer dot-product for each fragment-ion and sum element-wise
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  return s;
}

SimMatrix SumOuterProdL2NormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  // L2 normalize each vector of vector.
  std::vector<std::vector<double>> d1_new = L2NormalizeVecOfVec(d1);
  std::vector<std::vector<double>> d2_new = L2NormalizeVecOfVec(d2);
  // Calculate outer dot-product for each fragment-ion and sum element-wise
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  return s;
}

SimMatrix SumOuterEuclMeanNormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  // Mean-normalize each vector of vector.
  std::vector<std::vector<double>> d1_new = meanNormalizeVecOfVec(d1);
  std::vector<std::vector<double>> d2_new = meanNormalizeVecOfVec(d2);
  // Calculate outer-euclidean distance for each sample.
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  // Take sqrt to get eucledian distance from the sum of squared-differences.
  // TODO std::ptr_fun<double, double> wHY?
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::ptr_fun<double, double>(sqrt));
  // Convert distance into similarity.
  distToSim(s, 1.0, 1.0); // similarity = Numerator/(offset + distance)
  return s;
}

SimMatrix SumOuterEuclL2NormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  // L2 normalize each vector of vector.
  std::vector<std::vector<double>> d1_new = L2NormalizeVecOfVec(d1);
  std::vector<std::vector<double>> d2_new = L2NormalizeVecOfVec(d2);
  // Calculate outer-euclidean distance for each sample.
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  // Take sqrt to get eucledian distance from the sum of squared-differences.
  // TODO std::ptr_fun<double, double> wHY?
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::ptr_fun<double, double>(sqrt));
  // Convert distance into similarity.
  distToSim(s, 1.0, 1.0); // similarity = Numerator/(offset + distance)
  return s;
}

SimMatrix SumOuterCosineL2NormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  // L2 normalize each vector of vector.
  std::vector<std::vector<double>> d1_new = L2NormalizeVecOfVec(d1);
  std::vector<std::vector<double>> d2_new = L2NormalizeVecOfVec(d2);
  return s;
}

/***
vector<vector<double> > stuff;
stuff.push_back({1,3,2});
stuff.push_back({0,0,0});
stuff.push_back({4,4,4});
***/

