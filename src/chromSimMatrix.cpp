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

std::vector<double> perSampleEucLenVecOfVec(const std::vector<std::vector<double>>& vec){
  std::vector<double> mag;
  mag.resize(vec[0].size(), 0.0);
  int n_frag = vec.size();
  for (int i = 0; i < mag.size(); i++){
    for (int fragIon = 0; fragIon < n_frag; fragIon++){
      mag[i] += vec[fragIon][i] * vec[fragIon][i];
    }
    mag[i] = std::sqrt(mag[i]);
  }
  return mag;
}

std::vector<double> perSampleSqrSumVecOfVec(const std::vector<std::vector<double>>& vec){
  std::vector<double> mag;
  mag.resize(vec[0].size(), 0.0);
  int n_frag = vec.size();
  for (int i = 0; i < mag.size(); i++){
    for (int fragIon = 0; fragIon < n_frag; fragIon++){
      mag[i] += vec[fragIon][i] * vec[fragIon][i];
    }
  }
  return mag;
}

std::vector<double> perSampleMeanVecOfVec(const std::vector<std::vector<double>>& vec){
  std::vector<double> mean;
  mean.resize(vec[0].size(), 0.0);
  int n_frag = vec.size();
  for (int i = 0; i < mean.size(); i++){
    for (int fragIon = 0; fragIon < n_frag; fragIon++){
      mean[i] += vec[fragIon][i];
    }
    mean[i] = mean[i]/n_frag;
  }
  return mean;
}

std::vector<double> perSampleSumVecOfVec(const std::vector<std::vector<double>>& vec){
  std::vector<double> sum;
  sum.resize(vec[0].size(), 0.0);
  int n_frag = vec.size();
  for (int i = 0; i < sum.size(); i++){
    for (int fragIon = 0; fragIon < n_frag; fragIon++){
      sum[i] += vec[fragIon][i];
    }
  }
  return sum;
}

void distToSim(SimMatrix& s, double offset, double Numerator){
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::bind(std::plus<double>(), std::placeholders::_1, offset));
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::bind(std::divides<double>(), Numerator, std::placeholders::_1));
}

void clamp(std::vector<double>& vec, double minValue, double maxValue){
  for (auto& i : vec) {i = (i > maxValue) ? maxValue : i;
    i = (i < minValue) ? minValue : i;}
}

double getQuantile(std::vector<double> vec, float quantile){
  int n = vec.size();
  float p = quantile;
  float m = 1-p; // Type 7 definition as implemented in R.
  int j = floor(n*p + m);
  float g = n*p + m - j;
  float gamma = g;
  int idx = floor((1.0-quantile)*vec.size()); // make it ceil
  std::nth_element(vec.begin(), vec.begin()+1, vec.end(), std::greater<double>());
  float sampleQuant = (1.0 - gamma)*vec[n-j] + gamma*vec[n-j-1];
  return sampleQuant;
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

void ElemWiseSumOuterProdMeanSub(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s, const std::vector<double>& mean1, const std::vector<double>& mean2){
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_col == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");
  int nrow = d1.size();
  int ncol = d2.size();
  for (int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      s.data[i*ncol + j] += (d1[i]-mean1[i])*(d2[j]-mean2[j]); // Summing outer product of vectors across fragment-ions.
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

void ElemWiseOuterCosine(const std::vector<double>& d1, const std::vector<double>& d2, const std::vector<double>& d1_mag, const std::vector<double>& d2_mag, SimMatrix& s){
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_col == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");
  PRECONDITION(d1_mag.size() == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(d2_mag.size() == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");
  int nrow = s.n_row;
  int ncol = s.n_col;
  for (int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      s.data[i*ncol + j] += d1[i]*d2[j]/(d1_mag[i]*d2_mag[j]); // Summing outer product of vectors across fragment-ions.
    }
  }
}

void SumOuterProd(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension (row)");
  PRECONDITION(s.n_col == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension (column)");
  std::vector<std::vector<double>> d1_new, d2_new;
  if(Normalization == "mean"){
    // Mean normalize each vector of vector.
    d1_new = meanNormalizeVecOfVec(d1);
    d2_new = meanNormalizeVecOfVec(d2);
  } else if(Normalization == "L2"){
    // L2 normalize each vector of vector.
    d1_new = L2NormalizeVecOfVec(d1);
    d2_new = L2NormalizeVecOfVec(d2);
  } else {
    d1_new = d1;
    d2_new = d2;
  }
  // Calculate outer dot-product for each fragment-ion and sum element-wise
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
}

void SumOuterCov(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension (row)");
  PRECONDITION(s.n_col == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension (column)");
  std::vector<std::vector<double>> d1_new, d2_new;
  if(Normalization == "mean"){
    // Mean normalize each vector of vector.
    d1_new = meanNormalizeVecOfVec(d1);
    d2_new = meanNormalizeVecOfVec(d2);
  } else if(Normalization == "L2"){
    // L2 normalize each vector of vector.
    d1_new = L2NormalizeVecOfVec(d1);
    d2_new = L2NormalizeVecOfVec(d2);
  } else {
    d1_new = d1;
    d2_new = d2;
  }
  std::vector<double> d1_mean = perSampleMeanVecOfVec(d1_new);
  std::vector<double> d2_mean = perSampleMeanVecOfVec(d2_new);
  // Calculate outer dot-product for each fragment-ion and sum element-wise
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProdMeanSub(d1_new[fragIon], d2_new[fragIon], s, d1_mean, d2_mean);
  }
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::bind(std::divides<double>(), std::placeholders::_1, n_frag));
}

void SumOuterCorr(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension (row)");
  PRECONDITION(s.n_col == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension (column)");
  std::vector<std::vector<double>> d1_new, d2_new;
  if(Normalization == "mean"){
    // Mean normalize each vector of vector.
    d1_new = meanNormalizeVecOfVec(d1);
    d2_new = meanNormalizeVecOfVec(d2);
  } else if(Normalization == "L2"){
    // L2 normalize each vector of vector.
    d1_new = L2NormalizeVecOfVec(d1);
    d2_new = L2NormalizeVecOfVec(d2);
  } else {
    d1_new = d1;
    d2_new = d2;
  }
  std::vector<double> d1_sum = perSampleSumVecOfVec(d1_new);
  std::vector<double> d2_sum = perSampleSumVecOfVec(d2_new);
  std::vector<double> d1_squareSum = perSampleSqrSumVecOfVec(d1_new);
  std::vector<double> d2_squareSum = perSampleSqrSumVecOfVec(d2_new);
  // Calculate outer dot-product for each fragment-ion and sum element-wise
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  double var1, var2 = 0.0;

  for (int i = 0; i < s.n_row; i++){
    for(int j = 0; j < s.n_col; j++){
      var1 = n_frag*d1_squareSum[i]-d1_sum[i]*d1_sum[i];
      var2 = n_frag*d2_squareSum[j]-d2_sum[j]*d2_sum[j];
      if(var1 < 0.0 | var2 <= 0.0){
        Rcpp::Rcout << "In SumOuterCorr the standard deviation is zero" << std::endl;
        s.data[i*s.n_col+j] = -20; // TODO: What to output in this case?
      }
      else
        s.data[i*s.n_col+j] = (n_frag*s.data[i*s.n_col+j] - d1_sum[i]*d2_sum[j])/sqrt(var1*var2); // Summing outer product of vectors across fragment-ions.
    }
  }
}

void SumOuterEucl(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension (row)");
  PRECONDITION(s.n_col == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension (column)");
  std::vector<std::vector<double>> d1_new, d2_new;
  if(Normalization == "mean"){
    // Mean normalize each vector of vector.
    d1_new = meanNormalizeVecOfVec(d1);
    d2_new = meanNormalizeVecOfVec(d2);
  } else if(Normalization == "L2"){
    // L2 normalize each vector of vector.
    d1_new = L2NormalizeVecOfVec(d1);
    d2_new = L2NormalizeVecOfVec(d2);
  } else {
    d1_new = d1;
    d2_new = d2;
  }
  // Calculate outer-euclidean distance for each sample.
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseSumOuterProd(d1_new[fragIon], d2_new[fragIon], s);
  }
  // Take sqrt to get eucledian distance from the sum of squared-differences.
  // TODO std::ptr_fun<double, double> Why? Effectively calls std::pointer_to_unary_function<Arg,Result>(f)
  std::transform(s.data.begin(), s.data.end(), s.data.begin(), std::ptr_fun<double, double>(sqrt));
  // Convert distance into similarity.
  distToSim(s, 1.0, 1.0); // similarity = Numerator/(offset + distance)
}

void SumOuterCosine(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s){
  PRECONDITION(!d1.empty(), "Vector of vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension (row)");
  PRECONDITION(s.n_col == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension (column)");
  std::vector<std::vector<double>> d1_new, d2_new;
  if(Normalization == "mean"){
    // Mean normalize each vector of vector.
    d1_new = meanNormalizeVecOfVec(d1);
    d2_new = meanNormalizeVecOfVec(d2);
  }
  else if(Normalization == "L2"){
    // L2 normalize each vector of vector.
    d1_new = L2NormalizeVecOfVec(d1);
    d2_new = L2NormalizeVecOfVec(d2);
  }
  else {
    d1_new = d1;
    d2_new = d2;
  }
  // No normalization needed for calculating cosine similarity.
  std::vector<double> d1_mag = perSampleEucLenVecOfVec(d1_new);
  std::vector<double> d2_mag = perSampleEucLenVecOfVec(d2_new);
  int n_frag = d1.size();
  for (int fragIon = 0; fragIon < n_frag; fragIon++){
    ElemWiseOuterCosine(d1_new[fragIon], d2_new[fragIon], d1_mag, d2_mag, s);
  }
  clamp(s.data, -1.0, 1.0); // Clamp the cosine similarity between -1.0 and 1.0
}

SimMatrix getSimilarityMatrix(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, \
                              const std::string Normalization, const std::string SimType, double cosAngleThresh, \
                              double dotProdThresh){
  SimMatrix s;
  s.n_row = d1[0].size();
  s.n_col = d2[0].size();
  s.data.resize(s.n_row*s.n_col, 0.0);
  if (SimType == "dotProductMasked"){
    SumOuterProd(d1, d2, Normalization, s);
    SimMatrix s2;
    s2.n_row = s.n_row;
    s2.n_col = s.n_col;
    s2.data.resize(s.n_row*s.n_col, 0.0);
    SumOuterCosine(d1, d2, Normalization, s2);
    for(auto& i : s2.data) i = std::cos(2*std::acos(i));
    double Quant = getQuantile(s.data, dotProdThresh);
    // Rcpp::Rcout << Quant << std::endl;
    std::vector<double> MASK;
    MASK.resize(s.n_row*s.n_col, 0.0);
    for(int i = 0; i < MASK.size(); i++){
      MASK[i] = (s.data[i] < Quant) ? 0.0 : 1.0;
      MASK[i] = (MASK[i]*s2.data[i] + (1.0-MASK[i]) > cosAngleThresh) ? 1.0 : 0.0;
      s.data[i] = s.data[i] * MASK[i];
    }
  }
  else if (SimType == "dotProduct")
    SumOuterProd(d1, d2, Normalization, s);
  else if(SimType == "cosineAngle")
    SumOuterCosine(d1, d2, Normalization, s);
  else if(SimType == "cosine2Angle"){
    SumOuterCosine(d1, d2, Normalization, s);
    for(auto& i : s.data) i = std::cos(2*std::acos(i));
    clamp(s.data, -1.0, 1.0); // Clamp the cosine similarity between -1.0 and 1.0
  }
  else if(SimType == "euclidianDist")
    SumOuterEucl(d1, d2, Normalization, s);
  else if(SimType == "covariance")
    SumOuterCov(d1, d2, Normalization, s);
  else if(SimType == "correlation")
    SumOuterCorr(d1, d2, Normalization, s);
  else
    Rcpp::Rcout << "getChromSimMat should have value from given choices only!" << std::endl;
  return s;
}

/***
vector<vector<double> > stuff;
stuff.push_back({1,3,2});
stuff.push_back({0,0,0});
stuff.push_back({4,4,4});
***/

