#include "simpleFcn.h"

// Outputs a NumericMatrix of given row and column size.
NumericMatrix initializeMatrix(float initVal, int ROW_SIZE, int COL_SIZE){
  NumericMatrix s(ROW_SIZE, COL_SIZE);
  for(int i = 0; i < ROW_SIZE; i++){
    for(int j = 0; j < COL_SIZE; j++){
      s(i, j) = initVal;
    }
  }
  return s;
}

double getQuantile(std::vector<double> vec, double quantile){
  int n = vec.size();
  double p = quantile;
  double m = 1-p; // Type 7 definition as implemented in R.
  int j = floor(n*p + m);
  double g = n*p + m - j; // Replacing n*p with 24884.16 outputs correct result.
  double gamma = g;
  // TODO : Precision and n_th element
  // long double jh = 24884.16;
  // Rcpp::Rcout << "g  = " << g << std::endl;
  //Rcpp::Rcout << "jh  = " << std::setprecision(8) << (double)n*(double)p + m - j << std::endl;
  //std::nth_element(vec.begin(), vec.begin()+1, vec.end(), std::greater<double>());
  sort(vec.begin(), vec.end());
  double sampleQuant = (1.0 - gamma)*vec[j] + gamma*vec[j+1];
  return sampleQuant;
}
