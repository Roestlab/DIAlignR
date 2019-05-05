#include "interface.h"

std::vector<std::vector<double> > list2VecOfVec (Rcpp::List l){
  int len = l.size();
  std::vector<std::vector<double> > VecOfVec;
  for (int i = 0; i < len; i++){
    VecOfVec.push_back(as<std::vector<double>>(l[i]));
  }
  return VecOfVec;
}

void printVecOfVec(Rcpp::List l){
  // Printing output of list2VecOfVec function
  std::vector<std::vector<double> > VecOfVec = list2VecOfVec(l);
  for(int j = 0; j < VecOfVec.size(); j++){
    for (const auto& i : VecOfVec[j]) Rcpp::Rcout<< i << " ";
    Rcpp::Rcout<< std::endl;
  }
  // CharacterVector x = l.names();
  // Rcpp::Rcout<< x[0]<< std::endl; // Error if NULL
  // Rcpp::Rcout << typeid(l).name() << std::endl; // N4Rcpp6VectorILi19ENS_15PreserveStorageEEE
}

NumericMatrix Vec2NumericMatrix(std::vector<double> vec, int nrow, int ncol){
  NumericMatrix mat(nrow, ncol, vec.begin());
  // STL vector has matrix stored as rowwise, but, NumericMatrix are filled columnwise.
  // Hence, matrix transposition is needed for similar representation.
  mat = transpose(mat);
  return mat;
}
