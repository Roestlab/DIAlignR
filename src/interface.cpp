#include "interface.h"

namespace DIAlign
{
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
  for(unsigned int j = 0; j < VecOfVec.size(); j++){
    for (const auto& i : VecOfVec[j]) Rcpp::Rcout<< i << " ";
    Rcpp::Rcout<< std::endl;
  }
  // CharacterVector x = l.names();
  // Rcpp::Rcout<< x[0]<< std::endl; // Error if NULL
  // Rcpp::Rcout << typeid(l).name() << std::endl; // N4Rcpp6VectorILi19ENS_15PreserveStorageEEE
}

SimMatrix NumericMatrix2Vec(Rcpp::NumericMatrix mat){
  SimMatrix s;
  s.n_row = mat.nrow();
  s.n_col = mat.ncol();
  // Cast to std vector
  mat = transpose(mat);
  s.data = Rcpp::as<std::vector<double> >(mat);
  return s;
}
} // namespace DIAlign
