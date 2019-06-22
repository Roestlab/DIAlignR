#ifndef INTERFACE_H
#define INTERFACE_H

#include <Rcpp.h>
#include <vector>
#include "simpleFcn.h"
using namespace Rcpp;

std::vector<std::vector<double> > list2VecOfVec (Rcpp::List l);

void printVecOfVec(Rcpp::List l);

template<class T>
NumericMatrix Vec2NumericMatrix(std::vector<T> vec, int nrow, int ncol){
  NumericMatrix mat(ncol, nrow, vec.begin());
  // STL vector has matrix stored as rowwise, but, NumericMatrix are filled columnwise.
  // Therefore, we first switch rows and column and then take transposition.
  mat = transpose(mat);
  return mat;
}

SimMatrix NumericMatrix2Vec(Rcpp::NumericMatrix mat);

#endif // INTERFACE_H
