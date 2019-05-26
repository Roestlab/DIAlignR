#ifndef INTERFACE_H
#define INTERFACE_H

#include <Rcpp.h>
#include <vector>
#include "simpleFcn.h"
using namespace Rcpp;

std::vector<std::vector<double> > list2VecOfVec (Rcpp::List l);

void printVecOfVec(Rcpp::List l);

NumericMatrix Vec2NumericMatrix(std::vector<double> vec, int nrow, int ncol);

SimMatrix NumericMatrix2Vec(Rcpp::NumericMatrix mat);

#endif // INTERFACE_H
