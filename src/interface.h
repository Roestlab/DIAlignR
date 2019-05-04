#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

std::vector<std::vector<double> > list2VecOfVec (Rcpp::List l);

void printVecOfVec(Rcpp::List l);

NumericMatrix Vec2NumericMatrix(std::vector<double>, int nrow, int ncol);
