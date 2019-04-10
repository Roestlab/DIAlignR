#ifndef SIMPLEFCN_H
#define SIMPLEFCN_H

#include <Rcpp.h>
#include <vector>
#include "affinealignobj.h"
using namespace Rcpp;

//' Initialize a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-05
//' @param initVal (char) Matrix intialization value
//' @param ROW_SIZE (int) Number of rows
//' @param COL_SIZE (int) Number of columns
//' @return s (matrix) a similarity matrix
//' @export
// [[Rcpp::export]]
NumericMatrix initializeMatrix(float initVal, int ROW_SIZE, int COL_SIZE);
// s <- initializeMatrix(0, 4, 5)

// Template definitions should always be in header file.
template<class T>
void printMatrix(T Mat, int ROW_SIZE, int COL_SIZE){
  for(int i = 0; i < ROW_SIZE; i++){
    for(int j = 0; j < COL_SIZE; j++){
      Rcpp::Rcout << Mat[i*COL_SIZE + j] << " ";
    }
    Rcpp::Rcout << std::endl;
  }
}

#endif // SIMPLEFCN_H
