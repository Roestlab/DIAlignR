#ifndef SIMPLEFCN_H
#define SIMPLEFCN_H

#include <vector>
#include "affinealignobj.h"
#include "similarityMatrix.h"
#include "utils.h"

#ifdef DIALIGN_USE_Rcpp
#include <Rcpp.h>
#endif // DIALIGN_USE_RCPP

namespace DIAlign
{

#ifdef DIALIGN_USE_Rcpp
using namespace Rcpp;

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

#endif // DIALIGN_USE_RCPP

/// Returns similatrity score matrix between seq1 and seq2.
SimMatrix getseqSim(std::string seq1, std::string seq2, double match, double misMatch);

} // namespace DIAlign

#endif // SIMPLEFCN_H
