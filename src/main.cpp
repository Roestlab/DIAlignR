#include <Rcpp.h>
using namespace Rcpp;

//' Calculate similarity matrix for two sequences
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
NumericMatrix initializeMatrix(float initVal, int ROW_SIZE, int COL_SIZE){
  NumericMatrix s(ROW_SIZE, COL_SIZE);
  for(int i = 0; i < ROW_SIZE; i++){
    for(int j = 0; j < COL_SIZE; j++){
      s(i, j) = initVal;
    }
  }
  return s;
}
// s <- initializeMatrix(0, 4, 5)

//' Calculate similarity matrix for two sequences
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-05
//' @param seq1 (char) A single string
//' @param seq2 (char) A single string
//' @param Match (float) Score for character match
//' @param MisMatch (float) score for character mismatch
//' @param s (matrix) A similarity matrix
//' @return s (matrix) Updated similarity matrix
//' @export
// [[Rcpp::export]]
void getseqSimMat(std::string seq1, std::string seq2, float Match, float MisMatch, NumericMatrix s){
  int ROW_SIZE = seq1.size();
  int COL_SIZE = seq2.size();
  for(int j = 0; j < COL_SIZE; j++){
    for(int i = 0; i < ROW_SIZE; i++){
      seq1[i] == seq2[j] ?  s(i, j) = Match : s(i, j) = MisMatch;
    }
  }
}
// Match=10; MisMatch=-2; go=22; ge=7; gap=go
// seq1 = "GCAT"; seq2 = "CAGTG"
// getseqSimMat(seq1, seq2, Match, MisMatch, s)
