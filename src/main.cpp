#include <Rcpp.h>
#include "simpleFcn.h"
#include "alignment.h"
#include "affinealignobj.h"
#include "affinealignment.h"
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

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

//' Initialize a S4 object AffineAlignObj1
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param seq1Len (int) Length of sequence1
//' @param seq2Len (int) Length of sequence2
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @export
// [[Rcpp::export]]
S4 setAffineAlignObj_S4(int ROW_SIZE, int COL_SIZE){
  AffineAlignObj obj(ROW_SIZE, COL_SIZE);
  // Creating an object of Person class
  S4 x("AffineAlignObj");
  // Setting values to the slots
  x.slot("M")  = obj.M;
  x.slot("A")  = obj.A;
  x.slot("B")  = obj.B;
  // std::vector<char> traceback;
  // traceback = EnumToChar(obj.Traceback);
  x.slot("Traceback")  = EnumToChar(obj.Traceback);
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  return(x);
}
// setClass("Person", representation(name="char", birth="Date"))


//' Initialize a S4 object AffineAlignObj1
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param seq1Len (int) Length of sequence1
//' @param seq2Len (int) Length of sequence2
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @export
// [[Rcpp::export]]
S4 setAlignObj_S4(int ROW_SIZE, int COL_SIZE){
  AlignObj obj(ROW_SIZE, COL_SIZE);
  // Creating an object of Person class
  S4 x("AlignObj");
  // Setting values to the slots
  x.slot("M")  = obj.M;
  x.slot("Traceback")  = EnumToChar(obj.Traceback);
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  x.slot("indexA_aligned") = obj.indexA_aligned;
  x.slot("indexB_aligned") = obj.indexB_aligned;
  x.slot("score") = obj.score;
  return(x);
}


//' Initialize a S4 object AffineAlignObj1
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param seq1Len (int) Length of sequence1
//' @param seq2Len (int) Length of sequence2
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @export
// [[Rcpp::export]]
S4 doAlignment_S4(NumericMatrix s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment){
  AlignObj obj(signalA_len+1, signalB_len+1);
  obj = doAlignment(s, signalA_len, signalB_len, gap, OverlapAlignment);
  getAlignedIndices(obj);
  // Creating an object of Person class
  S4 x("AlignObj");
  // Setting values to the slots
  x.slot("M")  = obj.M;
  x.slot("Traceback")  = EnumToChar(obj.Traceback);
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  x.slot("indexA_aligned") = obj.indexA_aligned;
  x.slot("indexB_aligned") = obj.indexB_aligned;
  x.slot("score") = obj.score;
  return(x);
}

//' Initialize a S4 object AffineAlignObj
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param seq1Len (int) Length of sequence1
//' @param seq2Len (int) Length of sequence2
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @export
// [[Rcpp::export]]
S4 doAffineAlignment_S4(NumericMatrix s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment){
  AffineAlignObj obj(signalA_len+1, signalB_len+1);
  obj = doAffineAlignment(s, signalA_len, signalB_len, go, ge, OverlapAlignment);
  // printMatrix(obj.M, signalA_len+1, signalB_len+1);
  // getAlignedIndices(obj);
  // Creating an object of Person class
  S4 x("AffineAlignObj");
  // Setting values to the slots
  x.slot("M")  = obj.M;
  x.slot("A")  = obj.A;
  x.slot("B")  = obj.B;
  x.slot("Traceback")  = EnumToChar(obj.Traceback);
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  x.slot("indexA_aligned") = obj.indexA_aligned;
  x.slot("indexB_aligned") = obj.indexB_aligned;
  x.slot("score") = obj.score;
  return(x);
}

//' Initialize a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param seq1Len (int) Length of sequence1
//' @param seq2Len (int) Length of sequence2
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @export
// [[Rcpp::export]]
S4 rcpp_s4(std::string Name){
  // Creating an object of Person class
  S4 x("Person");
  // Setting values to the slots
  x.slot("name")  = Name;
  x.slot("birth") = Date("1889-12-21");
  return(x);
}

// Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// s <- matrix(NA, nrow = 4, ncol = 5)
// s[, 1] <- c(-2, 10, -2, -2)
// s[, 2] <- c(-2, -2, 10, -2)
// s[, 3] <- c(10, -2, -2, -2)
// s[, 4] <- c(-2, -2, -2, 10)
// s[, 5] <- c(10, -2, -2, -2)
// doAlignment_S4(s, 4, 5, 22, F)
