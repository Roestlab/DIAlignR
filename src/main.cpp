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
//' @return s (matrix) Numeric similarity matrix. Rows and columns expresses seq1 and seq2, respectively
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' getSeqSimMat(seq1, seq2, Match, MisMatch)
//' @export
// [[Rcpp::export]]
NumericMatrix getSeqSimMat(std::string seq1, std::string seq2, float Match, float MisMatch){
  int ROW_SIZE = seq1.size();
  int COL_SIZE = seq2.size();
  NumericMatrix s(ROW_SIZE, COL_SIZE);
  for(int j = 0; j < COL_SIZE; j++){
    for(int i = 0; i < ROW_SIZE; i++){
      seq1[i] == seq2[j] ?  s(i, j) = Match : s(i, j) = MisMatch;
    }
  }
  return(s);
}

//' Get a dummy S4 object of C++ class AffineAlignObj
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return affineAlignObj (S4class) A S4class dummy object from C++ AffineAlignObj struct
//' @examples
//' x <- setAffineAlignObj_S4(4, 5)
//' x@signalA_len # 3
//' @export
// [[Rcpp::export]]
S4 setAffineAlignObj_S4(int ROW_SIZE, int COL_SIZE){
  AffineAlignObj obj(ROW_SIZE, COL_SIZE); // Initializing C++ AffineAlignObj struct
  S4 x("AffineAlignObj"); // Creating an empty S4 object of AffineAlignObj class
  // Setting values to the slots
  x.slot("M")  = obj.M;
  x.slot("A")  = obj.A;
  x.slot("B")  = obj.B;
  x.slot("Traceback")  = EnumToChar(obj.Traceback); // EnumToChar adds 48 to get ASCII character value of single-digit numeral.
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  return(x);
}

//' Get a dummy S4 object of C++ class AlignObj
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return AlignObj (S4class) A S4class dummy object from C++ AlignObj struct
//' @examples
//' x <- setAlignObj_S4(4, 5)
//' x@signalA_len # 3
//' @export
// [[Rcpp::export]]
S4 setAlignObj_S4(int ROW_SIZE, int COL_SIZE){
  AlignObj obj(ROW_SIZE, COL_SIZE); // Initializing C++ AlignObj struct
  S4 x("AlignObj"); // Creating an empty S4 object of AlignObj class
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

//' Perform non-affine global and overlap alignment on a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param s (NumericMatrix) A numeric matrix with similarity values of two sequences or signals
//' @param signalA_len (int) Length of signalA or sequenceA. Expresses along the rows of s
//' @param signalB_len (int) Length of signalB or sequenceB. Expresses along the columns of s
//' @param gap (float) Penalty for introducing gaps in alignment
//' @param OverlapAlignment (bool) An inpul for gap-free alignment. False: Global alignment, True: Gap-free overlap alignment
//' @return AlignObj (S4class) An object from C++ class of AlignObj
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' s <- getSeqSimMat(seq1, seq2, Match, MisMatch)
//' obj_Global <- doAlignment_S4(s, 4, 5, 22, FALSE)
//' obj_Global@score # -2 -4 -6 4 -18
//' obj_Olap <- doAlignment_S4(s, 4, 5, 22, TRUE)
//' obj_Olap@score # 0 10 20 18 18 18
//' @export
// [[Rcpp::export]]
S4 doAlignment_S4(NumericMatrix s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment){
  AlignObj obj(signalA_len+1, signalB_len+1); // Initializing C++ AlignObj struct
  obj = doAlignment(s, signalA_len, signalB_len, gap, OverlapAlignment); // Performs alignment on s matrix and returns AlignObj struct
  getAlignedIndices(obj); // Performs traceback and fills aligned indices in AlignObj struct
  S4 x("AlignObj"); // Creating an empty S4 object of AlignObj class
  // Copying values to slots
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
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' s <- getSeqSimMat(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignment_S4(s, 4, 5, 22, 7, FALSE)
//' objAffine_Global@score # -2 -4 -6 4 -33
//' objAffine_Olap <- doAffineAlignment_S4(s, 4, 5, 22, 7, TRUE)
//' objAffine_Olap@score # 0 10 20 18 18 18
//' @export
// [[Rcpp::export]]
S4 doAffineAlignment_S4(NumericMatrix s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment){
  AffineAlignObj obj(signalA_len+1, signalB_len+1);
  obj = doAffineAlignment(s, signalA_len, signalB_len, go, ge, OverlapAlignment);
  getAffineAlignedIndices(obj);
  // printMatrix(obj.M, signalA_len+1, signalB_len+1);
  // getAlignedIndices(obj);
  S4 x("AffineAlignObj");  // Creating an empty S4 object of AffineAlignObj class
  // Copying values to slots
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
//' @return s (S4class) A similarity matrix
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
// Match=10; MisMatch=-2; go=22; ge=7; gap=go
// seq1 = "GCAT"; seq2 = "CAGTG"
// s <- getSeqSimMat(seq1, seq2, Match, MisMatch)
// s <- matrix(NA, nrow = 4, ncol = 5)
// s[, 1] <- c(-2, 10, -2, -2)
// s[, 2] <- c(-2, -2, 10, -2)
// s[, 3] <- c(10, -2, -2, -2)
// s[, 4] <- c(-2, -2, -2, 10)
// s[, 5] <- c(10, -2, -2, -2)
// doAlignment_S4(s, 4, 5, 22, F)
// library(Biostrings)
// mat <- nucleotideSubstitutionMatrix(match = Match, mismatch = MisMatch, baseOnly = TRUE)
// pairwiseAlignment(seq1, subject = seq2, type = "global", substitutionMatrix = mat, gapOpening = 0, gapExtension = 22)
// pairwiseAlignment(seq1, subject = seq2, type = "overlap", substitutionMatrix = mat, gapOpening = 0, gapExtension = 22)
