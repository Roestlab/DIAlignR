#include <Rcpp.h>
#include "simpleFcn.h"
#include "interface.h"
#include "chromSimMatrix.h"
#include "alignment.h"
#include "gapPenalty.h"
#include "affinealignobj.h"
#include "affinealignment.h"
#include "constrainMat.h"
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
//' @param Match (double) Score for character match
//' @param MisMatch (double) score for character mismatch
//' @return s (matrix) Numeric similarity matrix. Rows and columns expresses seq1 and seq2, respectively
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' getSeqSimMat(seq1, seq2, Match, MisMatch)
//' @export
// [[Rcpp::export]]
NumericMatrix getSeqSimMat(std::string seq1, std::string seq2, double Match, double MisMatch){
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

//' Calculate similarity matrix of two fragment-ion chromatogram group.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-05
//' @param l1 (list) A list of vectors. Length should be same as of l2.
//' @param l2 (list) A list of vectors. Length should be same as of l1.
//' @param Normalization (char) A character string. Normalization must be selected from (L2, mean).
//' @param SimType (char) A character string. Similarity type must be selected from (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclidianDist, covariance, correlation).
//' @return s (matrix) Numeric similarity matrix. Rows and columns expresses seq1 and seq2, respectively
//' @examples
//' # Get similarity matrix of dummy chromatograms
//' r1 <- list(c(1.0,3.0,2.0,4.0), c(0.0,0.0,0.0,1.0), c(4.0,4.0,4.0,5.0))
//' r2 <- list(c(1.4,2.0,1.5,4.0), c(0.0,0.5,0.0,0.0), c(2.0,3.0,4.0,0.9))
//' getChromSimMat(r1, r2, "L2", "dotProductMasked")
//' matrix(c(0.1251213, 0.1623915, 0.1437564, 0.2076481, 0.1863509, 0.2395940,
//' 0.2129724, 0.3128033, 0.2329386, 0.2728709, 0.2529048, 0.3460802, 0.1011619,
//' 0.2076481, 0.1544050, 0.2728709), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "dotProduct")
//' matrix(c(0.1251213, 0.1623915, 0.1437564, 0.2076481, 0.1863509,
//'  0.2395940, 0.2129724, 0.3128033, 0.2329386, 0.2728709, 0.2529048,
//'   0.3460802, 0.1011619, 0.2076481, 0.1544050, 0.2728709), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "cosineAngle")
//' matrix(c(0.9338568, 0.9994629, 0.9892035, 0.9859998, 0.9328152, 0.9889961,
//'  0.9828722, 0.9961742, 0.9935327, 0.9597374, 0.9945055, 0.9391117, 0.4495782,
//'  0.7609756, 0.6326436, 0.7715167), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "cosine2Angle")
//' matrix(c(0.7441769, 0.9978523, 0.9570470, 0.9443912, 0.7402886, 0.9562264, 0.9320755,
//' 0.9847260, 0.9742143, 0.8421918, 0.9780822, 0.7638617, -0.5957588, 0.1581678,
//' -0.1995241, 0.1904762), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "euclidianDist")
//' matrix(c(0.7387025, 0.7127694, 0.7250831, 0.6869622, 0.6984783, 0.6713737,
//' 0.6842335, 0.6413183, 0.6744739, 0.6568703, 0.6653819, 0.6296096, 0.7586910,
//' 0.6869622, 0.7179039, 0.6568703), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "covariance")
//' matrix(c(0.016564523, 0.018930883, 0.017747703, 0.018930883, 0.021445141, 0.022924117,
//' 0.022184629, 0.022924117, 0.036974382, 0.034016431, 0.035495406, 0.034016431,
//' -0.002514258, 0.018487191, 0.007986466, 0.018487191), 4, 4, byrow = F)
//'
//' getChromSimMat(r1, r2, "L2", "correlation")
//' matrix(c(0.87372107, 0.99853837, 0.97435470, 0.99853837, 0.92261291, 0.98624138, 0.99339927,
//' 0.98624138, 0.99053606, 0.91129318, 0.98974332, 0.91129318, -0.06486282, 0.47693252,
//' 0.21444787, 0.47693252), 4, 4, byrow = F)
//' @export
// [[Rcpp::export]]
NumericMatrix getChromSimMat(Rcpp::List l1, Rcpp::List l2, std::string Normalization, std::string SimType, double cosAngleThresh = 0.3, double dotProdThresh = 0.96){
  std::vector<std::vector<double> > r1 = list2VecOfVec(l1);
  std::vector<std::vector<double> > r2 = list2VecOfVec(l2);
  // printVecOfVec(l1);
  SimMatrix s = getSimilarityMatrix(r1, r2, Normalization, SimType, cosAngleThresh, dotProdThresh);
  // printMatrix(s.data, s.n_row, s.n_col);
  NumericMatrix simMat = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
  return simMat;
}

//' Get a mask for constraining similarity matrix.
//'
//' This function takes in timeVectors from both runs, a global-fit object and
//' sample-length of window of no constraining. Outside of window, all elements
//' of matrix are either equally weighted or weighted proportional to distance
//' from window-boundry.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return affineAlignObj (S4class) A S4class dummy object from C++ AffineAlignObj struct
//' @examples
//' tA <- c(3353.2, 3356.6, 3360.0, 3363.5)
//' tB <- c(3325.9, 3329.3, 3332.7, 3336.1)
//' B1p <- 3325.751; B2p <- 3336.119
//' noBeef <- 1
//' m <- getGlobalAlignMask(tA, tB, B1p, B2p, noBeef, FALSE)
//' matrix(c(0.0000, 0.0000, 0.7071, 1.4142, 0.0000, 0.0000, 0.0000, 0.7071, 0.7071, 0.0000,
//' 0.0000, 0.0000, 1.4142, 0.7071, 0.0000, 0.0000), 4, 4, byrow = F)
//' @export
// [[Rcpp::export]]
NumericMatrix getGlobalAlignMask(const std::vector<double>& tA, const std::vector<double>& tB, double B1p, double B2p, int noBeef = 50, bool hardConstrain = false){
  SimMatrix MASK;
  MASK.n_row = tA.size();
  MASK.n_col = tB.size();
  MASK.data.resize(MASK.n_row*MASK.n_col, 0.0);
  double A1 = tA[0], A2 = tA[MASK.n_row-1];
  double B1 = tB[0], B2 = tB[MASK.n_col-1];
  calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, hardConstrain);
  NumericMatrix mask = Vec2NumericMatrix(MASK.data, MASK.n_row, MASK.n_col);
  return mask;
}

//' Constrain similarity matrix using mask and constrainValue.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return affineAlignObj (S4class) A S4class dummy object from C++ AffineAlignObj struct
//' @examples
//'
//' @export
// [[Rcpp::export]]
NumericMatrix constrainSimMain(const NumericMatrix& sim, const NumericMatrix& MASK, double samples4gradient = 100.0){
  SimMatrix s = NumericMatrix2Vec(sim);
  SimMatrix mask = NumericMatrix2Vec(MASK);
  auto maxIt = max_element(std::begin(s.data), std::end(s.data));
  double maxVal = *maxIt;
  constrainSimilarity(s, mask, -2.0*maxVal/samples4gradient);
  // sim = Vec2NumericMatrix(s.data, s.n_row, s.n_col); // This code doesn't update sim matrix. Why?
  NumericMatrix s1 = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
  return s1;
}

//' Calculates gap penalty for dynamic programming based alignment.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return affineAlignObj (S4class) A S4class dummy object from C++ AffineAlignObj struct
//' @examples
//'
//' @export
// [[Rcpp::export]]
double getGapPenaltyMain(const NumericMatrix& sim, std::string SimType, double gapQuantile = 0.5){
  SimMatrix s = NumericMatrix2Vec(sim);
  double gapPenalty = getGapPenalty(s, gapQuantile, SimType);
  return gapPenalty;
}

//' Outputs aligned chromatograms.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param ROW_SIZE (int) Number of rows of a matrix
//' @param COL_SIZE (int) Number of columns of a matrix
//' @return affineAlignObj (S4class) A S4class dummy object from C++ AffineAlignObj struct
//' @examples
//' simMeasure <- "dotProductMasked"
//' run_pair <- c("run1", "run2")
//' peptide <- peptides[1]
//' r1 <- lapply(StrepChroms[[run_pair[1]]][[peptide]], `[[`, 2)
//' r2 <- lapply(StrepChroms[[run_pair[2]]][[peptide]], `[[`, 2)
//' tRunAVec <- StrepChroms[[run_pair[1]]][[peptide]][[1]][["time"]]
//' tRunBVec <- StrepChroms[[run_pair[2]]][[peptide]][[1]][["time"]]
//' noBeef <- 6
//' B1p <- predict(Loess.fit, tRunAVec[1]); B2p <- predict(Loess.fit, tRunAVec[length(tRunAVec)])
//' Alignobj <- alignChromatograms_cpp(r1, r2, "hybrid", tRunAVec, tRunBVec, "mean", simMeasure, B1p, B2p, noBeef)
//'
//' @export
// [[Rcpp::export]]
S4 alignChromatograms_cpp(Rcpp::List l1, Rcpp::List l2, std::string alignType,
                            const std::vector<double>& tA, const std::vector<double>& tB,
                            std::string Normalization, std::string SimType,
                            double B1p = 0.0, double B2p =0.0, int noBeef = 0,
                            double goFactor = 0.125, double geFactor = 40,
                            double cosAngleThresh = 0.3, bool OverlapAlignment = true,
                            double dotProdThresh = 0.96, double gapQuantile = 0.5,
                            bool hardConstrain = false, double samples4gradient = 100.0){
  std::vector<std::vector<double> > r1 = list2VecOfVec(l1);
  std::vector<std::vector<double> > r2 = list2VecOfVec(l2);
  SimMatrix s = getSimilarityMatrix(r1, r2, Normalization, SimType, cosAngleThresh, dotProdThresh);
  double gapPenalty = getGapPenalty(s, gapQuantile, SimType);
  if (alignType == "hybrid"){
    SimMatrix MASK;
    MASK.n_row = tA.size();
    MASK.n_col = tB.size();
    MASK.data.resize(MASK.n_row*MASK.n_col, 0.0);
    double A1 = tA[0], A2 = tA[MASK.n_row-1];
    double B1 = tB[0], B2 = tB[MASK.n_col-1];
    calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, hardConstrain);
    auto maxIt = max_element(std::begin(s.data), std::end(s.data));
    double maxVal = *maxIt;
    constrainSimilarity(s, MASK, -2.0*maxVal/samples4gradient);
  }
  AffineAlignObj obj(s.n_row+1, s.n_col+1); // Initializing C++ AffineAlignObj struct
  obj = doAffineAlignment(s, s.n_row, s.n_col, gapPenalty*goFactor, gapPenalty*geFactor, OverlapAlignment); // Performs alignment on s matrix and returns AffineAlignObj struct
  getAffineAlignedIndices(obj); // Performs traceback and fills aligned indices in AffineAlignObj struct
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
//' @param sim (NumericMatrix) A numeric matrix with similarity values of two sequences or signals
//' @param signalA_len (int) Length of signalA or sequenceA. Expresses along the rows of s
//' @param signalB_len (int) Length of signalB or sequenceB. Expresses along the columns of s
//' @param gap (double) Penalty for introducing gaps in alignment
//' @param OverlapAlignment (bool) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment
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
S4 doAlignment_S4(NumericMatrix sim, int signalA_len, int signalB_len, double gap, bool OverlapAlignment){
  AlignObj obj(signalA_len+1, signalB_len+1); // Initializing C++ AlignObj struct
  SimMatrix s = NumericMatrix2Vec(sim);
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

//' Perform affine global and overlap alignment on a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param sim (NumericMatrix) A numeric matrix with similarity values of two sequences or signals
//' @param signalA_len (int) Length of signalA or sequenceA. Expresses along the rows of s
//' @param signalB_len (int) Length of signalB or sequenceB. Expresses along the columns of s
//' @param go (double) Penalty for introducing first gap in alignment
//' @param ge (double) Penalty for introducing subsequent gaps in alignment
//' @param OverlapAlignment (bool) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' s <- getSeqSimMat(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignment_S4(s, 4, 5, 22, 7, FALSE)
//' objAffine_Global@score # -2  -4  -6  4 -18
//' objAffine_Olap <- doAffineAlignment_S4(s, 4, 5, 22, 7, TRUE)
//' objAffine_Olap@score # 0 10 20 18 18 18
//'
//' seq1 = "CAT"; seq2 = "CAGTG"
//' s <- getSeqSimMat(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignment_S4(s, 3, 5, 22, 7, FALSE)
//' objAffine_Global@score # 10  20  -2  -9 -11
//' objAffine_Olap <- doAffineAlignment_S4(s, 3, 5, 22, 7, TRUE)
//' objAffine_Olap@score # 10 20 18 18 18
//' @export
// [[Rcpp::export]]
S4 doAffineAlignment_S4(NumericMatrix sim, int signalA_len, int signalB_len, double go, double ge, bool OverlapAlignment){
  AffineAlignObj obj(signalA_len+1, signalB_len+1); // Initializing C++ AffineAlignObj struct
  SimMatrix s = NumericMatrix2Vec(sim);
  // printMatrix(s.data, s.n_row, s.n_col);
  obj = doAffineAlignment(s, signalA_len, signalB_len, go, ge, OverlapAlignment);  // Performs alignment on s matrix and returns AffineAlignObj struct
  getAffineAlignedIndices(obj); // Performs traceback and fills aligned indices in AffineAlignObj struct
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

// Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// Match=10; MisMatch=-2; go=22; ge=7; gap=go
// seq1 = "GCAT"; seq2 = "CAGTG"
// library(Biostrings)
// mat <- nucleotideSubstitutionMatrix(match = Match, mismatch = MisMatch, baseOnly = TRUE)
// pairwiseAlignment(seq1, subject = seq2, type = "global", substitutionMatrix = mat, gapOpening = 0, gapExtension = 22)
// pairwiseAlignment(seq1, subject = seq2, type = "overlap", substitutionMatrix = mat, gapOpening = 0, gapExtension = 22)
// pairwiseAlignment(seq1, subject = seq2, type = "global", substitutionMatrix = mat, gapOpening = 15, gapExtension = 7)
// pairwiseAlignment(seq1, subject = seq2, type = "overlap", substitutionMatrix = mat, gapOpening = 15, gapExtension = 7)

/***
MeanNormA <- sapply(r1, function(x) sum(x)/4)
MeanNormA <- mean(MeanNormA)
MeanNormB <- sapply(r2, function(x) sum(x)/4)
MeanNormB <- mean(MeanNormB)
L2NormA <- sapply(r1, function(x) x)
L2NormA <- sqrt(rowSums(L2NormA^2))
L2NormB <- sapply(r2, function(x) x)
L2NormB <- sqrt(rowSums(L2NormB^2))
outerProdList <- list()
for (i in 1:3){
  NormIntensityA <- r1[[i]]/L2NormA
  NormIntensityB <- r2[[i]]/L2NormB
  outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
  }
add <- function(x) Reduce("+", x)
add(outerProdList)
s1 <- getChromSimMat(r1, r2, "L2", "dotProduct")
s2 <- getChromSimMat(r1, r2, "L2", "cosine2Angle")
MASK <- (s1 > quantile(s1, 0.96))
AngleGreat <- (((1*MASK)*s2) + (1-MASK)) > 0.3
s <- s1*(1*AngleGreat)
***/

// Comparing the alignment with
/***
gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
simMeasure <- "dotProductMasked"
run_pair <- c("run1", "run2")
Err <- matrix(NA, nrow = length(peptides), ncol = 1)
rownames(Err) <- peptides
  for(peptide in peptides){
    r1 <- lapply(StrepChroms[[run_pair[1]]][[peptide]], `[[`, 2)
    r2 <- lapply(StrepChroms[[run_pair[2]]][[peptide]], `[[`, 2)
    tRunAVec <- StrepChroms[[run_pair[1]]][[peptide]][[1]][["time"]]
    tRunBVec <- StrepChroms[[run_pair[2]]][[peptide]][[1]][["time"]]
    noBeef <- ceiling(RSEdistFactor*min(globalStrep["RSE", pair], meanRSE)/samplingTime)
    B1p <- predict(Loess.fit, tRunAVec[1]); B2p <- predict(Loess.fit, tRunAVec[length(tRunAVec)])
    # Alignobj <- alignChromatograms_cpp(r1, r2, tRunAVec, tRunBVec, "mean", simMeasure, B1p, B2p, noBeef)
    s1 <- getChromSimMat(r1, r2, Normalization = "mean", SimType = simMeasure)
    gapPenalty <- getGapPenalty(s1, gapQuantile, type = simMeasure)
    Mask <- getGlobalAlignMask(tRunAVec, tRunBVec, B1p, B2p, noBeef, FALSE)
    s1 <- constrainSimMain(s1, Mask, samples4gradient)
    Alignobj <- doAffineAlignment_S4(s, nrow(s), ncol(s), go = gapPenalty*goFactor, ge = gapPenalty*geFactor, OverlapAlignment = TRUE)
    AlignedIndices <- cbind(Alignobj@indexA_aligned, Alignobj@indexB_aligned, Alignobj@score)
    colnames(AlignedIndices) <- c("indexA_aligned", "indexB_aligned", "score")
    AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
    tA <- StrepChroms[[run_pair[1]]][[peptide]][[1]][["time"]]
    tB <- StrepChroms[[run_pair[2]]][[peptide]][[1]][["time"]]
    tA.aligned <- mapIdxToTime(tA, AlignedIndices[,"indexA_aligned"])
    tB.aligned <- mapIdxToTime(tB, AlignedIndices[,"indexB_aligned"])
    predictTime <- tB.aligned[which.min(abs(tA.aligned - StrepAnnot[peptide, run_pair[1]]))]
    deltaT <- predictTime - StrepAnnot[peptide, run_pair[2]]
    Err[peptide, 1] <- deltaT
  }
sum(abs(Err[,1])) # 49.2
***/
