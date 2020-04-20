#include <Rcpp.h>
#include "simpleFcn.h"
#include "interface.h"
#include "chromSimMatrix.h"
#include "alignment.h"
#include "gapPenalty.h"
#include "affinealignobj.h"
#include "affinealignment.h"
#include "constrainMat.h"
#include "integrateArea.h"
#include "PeakIntegrator.h"
#include "MSChromatogram.h"
#include "ChromatogramPeak.h"
#include "DPosition.h"
using namespace Rcpp;
using namespace DIAlign;
using namespace AffineAlignment;
using namespace Alignment;
using namespace SimilarityMatrix;
using namespace ConstrainMatrix;
using namespace Utils;
using namespace Traceback;
using namespace PeakGroupIntensity;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//' Calculates similarity matrix for two sequences
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-05
//' @param seq1 (char) A single string.
//' @param seq2 (char) A single string.
//' @param match (double) Score for character match.
//' @param misMatch (double) score for character mismatch.
//' @return s (matrix) Numeric similarity matrix. Rows and columns expresses seq1 and seq2, respectively.
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' matrix(c(-2, 10, -2, -2, -2, -2, 10, -2, 10, -2, -2, -2, -2, -2, -2, 10, 10, -2, -2, -2),
//'  4, 5, byrow = FALSE)
//' @export
// [[Rcpp::export]]
NumericMatrix getSeqSimMatCpp(std::string seq1, std::string seq2, double match, double misMatch){
  SimMatrix s = getseqSim(seq1, seq2, match, misMatch);
  NumericMatrix sim = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
  return(sim);
}

//' Calculates similarity matrix of two fragment-ion chromatogram groups or extracted-ion chromatograms(XICs)
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-05
//' @param l1 (list) A list of vectors. Length should be same as of l2.
//' @param l2 (list) A list of vectors. Length should be same as of l1.
//' @param normalization (char) A character string. Normalization must be selected from (L2, mean or none).
//' @param simType (char) A character string. Similarity type must be selected from (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation).\cr
//' Mask = s > quantile(s, dotProdThresh)\cr
//' AllowDotProd= [Mask × cosine2Angle + (1 - Mask)] > cosAngleThresh\cr
//' s_new= s × AllowDotProd
//' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
//' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
//' @return s (matrix) Numeric similarity matrix. Rows and columns expresses seq1 and seq2, respectively.
//' @examples
//' # Get similarity matrix of dummy chromatograms
//' r1 <- list(c(1.0,3.0,2.0,4.0), c(0.0,0.0,0.0,1.0), c(4.0,4.0,4.0,5.0))
//' r2 <- list(c(1.4,2.0,1.5,4.0), c(0.0,0.5,0.0,0.0), c(2.0,3.0,4.0,0.9))
//' round(getChromSimMatCpp(r1, r2, "L2", "dotProductMasked"), 3)
//' matrix(c(0.125, 0.162, 0.144, 0.208, 0.186, 0.240,
//' 0.213, 0.313, 0.233, 0.273, 0.253, 0.346, 0.101, 0.208, 0.154, 0.273), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "L2", "dotProduct"), 3)
//' matrix(c(0.125, 0.162, 0.144, 0.208, 0.186,0.240, 0.213, 0.313, 0.233,
//' 0.273, 0.253, 0.346, 0.101, 0.208, 0.154, 0.273), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "L2", "cosineAngle"), 3)
//' matrix(c(0.934, 0.999, 0.989, 0.986, 0.933, 0.989,
//'  0.983, 0.996, 0.994, 0.960, 0.995, 0.939, 0.450,
//'  0.761, 0.633, 0.772), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "L2", "cosine2Angle"), 3)
//' matrix(c(0.744, 0.998, 0.957, 0.944, 0.740, 0.956, 0.932,
//' 0.985, 0.974, 0.842, 0.978, 0.764, -0.596, 0.158,
//' -0.200, 0.190), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "mean", "euclideanDist"), 3)
//' matrix(c(0.608, 0.614, 0.680, 0.434, 0.530, 0.742,
//' 0.659, 0.641, 0.520, 0.541, 0.563, 0.511, 0.298,
//' 0.375, 0.334, 0.355), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "L2", "covariance"), 3)
//' matrix(c(0.025, 0.028, 0.027, 0.028, 0.032, 0.034,
//' 0.033, 0.034, 0.055, 0.051, 0.053, 0.051,
//' -0.004, 0.028, 0.012, 0.028), 4, 4, byrow = FALSE)
//'
//' round(getChromSimMatCpp(r1, r2, "L2", "correlation"), 3)
//' matrix(c(0.874, 0.999, 0.974, 0.999, 0.923, 0.986, 0.993,
//' 0.986, 0.991, 0.911, 0.990, 0.911, -0.065, 0.477,
//' 0.214, 0.477), 4, 4, byrow = FALSE)
//' @export
// [[Rcpp::export]]
NumericMatrix getChromSimMatCpp(Rcpp::List l1, Rcpp::List l2, std::string normalization, std::string simType, double cosAngleThresh = 0.3, double dotProdThresh = 0.96){
  // C++ code is compatible with vector of vector input, therefore, converting list to vector of vectors.
  std::vector<std::vector<double> > r1 = list2VecOfVec(l1);
  std::vector<std::vector<double> > r2 = list2VecOfVec(l2);
  SimMatrix s = getSimilarityMatrix(r1, r2, normalization, simType, cosAngleThresh, dotProdThresh);
  // printVecOfVec(l1);
  // printMatrix(s.data, s.n_row, s.n_col);
  // Proceessing in R is done with matrix-format, therefore, converting STL vector into NumericMatrix.
  NumericMatrix simMat = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
  return simMat;
}

//' Outputs a mask for constraining similarity matrix
//'
//' This function takes in timeVectors from both runs, global-fit mapped values
//' of end-points of first time vector and sample-length of window of no constraining.
//' Outside of window, all elements of matrix are either equally weighted or weighted
//' proportional to distance from window-boundry.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param tA (numeric) A numeric vector. This vector has equally spaced timepoints of XIC A.
//' @param tB (numeric) A numeric vector. This vector has equally spaced timepoints of XIC B.
//' @param B1p (numeric) Timepoint mapped by global fit for tA[1].
//' @param B2p (numeric) Timepoint mapped by global fit for tA[length(tA)].
//' @param noBeef (integer) It defines the distance from the global fit, upto which no penalization is performed.\cr
//' The window length = 2*noBeef.
//' @param hardConstrain (logical) if false; indices farther from noBeef distance are filled with distance from linear fit line.
//' @return mask (matrix) A numeric matrix.
//' @examples
//' tA <- c(3353.2, 3356.6, 3360.0, 3363.5)
//' tB <- c(3325.9, 3329.3, 3332.7, 3336.1)
//' B1p <- 3325.751; B2p <- 3336.119
//' noBeef <- 1
//' mask <- getGlobalAlignMaskCpp(tA, tB, B1p, B2p, noBeef, FALSE)
//' round(mask, 3)
//' matrix(c(0.000, 0.000, 0.707, 1.414, 0.000, 0.000, 0.000, 0.707, 0.707, 0.000,
//' 0.000, 0.000, 1.414, 0.707, 0.000, 0.000), 4, 4, byrow = FALSE)
//' @export
// [[Rcpp::export]]
NumericMatrix getGlobalAlignMaskCpp(const std::vector<double>& tA, const std::vector<double>& tB,
                                 double B1p, double B2p, int noBeef = 50, bool hardConstrain = false){
  SimMatrix MASK; // Initializing MASK
  MASK.n_row = tA.size();
  MASK.n_col = tB.size();
  MASK.data.resize(MASK.n_row*MASK.n_col, 0.0);
  double A1 = tA[0], A2 = tA[MASK.n_row-1];
  double B1 = tB[0], B2 = tB[MASK.n_col-1];
  calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, hardConstrain);
  // Proceessing in R is done with matrix-format, therefore, converting STL vector into NumericMatrix.
  NumericMatrix mask = Vec2NumericMatrix(MASK.data, MASK.n_row, MASK.n_col);
  return mask;
}

//' Constrain similarity matrix with a mask
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param sim (matrix) A numeric matrix. Input similarity matrix.
//' @param MASK (matrix) A numeric matrix. Masked indices have non-zero values.
//' @param samples4gradient (numeric) This paarameter modulates penalization of masked indices.
//' @return s_new (matrix) A constrained similarity matrix.
//' @examples
//' sim <- matrix(c(-2, 10, -2, -2, -2, -2, 10, -2, 10, -2, -2, -2, -2, -2, -2, 10, 10, -2,-2, -2),
//'  4, 5, byrow = FALSE)
//' MASK <- matrix(c(0.000, 0.000, 0.707, 1.414, 0.000, 0.000, 0.000, 0.707, 0.707, 0.000,
//' 0.000, 0.000, 1.414, 0.707, 0, 0, 2.121, 1.414, 0, 0), 4, 5, byrow = FALSE)
//' constrainSimCpp(sim, MASK, 10)
//' matrix(c(-2, 10, -3.414, -4.828, -2, -2, 10, -3.414, 8.586, -2, -2, -2, -4.828,
//' -3.414, -2, 10, 5.758, -4.828, -2, -2), 4, 5, byrow = FALSE)
//' @export
// [[Rcpp::export]]
NumericMatrix constrainSimCpp(const NumericMatrix& sim, const NumericMatrix& MASK, double samples4gradient = 100.0){
  SimMatrix s = NumericMatrix2Vec(sim); // converting NumericMatrix to STL vector because of C++ compatibility.
  SimMatrix mask = NumericMatrix2Vec(MASK);
  // Calculating maximum value of similarity matrix.
  auto maxIt = max_element(std::begin(s.data), std::end(s.data),
                           [](double const x, double const y) {return std::abs(x) < std::abs(y);});
  // This is a lambda function because it doesn't have a name and is anonymous. [] is a capture group.
  // [=] defines capture by value. [&] defines capture by reference.
  double maxVal = *maxIt;
  constrainSimilarity(s, mask, -2.0*maxVal/samples4gradient);
  // Proceessing in R is done with matrix-format, therefore, converting STL vector into NumericMatrix.
  NumericMatrix s_new = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
  return s_new;
}

//' Calculates gap penalty for dynamic programming based alignment.
//'
//' This function outputs base gap-penalty depending on SimType used. In case of getting base gap-penalty
//' from similarity matrix distribution, gapQuantile will be used to pick the value.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param sim (matrix) A numeric matrix. Input similarity matrix.
//' @param SimType (char) A character string. Similarity type must be selected from (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation).
//' @param gapQuantile (numeric) Must be between 0 and 1.
//' @return baseGapPenalty (numeric).
//' @examples
//' sim <- matrix(c(-12, 1.0, 12, -2.3, -2, -2, 1.07, -2, 1.80, 2, 22, 42, -2, -1.5, -2, 10), 4, 4,
//'  byrow = FALSE)
//' getBaseGapPenaltyCpp(sim, "dotProductMasked", 0.5) # -0.25
//' @export
// [[Rcpp::export]]
double getBaseGapPenaltyCpp(const NumericMatrix& sim, std::string SimType, double gapQuantile = 0.5){
  SimMatrix s = NumericMatrix2Vec(sim); // converting NumericMatrix to STL vector because of C++ compatibility.
  double gapPenalty = getGapPenalty(s, gapQuantile, SimType);
  return gapPenalty;
}

//' Calculates area between signal-boundaries.
//'
//' This function sums all the intensities between left-index and right-index.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param l1 (list) A list of time vectors.
//' @param l2 (list) A list of intensity vectors.
//' @param left (numeric) left boundary of the peak.
//' @param right (numeric) right boundary of the peak.
//' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
//'  from "intensity_sum", "trapezoid", "simpson".
//' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
//'  from "base_to_base", "vertical_division_min", "vertical_division_max".
//' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
//' @return area (numeric).
//' @examples
//' data("XIC_QFNNTDIVLLEDFQK_3_DIAlignR", package = "DIAlignR")
//' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
//' l1 <- lapply(XICs, `[[`, 1)
//' l2 <- lapply(XICs, `[[`, 2)
//' areaIntegrator(l1, l2, left = 5203.7, right = 5268.5, "intensity_sum", "base_to_base", FALSE)
//' # 224.9373
//' @export
// [[Rcpp::export]]
double areaIntegrator(Rcpp::List l1, Rcpp::List l2, double left, double right,
                      std::string integrationType, std::string baselineType, bool fitEMG){
  std::vector<std::vector<double> > vov1 = list2VecOfVec(l1);
  std::vector<std::vector<double> > vov2 = list2VecOfVec(l2);
  if(std::isnan(left) or std::isnan(right)) return NumericVector::get_na();
  std::vector<std::vector<double> > set = peakGroupArea(vov1, vov2, left, right, integrationType, baselineType, fitEMG= false);
  double area = 0.0;
  for(auto const& i : set[0]) area+= i;
  return area;
}

//' Aligns MS2 extracted-ion chromatograms(XICs) pair.
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param l1 (list) A list of numeric vectors. l1 and l2 should have same length.
//' @param l2 (list) A list of numeric vectors. l1 and l2 should have same length.
//' @param alignType (char) A character string. Available alignment methods are "global", "local" and "hybrid".
//' @param tA (numeric) A numeric vector. This vector has equally spaced timepoints of XIC A.
//' @param tB (numeric) A numeric vector. This vector has equally spaced timepoints of XIC B.
//' @param normalization (char) A character string. Normalization must be selected from (L2, mean or none).
//' @param simType (char) A character string. Similarity type must be selected from (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation).\cr
//' Mask = s > quantile(s, dotProdThresh)\cr
//' AllowDotProd= [Mask × cosine2Angle + (1 - Mask)] > cosAngleThresh\cr
//' s_new= s × AllowDotProd
//' @param B1p (numeric) Timepoint mapped by global fit for tA[1].
//' @param B2p (numeric) Timepoint mapped by global fit for tA[length(tA)].
//' @param noBeef (integer) It defines the distance from the global fit, upto which no penalization is performed.\cr
//' The window length = 2*noBeef.
//' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
//' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
//' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
//' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
//' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
//' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
//' @param hardConstrain (logical) if false; indices farther from noBeef distance are filled with distance from linear fit line.
//' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
//' @param objType (char) A character string. Must be either light, medium or heavy.
//' @return affineAlignObj (S4class) A S4class object from C++ AffineAlignObj struct.
//' @examples
//' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
//' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
//' data(oswFiles_DIAlignR, package="DIAlignR")
//' oswFiles <- oswFiles_DIAlignR
//' XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
//' XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
//' tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
//' tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
//' B1p <- 4964.752
//' B2p <- 5565.462
//' noBeef <- 77.82315/3.414
//' l1 <- lapply(XICs.ref, `[[`, 2)
//' l2 <- lapply(XICs.eXp, `[[`, 2)
//' AlignObj <- alignChromatogramsCpp(l1, l2, alignType = "hybrid", tA = tVec.ref, tB = tVec.eXp,
//'  normalization = "mean", simType = "dotProductMasked", B1p = B1p, B2p = B2p, noBeef = noBeef,
//'  goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3, OverlapAlignment = TRUE,
//'  dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE, samples4gradient = 100,
//'  objType = "light")
//' @export
// [[Rcpp::export]]
S4 alignChromatogramsCpp(Rcpp::List l1, Rcpp::List l2, std::string alignType,
                            const std::vector<double>& tA, const std::vector<double>& tB,
                            std::string normalization, std::string simType,
                            double B1p = 0.0, double B2p =0.0, int noBeef = 0,
                            double goFactor = 0.125, double geFactor = 40,
                            double cosAngleThresh = 0.3, bool OverlapAlignment = true,
                            double dotProdThresh = 0.96, double gapQuantile = 0.5,
                            bool hardConstrain = false, double samples4gradient = 100.0,
                            std::string objType = "heavy"){
  std::vector<std::vector<double> > r1 = list2VecOfVec(l1);
  std::vector<std::vector<double> > r2 = list2VecOfVec(l2);
  SimMatrix s = getSimilarityMatrix(r1, r2, normalization, simType, cosAngleThresh, dotProdThresh);
  double gapPenalty = getGapPenalty(s, gapQuantile, simType);
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
  doAffineAlignment(obj, s, gapPenalty*goFactor, gapPenalty*geFactor, OverlapAlignment); // Performs alignment on s matrix and returns AffineAlignObj struct
  getAffineAlignedIndices(obj, 9); // Performs traceback and fills aligned indices in AffineAlignObj struct

  if(objType == "light"){
    S4 x("AffineAlignObjLight");  // Creating an empty S4 object of AffineAlignObj class
    // Copying values to slots
    x.slot("indexA_aligned") = obj.indexA_aligned;
    x.slot("indexB_aligned") = obj.indexB_aligned;
    x.slot("score") = obj.score;
    return(x);
  } else if(objType == "medium"){
    S4 x("AffineAlignObjMedium");  // Creating an empty S4 object of AffineAlignObj class
    // Copying values to slots
    x.slot("s") = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
    x.slot("path") = transpose(NumericMatrix(s.n_col+1, s.n_row+1, obj.Path));
    x.slot("indexA_aligned") = obj.indexA_aligned;
    x.slot("indexB_aligned") = obj.indexB_aligned;
    x.slot("score") = obj.score;
    return(x);
  } else {
    S4 x("AffineAlignObj");  // Creating an empty S4 object of AffineAlignObj class
    // Copying values to slots
    x.slot("s") = Vec2NumericMatrix(s.data, s.n_row, s.n_col);
    x.slot("M")  = transpose(NumericMatrix(s.n_col+1, s.n_row+1, obj.M));
    x.slot("A")  = transpose(NumericMatrix(s.n_col+1, s.n_row+1, obj.A));
    x.slot("B")  = transpose(NumericMatrix(s.n_col+1, s.n_row+1, obj.B));
    std::vector<TracebackType> tmp(obj.Traceback, obj.Traceback + 3*(s.n_col+1) *(s.n_row+1) );
    x.slot("Traceback")  = EnumToChar(tmp);
    x.slot("path") = transpose(NumericMatrix(s.n_col+1, s.n_row+1, obj.Path));
    x.slot("signalA_len") = obj.signalA_len;
    x.slot("signalB_len") = obj.signalB_len;
    x.slot("GapOpen") = obj.GapOpen;
    x.slot("GapExten") = obj.GapExten;
    x.slot("FreeEndGaps") = obj.FreeEndGaps;
    x.slot("indexA_aligned") = obj.indexA_aligned;
    x.slot("indexB_aligned") = obj.indexB_aligned;
    x.slot("score") = obj.score;
    x.slot("simScore_forw") = getForwardSim(s, obj.simPath);
    x.slot("nGaps") = obj.nGaps;
    return(x);
  }
}

//' Perform non-affine global and overlap alignment on a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param sim (NumericMatrix) A numeric matrix with similarity values of two sequences or signals.
//' @param gap (double) Penalty for introducing gaps in alignment.
//' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
//' @return AlignObj (S4class) An object from C++ class of AlignObj.
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' s <- getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' obj_Global <- doAlignmentCpp(s, 22, FALSE)
//' slot(obj_Global, "score") # -2 -4 -6 4 -18
//' obj_Olap <- doAlignmentCpp(s, 22, TRUE)
//' slot(obj_Olap, "score") # 0 10 20 18 18 18
//'
//' Match=1; MisMatch=-1
//' seq1 = "TTTC"; seq2 = "TGC"
//' s <- getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' obj_Global <- doAlignmentCpp(s, 2, FALSE)
//' slot(obj_Global, "optionalPaths")
//' matrix(data = c(1,1,1,1,1,1,1,1,1,2,1,2,1,3,3,1,1,3,6,3), nrow = 5, ncol =4, byrow = TRUE)
//' slot(obj_Global, "M_forw")
//' matrix(data = c(0,-2,-4,-6,-2,-7,-22,-45,-4,-20,-72,-184,-6,-41,-178,-547,-8,-72,-366,-1274),
//'  nrow = 5, ncol =4, byrow = TRUE)
//' @export
// [[Rcpp::export]]
S4 doAlignmentCpp(NumericMatrix sim, double gap, bool OverlapAlignment){
  int signalA_len = sim.nrow(); // Length of signalA or sequenceA. Expresses along the rows of s
  int signalB_len = sim.ncol(); // Length of signalB or sequenceB. Expresses along the columns of s
  AlignObj obj(signalA_len+1, signalB_len+1); // Initializing C++ AlignObj struct
  SimMatrix s = NumericMatrix2Vec(sim);
  obj = doAlignment(s, gap, OverlapAlignment); // Performs alignment on s matrix and returns AlignObj struct
  getAlignedIndices(obj); // Performs traceback and fills aligned indices in AlignObj struct
  S4 x("AlignObj"); // Creating an empty S4 object of AlignObj class
  // Copying values to slots
  x.slot("s") = sim;
  x.slot("M") = Vec2NumericMatrix(obj.M, signalA_len+1, signalB_len+1);
  x.slot("Traceback") = Vec2NumericMatrix(EnumToChar(obj.Traceback), signalA_len+1, signalB_len+1);
  x.slot("path") = Vec2NumericMatrix(obj.Path, signalA_len+1, signalB_len+1);
  x.slot("optionalPaths") = Vec2NumericMatrix(obj.OptionalPaths, signalA_len+1, signalB_len+1);
  x.slot("M_forw") = Vec2NumericMatrix(obj.M_forw, signalA_len+1, signalB_len+1);
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  x.slot("indexA_aligned") = obj.indexA_aligned;
  x.slot("indexB_aligned") = obj.indexB_aligned;
  x.slot("score") = obj.score;
  x.slot("score_forw") = obj.score_forw;
  x.slot("nGaps") = obj.nGaps;
  return(x);
}


//' Perform affine global and overlap alignment on a similarity matrix
//'
//' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
//' ORCID: 0000-0003-3500-8152
//' License: (c) Author (2019) + MIT
//' Date: 2019-03-08
//' @param sim (NumericMatrix) A numeric matrix with similarity values of two sequences or signals.
//' @param go (numeric) Penalty for introducing first gap in alignment.
//' @param ge (numeric) Penalty for introducing subsequent gaps in alignment.
//' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
//' @return affineAlignObj (S4class) An object from C++ class of AffineAlignObj.
//' @examples
//' # Get sequence similarity of two DNA strings
//' Match=10; MisMatch=-2
//' seq1 = "GCAT"; seq2 = "CAGTG"
//' s <- getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignmentCpp(s, 22, 7, FALSE)
//' slot(objAffine_Global, "score") # -2  -4  -6  4 -18
//' objAffine_Olap <- doAffineAlignmentCpp(s, 22, 7, TRUE)
//' slot(objAffine_Olap, "score") # 0 10 20 18 18 18
//'
//' Match=10; MisMatch=-2
//' seq1 = "CAT"; seq2 = "CAGTG"
//' s <- getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignmentCpp(s, 22, 7, FALSE)
//' slot(objAffine_Global, "score") # 10  20  -2  -9 -11
//' objAffine_Olap <- doAffineAlignmentCpp(s, 22, 7, TRUE)
//' slot(objAffine_Olap, "score") # 10 20 18 18 18
//'
//' Match=10; MisMatch=-2
//' seq1 = "CA"; seq2 = "AG"
//' s <- getSeqSimMatCpp(seq1, seq2, Match, MisMatch)
//' objAffine_Global <- doAffineAlignmentCpp(s, 22, 7, FALSE)
//' slot(objAffine_Global, "simScore_forw") # -4
//' @export
// [[Rcpp::export]]
S4 doAffineAlignmentCpp(NumericMatrix sim, double go, double ge, bool OverlapAlignment){
  int signalA_len = sim.nrow(); // Length of signalA or sequenceA. Expresses along the rows of s.
  int signalB_len = sim.ncol(); // Length of signalB or sequenceB. Expresses along the columns of s.
  AffineAlignObj obj(signalA_len+1, signalB_len+1); // Initializing C++ AffineAlignObj struct
  SimMatrix s = NumericMatrix2Vec(sim);
  doAffineAlignment(obj, s, go, ge, OverlapAlignment);  // Performs alignment on s matrix and returns AffineAlignObj struct
  getAffineAlignedIndices(obj); // Performs traceback and fills aligned indices in AffineAlignObj struct
  S4 x("AffineAlignObj");  // Creating an empty S4 object of AffineAlignObj class
  // Copying values to slots
  x.slot("s") = sim;
  x.slot("M") = transpose(NumericMatrix(signalB_len+1, signalA_len+1, obj.M));
  x.slot("A") = transpose(NumericMatrix(signalB_len+1, signalA_len+1, obj.A));
  x.slot("B") = transpose(NumericMatrix(signalB_len+1, signalA_len+1, obj.B));
  std::vector<TracebackType> tmp(obj.Traceback, obj.Traceback + 3*(signalB_len+1) *(signalA_len+1) );
  x.slot("Traceback")  = EnumToChar(tmp);
  x.slot("path") = transpose(NumericMatrix(signalB_len+1, signalA_len+1, obj.Path));
  x.slot("signalA_len") = obj.signalA_len;
  x.slot("signalB_len") = obj.signalB_len;
  x.slot("GapOpen") = obj.GapOpen;
  x.slot("GapExten") = obj.GapExten;
  x.slot("FreeEndGaps") = obj.FreeEndGaps;
  x.slot("indexA_aligned") = obj.indexA_aligned;
  x.slot("indexB_aligned") = obj.indexB_aligned;
  x.slot("score") = obj.score;
  x.slot("simScore_forw") = getForwardSim(s, obj.simPath);
  x.slot("nGaps") = obj.nGaps;
  return(x);
}

// gnu -> gcc -> g++ compiler
// -I means include path. DNDEBUG includes debug symbols. Position-independent code (PIC): E.g. jumps would be generated as relative rather than absolute.
// -02 : Maximum optimization. -W warnings, -L path to the library,
// Wl,-z,relro Where to define read-only. Writing in read-only would be segmentation fault. Mostly it is after stack.
// -D_FORTIFY_SOURCE=2 To protect against buffer overflow.
