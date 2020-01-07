#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "alignment.h"

#include "chromSimMatrix.h"
#include "alignment.h"
#include "gapPenalty.h"
#include "affinealignobj.h"
#include "affinealignment.h"
#include "constrainMat.h"


namespace DIAlign
{

  /*
   * Align two pairs of chromatograms derived from two LC-MS/MS experiments A and B. 
   *
   * @param obj A object of type AffineAlignObj which needs to have at least capacity of the chromatogram length
   * @param r1 The first set of chromatograms, where each chromatogram is a vector of intensities. These chromatograms are XICs derived from a single analyte in run A
   * @param r2 The second set of chromatograms, where each chromatogram is a vector of intensities. These chromatograms are XICs derived from a single analyte in run B
   * @param alignType Available alignment methods, must be one of "global", "local" and "hybrid".
   * @param tA This vector has equally spaced timepoints of the XICs from run A.
   * @param tB This vector has equally spaced timepoints of the XICs from run B.
   * @param normalization Normalization method, must be one of "L2", "mean" or "none".
   * @param simType Similarity computation method, defines how similarity between two chromatographic points is computed. Must be one of (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation).\cr
   * @param B1p Timepoint mapped by global fit for tA[0].
   * @param B2p Timepoint mapped by global fit for tA[tA.size()-1].
   * @param noBeef It defines the distance from the global fit, upto which no penalization is performed. 
   * The window length will be chosen as = 2*noBeef.
   * @param goFactor Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
   * @param geFactor Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
   * @param cosAngleThresh In simType == dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
   * @param OverlapAlignment If alignment should be performed with free end-gaps. False: Global alignment, True: overlap alignment.
   * @param dotProdThresh In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
   * @param gapQuantile Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
   * @param hardConstrain Whether a hard constraint should be used for the similarity matrix. If false; indices farther from noBeef distance are filled with distance from linear fit line.
   * @param samples4gradient This parameter modulates penalization of masked indices.
   *
   * @note r1 and r2 need to have the same length and the chromatograms need to be in the same order (e.g. r1[0] needs to be the same XIC trace as r2[0])
   * @note If no global fit is available, set B1p to tB.front() and B2p to tB.back()
*/
  void alignChromatogramsCpp( AffineAlignObj& obj,
                              const std::vector<std::vector<double> > & r1,
                              const std::vector<std::vector<double> > & r2,
                              std::string alignType,
                              const std::vector<double>& tA, const std::vector<double>& tB,
                              const std::string & normalization, const std::string& simType,
                              double B1p = 0.0, double B2p =0.0,
                              int noBeef = 0,
                              double goFactor = 0.125, double geFactor = 40,
                              double cosAngleThresh = 0.3, bool OverlapAlignment = true,
                              double dotProdThresh = 0.96, double gapQuantile = 0.5,
                              bool hardConstrain = false, double samples4gradient = 100.0)
  {
    SimMatrix s = SimilarityMatrix::getSimilarityMatrix(r1, r2, normalization, simType, cosAngleThresh, dotProdThresh);
    obj.reset(s.n_row + 1, s.n_col + 1);

    double gapPenalty = getGapPenalty(s, gapQuantile, simType);
    if (alignType == "hybrid")
    {
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
    doAffineAlignment(obj, s, gapPenalty*goFactor, gapPenalty*geFactor, OverlapAlignment); // Performs alignment on s matrix and returns AffineAlignObj struct
    getAffineAlignedIndices(obj); // Performs traceback and fills aligned indices in AffineAlignObj struct
  }

  /*
   * Calculates similarity matrix of two fragment-ion chromatogram groups or extracted-ion chromatograms (XICs) derived from two LC-MS/MS experiments A and B. 
   *
   * @param d1 The first set of chromatograms, where each chromatogram is a vector of intensities. These chromatograms are XICs derived from a single analyte in run A
   * @param d2 The second set of chromatograms, where each chromatogram is a vector of intensities. These chromatograms are XICs derived from a single analyte in run B
   * @param Normalization Normalization method, must be one of "L2", "mean" or "none".
   * @param SimType Similarity computation method, defines how similarity between two chromatographic points is computed. Must be one of (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation).\cr
   * @param cosAngleThresh In simType == dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
   * @param dotProdThresh In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
   * @return Returns a similarity matrix
*/
  SimMatrix getSimilarityMatrix(const std::vector<std::vector<double>>& d1,
      const std::vector<std::vector<double>>& d2,
      const std::string& Normalization,
      const std::string& SimType,
      double cosAngleThresh,
      double dotProdThresh)
  {
    return SimilarityMatrix::getSimilarityMatrix(d1, d2, normalization, simType, cosAngleThresh, dotProdThresh);
  }

  /*
   * Performs affine alignment on a given alignment object. 
   *
   * @note: Does not perform back-tracking, please call getAffineAlignedIndices afterwards.
   *
   * @param obj A object of type AffineAlignObj which needs to have at least capacity of the chromatogram length
   * @param s A previously computed similarity matrix
   * @param goPenalty Penalty for introducing first gap in alignment.
   * @param gePenalty Penalty for introducing subsequent gaps in alignment.
   * @param OverlapAlignment If alignment should be performed with free end-gaps. False: Global alignment, True: overlap alignment.
  */
  void doAffineAlignment(AffineAlignObj& obj, const SimMatrix& s, double goPenalty, double gePenalty, bool OverlapAlignment)
  {
    return AffineAlignment::doAffineAlignment(obj, s, goPenalty, gePenalty, OverlapAlignment);
  }

  /*
   * Compute path using AffineAlignObj
   *
   * @param obj A object of type AffineAlignObj which needs to have matrices A, B and M filled.
   * @param bandwith Compute path for multiple matrix cells (bandwith is number of cells to be computed)
   *
   * @note: this assumes that matrices A, B and M have been computed. Use after calling doAffineAlignment.
   * @note: this will compute the following members: indexA_aligned, indexB_aligned, score.
  */
  void getAffineAlignedIndices(AffineAlignObj &obj, int bandwidth = 0)
  {
    getAffineAlignedIndices(obj, bandwidth);
  }
}

