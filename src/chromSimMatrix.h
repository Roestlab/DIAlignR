#ifndef CHROMSIMMATRIX_H
#define CHROMSIMMATRIX_H

#include <vector>
#include <numeric>
#include "utils.h"
#include "similarityMatrix.h"

/**
 * @namespace DIAlign
 * @brief Generic namespace for all classes and functions of DIAlign
 */
namespace DIAlign
{

/** Similarity matrix namespace */
namespace SimilarityMatrix
{
  /// Calculates sum of previous result and square of current element (useful for sum of squares).
  template<typename T>
  struct square
  {
    T operator()(const T& Left, const T& Right) const
    {
      // We use this struct as binary operation function object. It should take current accumulation value (Left) and value of current element (Right).
      return (Left + Right*Right);
    }
  };

  /// Returns the average value of vector of vectors.
  double meanVecOfVec(const std::vector<std::vector<double>>& vov);

  /// Returns the sum of squares of all elemets of vector of vectors.
  double eucLenVecOfVec(const std::vector<std::vector<double>>& vov);

  /// Returns Eucledian length of the vector at each index projected by vectors in VoV.
  /// Eucledian length is calculated as the square-root of sum of squares of vector intensities at an index.
  /// All vectors must be of the same length.
  std::vector<double> perSampleEucLenVecOfVec(const std::vector<std::vector<double>>& vec);

  /// Returns the sum of squares of vector intensities at each index (All vectors must be of the same length).
  std::vector<double> perSampleSqrSumVecOfVec(const std::vector<std::vector<double>>& vec);

  /// Returns the mean of vector intensities at each index (All vectors must be of the same length).
  std::vector<double> perSampleMeanVecOfVec(const std::vector<std::vector<double>>& vec);

  /// Returns the sum of vector intensities at each index (All vectors must be of the same length).
  std::vector<double> perSampleSumVecOfVec(const std::vector<std::vector<double>>& vec);

  /// Calculates distance as Distance = Numerator/(Similarity score + offset).
  void distToSim(SimMatrix& s, double offset, double Numerator);

  /// Limits values between minValue and maxValue.
  void clamp(std::vector<double>& vec, double minValue, double maxValue);

  /// Returns a vector of vector with values divided by the output of meanVecOfVec().
  std::vector<std::vector<double>> meanNormalizeVecOfVec(const std::vector<std::vector<double>>& d);

  /// Returns a vector of vector with values divided by the output of eucLenVecOfVec().
  std::vector<std::vector<double>> L2NormalizeVecOfVec(const std::vector<std::vector<double>>& d);

  /// Returns a vector of vector with values divided by num.
  std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& vov, double num);

  /// Adds cross-correlation of a kernel(-+ halfLen) between d1 and d2 in similarity matrix s.
  void ElemWiseSumXcorr(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s, int halfKer);

  /// Adds outer prodict of d1 and d2 in similarity matrix s.
  void ElemWiseSumOuterProd(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

  /// Adds outer prodict of (d1-mean1) and (d2-mean2) in similarity matrix s.
  void ElemWiseSumOuterProdMeanSub(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s, const std::vector<double>& mean1, const std::vector<double>& mean2);

  /// Adds outer prodict of (d1-d2)*(d1-d2) in similarity matrix s.
  void ElemWiseSumOuterEucl(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

  /// Adds outer prodict of cosAng(d1,d2) in similarity matrix s.
  void ElemWiseOuterCosine(const std::vector<double>& d1, const std::vector<double>& d2, const std::vector<double>& d1_mag, const std::vector<double>& d2_mag, SimMatrix& s);

  /// Given Normalization modifies d1 and d2, and subsequently sums ElemWiseSumXcorr() of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumXcorr(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s, int kerLen);

  /// Given Normalization modifies d1 and d2, and subsequently sums ElemWiseSumOuterProd() of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumOuterProd(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

  /// Given Normalization modifies d1 and d2, and subsequently sums ElemWiseSumOuterProdMeanSub() of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumOuterCov(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

  /// Given Normalization modifies d1 and d2, and subsequently sums correlation coefficient of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumOuterCorr(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

  /// Given Normalization modifies d1 and d2, and subsequently sums ElemWiseSumOuterEucl() of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumOuterEucl(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

  /// Given Normalization modifies d1 and d2, and subsequently sums ElemWiseOuterCosine() of d1 vectors with d2 vectors (d1 and d2 must be of same length).
  void SumOuterCosine(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

  /// Returns a similarity matrix between d1 and d2 vector of vectors.
  ///
  /// First d1 and d2 are normalized, subsequently, similarity matrix is calculated with appropriate SimType.
  /// For SimType = dotProductMasked, matrix is further modified with cosAngleThresh and dotProdThresh parameters.
  /// For SimType == "cosine2Angle", matrix is constrained between -1.0 and 1.0.
  /// @param d1 corresponds to signal A. Must be of same size of d2.
  /// @param d2 corresponds to signal B. Must be of same size of d1.
  /// @param Normalization Must be from "mean", "L2", "None".
  /// @param SimType Must be from "dotProductMasked", "dotProduct", "cosineAngle", "cosine2Angle", "euclideanDist", "covariance", "correlation", "crossCorrelation".
  /// @param cosAngleThresh In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
  /// @param dotProdThresh In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
  /// @param kerLen In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
  SimMatrix getSimilarityMatrix(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2,
                                const std::string Normalization, const std::string SimType, double cosAngleThresh,
                                double dotProdThresh, int kerLen);

} // namespace SimilarityMatrix
} // namespace DIAlign

#endif // CHROMSIMMATRIX_H
