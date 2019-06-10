#include "gapPenalty.h"

double getGapPenalty(const SimMatrix& s, double gapQuantile, std::string SimType){
  double gapPenalty;
  if (SimType == "dotProductMasked")
    gapPenalty = getQuantile(s.data, gapQuantile);
  else if (SimType == "dotProduct")
    gapPenalty = getQuantile(s.data, gapQuantile);
  else if(SimType == "cosineAngle")
    gapPenalty = 0.95;
  else if(SimType == "cosine2Angle")
    gapPenalty = 0.95;
  else if(SimType == "euclideanDist")
    gapPenalty = getQuantile(s.data, gapQuantile);
  else if(SimType == "covariance")
    gapPenalty = getQuantile(s.data, gapQuantile);
  else if(SimType == "correlation")
    gapPenalty = getQuantile(s.data, gapQuantile);
  else
    Rcpp::Rcout << "getChromSimMat should have value from given choices only!" << std::endl;
  return gapPenalty;
}
