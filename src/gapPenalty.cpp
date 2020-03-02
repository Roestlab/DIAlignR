#include "gapPenalty.h"

namespace DIAlign 
{
double getGapPenalty(const SimMatrix& s, double gapQuantile, std::string SimType){
  double gapPenalty;
  if (SimType == "dotProductMasked")
    gapPenalty = Utils::getQuantile(s.data, gapQuantile);
    // gapPenalty = 0.96;
  else if (SimType == "dotProduct")
    gapPenalty = Utils::getQuantile(s.data, gapQuantile);
  else if(SimType == "cosineAngle")
    gapPenalty = 0.95;
  else if(SimType == "cosine2Angle")
    gapPenalty = 0.95;
  else if(SimType == "euclideanDist")
    gapPenalty = Utils::getQuantile(s.data, gapQuantile);
  else if(SimType == "covariance")
    gapPenalty = Utils::getQuantile(s.data, gapQuantile);
  else if(SimType == "correlation")
    gapPenalty = Utils::getQuantile(s.data, gapQuantile);
  else
  {
    // Rcpp::Rcout << "getChromSimMat should have value from given choices only!" << std::endl;
  }
  return gapPenalty;
}
} // namespace DIAlign
