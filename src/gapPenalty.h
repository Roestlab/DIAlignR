#ifndef GAPPENALTY_H
#define GAPPENALTY_H

#include <vector>
#include "utils.h"
#include "similarityMatrix.h"

namespace DIAlign
{
  /**
   * @brief returns a gap penalty from the distribution of similarity scores.
   *
   * For simType = "cosineAngle" and "cosine2Angle", gap penalty = 0.95
   * Other simType values, gap Penalty is the gapQuantile quantile from the score disctribution of s.
   * gapQuantile must be between 0 and 1.
   */
  double getGapPenalty(const SimMatrix& s, double gapQuantile, std::string SimType);
} // namespace DIAlign

#endif // GAPPENALTY_H
