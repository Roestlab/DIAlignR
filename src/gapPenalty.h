#ifndef GAPPENALTY_H
#define GAPPENALTY_H

#include <vector>
#include "utils.h"
#include "similarityMatrix.h"

namespace DIAlign 
{
  double getGapPenalty(const SimMatrix& s, double gapQuantile, std::string SimType);
} // namespace DIAlign

#endif // GAPPENALTY_H
