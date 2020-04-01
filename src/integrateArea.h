#ifndef INTEGRATEAREA_H
#define INTEGRATEAREA_H

#include <vector>
#include <cmath>
#include "utils.h"
#include "similarityMatrix.h"

namespace DIAlign
{
  /**
   * @brief returns the summation of signals between leftIdx and rightIdx from vov.
   *
   */
  double areaBwBoundaries(std::vector<std::vector<double> > vov, int leftIdx, int rightIdx);
} // namespace DIAlign

#endif // INTEGRATEAREA_H
