#ifndef CONSTRAINMAT_H
#define CONSTRAINMAT_H

#include <vector>
#include <cmath>
#include "utils.h"
#include "similarityMatrix.h"

namespace DIAlign 
{
void calcNoBeefMask(SimMatrix& MASK, double A1, double A2, double B1, double B2, double B1p, double B2p, int noBeef, bool hardConstrain);

void constrainSimilarity(SimMatrix& s, const SimMatrix& MASK, double constrainVal);
} // namespace DIAlign

#endif // CONSTRAINMAT_H
