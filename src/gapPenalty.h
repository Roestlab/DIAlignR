#ifndef GAPPENALTY_H
#define GAPPENALTY_H

#include "simpleFcn.h"

double getGapPenalty(const SimMatrix& s, double gapQuantile, std::string SimType);
#endif // GAPPENALTY_H
