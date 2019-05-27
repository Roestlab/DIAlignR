#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include "simpleFcn.h"
#include <limits>

// It performs affine alignment on similarity matrix and fills three matrices M, A and B, and corresponding traceback matrices.
AffineAlignObj doAffineAlignment(SimMatrix s, int signalA_len, int signalB_len, double go, double ge, bool OverlapAlignment);

void getAffineAlignedIndices(AffineAlignObj &affineAlignObj);

template<class T>
double getOlapAffineAlignStartIndices(T MatrixM, T MatrixA, T MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);

#endif // AFFINEALIGNMENT_H

