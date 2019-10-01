#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include "utils.h"
#include "similarityMatrix.h"
#include <limits>

namespace DIAlign
{
// It performs affine alignment on similarity matrix and fills three matrices M, A and B, and corresponding traceback matrices.
void doAffineAlignment(AffineAlignObj&,const SimMatrix& s, double go, double ge, bool OverlapAlignment);

void getAffineAlignedIndices(AffineAlignObj &affineAlignObj, int bandwidth = 0);

double getOlapAffineAlignStartIndices(double* MatrixM, double* MatrixA, double* MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);

void fillSimPath(bool* simPath, int bandwidth, int ROW_IDX, int COL_IDX, int ROW_SIZE, int COL_SIZE);

double getForwardSim(const SimMatrix& s, bool* simPath);
} // namespace DIAlign

#endif // AFFINEALIGNMENT_H

