#ifndef CONSTRAINMAT_H
#define CONSTRAINMAT_H

#include <vector>
#include <cmath>
#include "utils.h"
#include "similarityMatrix.h"

/**
 * @namespace DIAlign
 * @brief Generic namespace for all classes and functions of DIAlign
 */
namespace DIAlign
{
/** Namespace for constraing similarity score matrix. */
namespace ConstrainMatrix
{
/**
 * @brief Fill a diagonal strip with one through the matrix MASK.
 *
 * A1 and A2 are end-coordinates of signal A, similarly, B1 and B2 are end-coordinates for signal B.
 * B1P and B2p are mapping of A1 and A2 calculated through some external mapper.
 * The function fits a line between (A1, B1p) and (A2, B2p).
 * It further sets every value as 0 for the cells that are within distance noBeef from the line.
 * Other cells are either set as 1.0 (hardConstrain = true) or their distance from the line (hardConstrain = false).
 * In the given example below noBeef = 1, value = 0.0 is represented as X, non-zero values are represented as whitespace.
 * Only the region between (A1, B1) and (A2, B2) is in the MASK matrix.
 *
 * \verbatim
 *     |             XX| B2p
 *     |            XXX|
 *     |           XXX |
 *  B2 +---------------+
 *     |       XXX     |
 *     |      XXX      |
 *     |     XXX       |
 *     |    XXX        |
 *     |   XXX         |
 *     |  XXX          |
 *     | XXX           |
 *     |XX             | B1p
 * B1  +---------------+
 *    A1              A2
 * \endverbatim
 */
void calcNoBeefMask(SimMatrix& MASK, double A1, double A2, double B1, double B2, double B1p, double B2p, int noBeef, bool hardConstrain);

/**
 * @brief Applies the mask from calcNoBeefMask() on the similarity matrix.
 *
 * MASK and s must be of the same size. constrainVal is a penalizing factor for the MASK.
 */
void constrainSimilarity(SimMatrix& s, const SimMatrix& MASK, double constrainVal);
} // namespace ConstrainMatrix
} // namespace DIAlign

#endif // CONSTRAINMAT_H

