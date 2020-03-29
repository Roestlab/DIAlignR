#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include "utils.h"
#include "similarityMatrix.h"
#include <limits>

/**
 * @namespace DIAlign
 * @brief Generic namespace for all classes and functions of DIAlign
 *
 */
namespace DIAlign
{
/// Includes functions for the affine-alignment through a score-matrix using gap-opening and gap-closing penalties.
namespace AffineAlignment
{
/**
 * @brief Performs affine alignment on similarity-score matrix DIAlign::SimMatrix and fills M, A, B and corresponding Traceback matrices in AffineAlignObj.
 *
 * Affine alignment uses gap-opening and gap-extension penalties to control number of consecutive gaps.
 * It uses three cumulative score matrices M, A and B. At start, the first row and column of these matrices are initialized.
 * Afterward, remaining cells of the matrices are filled using method described in [Paste the link to the main page with equations.]
 * As cumulative score is calculated for each cell, movement is captured  in Traceback with its direction (Top, Left, Diagonal) and the source matrix (M, A, B).
 *
 * @param affineAlignObj An object of class AffineAlignObj. It must be initialized with appropriate ROW_SIZE and COL_SIZE.
 * @param s similarity score matrix must be of size (ROW_SIZE - 1) * (COL_SIZE - 1).
 * @param go gap opening penalty in alignment path.
 * @param ge gap extension penalty in alignment path.
 * @param OverlapAlignment If true, end gaps are not penalized. For a global-alignment, OverlapAlignment is set false which penalizes end-gaps as regular gaps.
 *
 */
void doAffineAlignment(AffineAlignObj&, const SimMatrix& s, double go, double ge, bool OverlapAlignment);

/**
 * @brief Calculates aligned indices for source signal A and B, additionaly, builds an alignment path matrix with true-hot encoding.
 *
 * For finding the alignment path it uses cumulative scores and Traceback matrices. For non-overlap alignment,
 * cell with the highest score is identified across rightmost cells of M, A and B. For overlap alignment,
 * last column and row of the matrices are searched for the highest score cell.
 * After identifying highest-score cell, Traceback is used to find the path. Simultaneously, number of gaps are
 * calculated and path matrix is filled with true-hot encoding.
 *
 * @param affineAlignObj An object of class AffineAlignObj. Must have been operated by doAffineAlignment() function before.
 * @param bandwidth Not needed. Will be removed.
 */
void getAffineAlignedIndices(AffineAlignObj &affineAlignObj, int bandwidth = 0);

/**
 * @brief Calculates the start indices and matrix for the alignment path and returns associated score.
 *
 * Initially, the maximum score is searched along the last row and column of the three matrices.
 * The highest-scoring cell is passed-by-reference.
 * @param MatrixM pointer for matrix of type AffineAlignObj::M.
 * @param MatrixA pointer for matrix of type AffineAlignObj::A.
 * @param MatrixB pointer for matrix of type AffineAlignObj::B.
 * @param ROW_SIZE number of rows in each matrix.
 * @param COL_SIZE number of columns in each matrix.
 * @param OlapStartRow row-index of the start-cell.
 * @param OlapStartCol column-index of the start-cell.
 * @param MatrixName Name of the matrix from which start-cell is selected.
 * @return The start-cell score that serves as the cumulative score of the alignment.
 *
 * @note MatrixM, MatrixA and MatrixB should have the same number of rows and columns.
 */
double getOlapAffineAlignStartIndices(double* MatrixM, double* MatrixA, double* MatrixB,
                                      int ROW_SIZE, int COL_SIZE, int &OlapStartRow,
                                      int &OlapStartCol, Traceback::tbJump &MatrixName);

/**
 * @brief Not needed will be removed.
 *
 */
void fillSimPath(bool* simPath, int bandwidth, int ROW_IDX, int COL_IDX, int ROW_SIZE, int COL_SIZE);

/**
 * @brief Not needed will be removed.
 *
 */
double getForwardSim(const SimMatrix& s, bool* simPath);
} // namespace AffineAlignment
} // namespace DIAlign

#endif // AFFINEALIGNMENT_H

