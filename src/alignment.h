#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <limits>

#include "utils.h"
#include "similarityMatrix.h"
#include "affinealignobj.h"

/**
 * @namespace DIAlign
 * @brief Generic namespace for all classes and functions of DIAlign
 */
namespace DIAlign
{

const int NA = 0;

/**
 * @brief Includes aligned indices of signal A and signal B with the cumulative score.
 *
 */
struct AlignedIndices{
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<double> score;
};

/**
 * @brief An alignment object (for non-affine alignment)
 *
 * It contains similarity matrix, matrix M to store cumulative-scores for dynamic programming.
 * Traceback matrices store source matrix 'M' and direction as matrices are filled with dynamic programming.
 * Path matrix encode alignment path that results in the highest cumulative score.
 * The aligned indices are also stored for signal A and signal B.
*/
struct AlignObj
{
    std::vector<double> s_data; ///< similarity score matrix.
    std::vector<Traceback::TracebackType> Traceback; ///< Traceback matrices store source matrix name and direction as matrices are filled with dynamic programming.
    std::vector<double> M; ///< Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.
    std::vector<double> M_forw; ///< Not needed, will be removed.
    std::vector<bool> Path; ///< Path matrix would represent alignment path through similarity matrix as binary-hot encoding.
    std::vector<bool> simPath; ///< Not needed, will be removed.
    std::vector<int> OptionalPaths; ///< Highlight the number of all optimal paths.
    int signalA_len; ///< Number of data-points in signal A.
    int signalB_len; ///< Number of data-points in signal B.
    double GapOpen; ///< Penalty for Gap opening. For n consecutive gaps: Penalty = n*GapOpen.
    double GapExten;  ///< Not needed, will be removed.
    bool FreeEndGaps; ///< True for Overlap alignment.
    std::vector<int> indexA_aligned; ///< Aligned signalA indices after alignment.
    std::vector<int> indexB_aligned; ///< Aligned signalB indices after alignment.
    std::vector<double> score; ///< Cumulative score along the aligned path.
    double score_forw; ///< Not needed, will be removed.
    double simScore_forw; ///< Summation of similarity score along the alignment path-band.
    double alterAlignScore; ///< Alignment score of 2nd best peak alignment.
    int nGaps; ///< Total number of gaps in the alignment path.

    /**
     * @brief Constructor for AlignObj.
     *
     * Allocates memory for s_data, M, Traceback, Path matrices.
     * @param ROW_SIZE Number of rows in matrix M.
     * @param COL_SIZE Number of columns in matrix M.
     */
    // Not a default constructor
    AlignObj(int ROW_SIZE, int COL_SIZE)
    {
      M.resize(ROW_SIZE * COL_SIZE, 0);
      M_forw.resize(ROW_SIZE * COL_SIZE, 0);
      Traceback.resize(ROW_SIZE * COL_SIZE, Traceback::SS);
      Path.resize(ROW_SIZE * COL_SIZE, false);
      simPath.resize(ROW_SIZE * COL_SIZE, false);
      OptionalPaths.resize(ROW_SIZE * COL_SIZE, 0);
      signalA_len = ROW_SIZE-1;
      signalB_len = COL_SIZE-1;
      s_data.resize(signalA_len * signalB_len, 0.0);
      GapOpen = 0.0;
      GapExten = 0.0;
      FreeEndGaps = true;
      score_forw = 0.0;
      simScore_forw = 0.0;
      alterAlignScore = 0.0;
      nGaps = 0;
    }

    // Rule 1 Copy constructor and Rule 2 Copy assignment operator are not needed.
    // C++ auto-create them and will do the right-thing compared to manual code.
    // Rule 3 Not a default destructor
    ~AlignObj()
    {  }
};

/**
 * @brief Contains some of the core alignment functions for non-affine alignment.
 *
 * See doAlignment() and getAlignedIndices().
 */
namespace Alignment
{

/**
 * @brief This function performs dynamic programming and calculates AlignObj::M and Traceback.
 * Traceback matrix keeps record of the path as we fill matrix M.
 */
AlignObj doAlignment(SimMatrix s, double gap, bool OverlapAlignment);

/**
 * @brief This tracebacks along the highest scoring path, preparing list of scores and aligned indices.
 *
 * It calculates row and column index pair associated with the highest scoring
 * path through similarity matrix. Output is a list of row-column index pairs
 * of all highest scoring traceback paths.
 *
*/
void getAlignedIndices(AlignObj &alignObj);

}
}

#endif // ALIGNMENT_H
