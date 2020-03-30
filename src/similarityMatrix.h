#ifndef SIMILARITY_MATRIX_H
#define SIMILARITY_MATRIX_H

#include <vector>

namespace DIAlign
{
  /**
     @brief Similarity matrix

     A similarity matrix is an n x m matrix that is computed from a vector of
     length n and a vector of length m using a given similarity measure. See getSimilarityMatrix().
  */
  struct SimMatrix
  {
    std::vector<double> data; ///< Similarity data
    int n_row;
    int n_col;
  };

  /**
   @brief Path matrix

   Boolean similarity matrix is used to store the alignment path as true-hot encoding.
   */
  // TODO: It should be called Path matrix.
  struct SimMatrix_bool
  {
    std::vector<bool> data;
    int n_row;
    int n_col;
  };
} // namespace DIAlign

#endif // SIMILARITY_MATRIX_H
