#ifndef SIMILARITY_MATRIX_H
#define SIMILARITY_MATRIX_H

#include <vector>

namespace DIAlign 
{
  struct SimMatrix
  {
    std::vector<double> data;
    int n_row;
    int n_col;
  };

  struct SimMatrix_bool
  {
    std::vector<bool> data;
    int n_row;
    int n_col;
  };
} // namespace DIAlign

#endif // SIMILARITY_MATRIX_H
