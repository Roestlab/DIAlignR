#include "simpleFcn.h"

namespace DIAlign 
{

// Outputs a NumericMatrix of given row and column size.
NumericMatrix initializeMatrix(float initVal, int ROW_SIZE, int COL_SIZE){
  NumericMatrix s(ROW_SIZE, COL_SIZE);
  for(int i = 0; i < ROW_SIZE; i++){
    for(int j = 0; j < COL_SIZE; j++){
      s(i, j) = initVal;
    }
  }
  return s;
}

} // namespace DIAlign
