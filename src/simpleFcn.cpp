#include "simpleFcn.h"

namespace DIAlign 
{

#ifdef USE_Rcpp
// TODO Write unit test for it.
// Outputs a NumericMatrix of given row and column size.
NumericMatrix initializeMatrix(double initVal, int ROW_SIZE, int COL_SIZE){
  NumericMatrix s(ROW_SIZE, COL_SIZE);
  for(int i = 0; i < ROW_SIZE; i++){
    for(int j = 0; j < COL_SIZE; j++){
      s(i, j) = initVal;
    }
  }
  return s;
}
#endif

SimMatrix getseqSim(std::string seq1, std::string seq2, double match, double misMatch){
  SimMatrix s;
  s.n_row = seq1.size();
  s.n_col = seq2.size();
  s.data.resize(seq1.size()*seq2.size(), 0.0);
  for(int i = 0; i < s.n_row; i++)
    for(int j = 0; j < s.n_col; j++)
      seq1[i] == seq2[j] ?  s.data[i*s.n_col+j] = match : s.data[i*s.n_col+j] = misMatch;
  return(s);
}

} // namespace DIAlign
