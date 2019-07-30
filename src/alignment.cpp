#include "alignment.h"

namespace DIAlign
{

// This function performs dynamic programming and calculates "M" and "Traceback". Traceback matrix keeps record of the path as we fill matrix M.
AlignObj doAlignment(SimMatrix s, double gap, bool OverlapAlignment){
  int signalA_len = s.n_row;
  int signalB_len = s.n_col;
  AlignObj alignObj(signalA_len+1, signalB_len+1); // initialize AlignObj struct
  alignObj.FreeEndGaps = OverlapAlignment;
  alignObj.GapOpen = gap;
  alignObj.GapExten = gap;
  alignObj.s_data = s.data;

  SimMatrix M;
  M.n_row = signalA_len+1;
  M.n_col = signalB_len+1;
  M.data.resize(M.n_row*M.n_col, 0.0);// get a NumericMatrix filled with zeros.
  SimMatrix M_forw;
  M_forw.n_row = signalA_len+1;
  M_forw.n_col = signalB_len+1;
  M_forw.data.resize(M_forw.n_row*M_forw.n_col, 0.0);// get a NumericMatrix filled with zeros.
  std::vector<TracebackType> Traceback;
  Traceback.resize((signalA_len+1)*(signalB_len+1), SS); // Fill Traceback matrix with SS.
  std::vector<int> OptionalPaths;
  OptionalPaths.resize((signalA_len+1)*(signalB_len+1), 0);

  // Initialize first row and first column for global and overlap alignment.
  if(alignObj.FreeEndGaps){
    // For Overlap alignment, First row and first column of M matrix is filled with zeros.
    for(int i = 0; i<=signalA_len; i++){
      M.data[i*M.n_col + 0] = 0;
      M_forw.data[i*M_forw.n_col + 0] = 0;
      Traceback[i*(signalB_len+1)+0] = TM; //Top. First column is filled with TM
      OptionalPaths[i*(signalB_len+1)+0] = 1;
    }
    for(int j = 0; j<=signalB_len; j++){
      M.data[0*M.n_col + j] = 0;
      M_forw.data[0*M_forw.n_col + j] = 0;
      Traceback[0*(signalB_len+1)+j] = LM; //Left. First row is filled with LM
      OptionalPaths[0*(signalB_len+1)+j] = 1;
    }
    Traceback[0*(signalB_len+1) + 0] = SS; //STOP. Top-Left cell of the traceback matrix indicates stop
  }
  else{
    // For global alignment, top-row and left-column cells are filled with values indicating distance from top-left corner.
    for(int i = 0; i<=signalA_len; i++){
      M.data[i*M.n_col + 0] = -i*gap;
      M_forw.data[i*M_forw.n_col + 0] = -i*gap;
      Traceback[i*(signalB_len+1)+0] = TM; //Top. First column is filled with TM
      OptionalPaths[i*(signalB_len+1)+0] = 1;
    }
    for(int j = 0; j<=signalB_len; j++){
      M.data[0*M.n_col + j] = -j*gap;
      M_forw.data[0*M_forw.n_col + j] = -j*gap;
      Traceback[0*(signalB_len+1)+j] = LM; //Left. First row is filled with LM
      OptionalPaths[0*(signalB_len+1)+j] = 1;
    }
    Traceback[0*(signalB_len+1) + 0] = SS; //STOP. Top-Left cell of the traceback matrix indicates stop
  }

  // Perform dynamic programming for alignment
  double Diago, gapInA, gapInB; // signalA is along the rows, signalB is along the columns
  for(int i=1; i<=signalA_len; i++ ){
    for(int j=1; j<=signalB_len; j++){
      Diago = M.data[(i-1)*M.n_col + j-1] + s.data[(i-1)*s.n_col + j-1];
      gapInB= M.data[(i-1)*M.n_col + j] - gap; // Travelling from Top. signalA is along the rows, signalB is along the columns
      gapInA = M.data[i*M.n_col + j-1] - gap; // Travelling from Left. signalA is along the rows, signalB is along the columns
      M_forw.data[i*M_forw.n_col + j] = Diago + gapInA + gapInB;
      int optimalPathCntr = 0;
      if(gapInA>=Diago && gapInA>=gapInB){
        Traceback[i*(signalB_len+1) + j] = LM; // L: Left
        M.data[i*M.n_col + j] = gapInA;
        optimalPathCntr += OptionalPaths[i*(signalB_len+1) + j-1];
      }
      if(gapInB>=Diago && gapInB>=gapInA){
        Traceback[i*(signalB_len+1) + j] = TM; // T: Top
        M.data[i*M.n_col + j] = gapInB;
        optimalPathCntr += OptionalPaths[(i-1)*(signalB_len+1) + j];
      }
      if(Diago>=gapInA && Diago>=gapInB){
        Traceback[i*(signalB_len+1) + j] = DM; // D: Diagonal
        M.data[i*M.n_col + j] = Diago;
        optimalPathCntr += OptionalPaths[(i-1)*(signalB_len+1) + j-1];
      }
      OptionalPaths[i*(signalB_len+1) + j] = optimalPathCntr;
    }
  }
  alignObj.Traceback = Traceback; // Copy traceback to alignObj
  alignObj.OptionalPaths = OptionalPaths; // Copy OptionalPaths to alignObj
  for (int i = 0; i < signalA_len+1; i++) {
    for (int j = 0; j < signalB_len+1; j++) {
      alignObj.M[i*(signalB_len+1) + j] = M.data[i*M.n_col + j]; // Copy NumericMatrix M to alignObj.M vector
      alignObj.M_forw[i*(signalB_len+1) + j] = M_forw.data[i*M_forw.n_col + j]; // Copy NumericMatrix M_forw to alignObj.M_forw vector
    }
  }
  // printMatrix(alignObj.M, signalA_len+1, signalB_len+1);
  // printMatrix(alignObj.Traceback, signalA_len+1, signalB_len+1);
  return alignObj;
}

// This tracebacks along the highest scoring path, preparing list of scores and aligned indices.
// It calculates row and column index pair associated with the highest scoring path through similarity
// matrix. Output is a list of row-column index pairs of all highest scoring traceback paths.
void getAlignedIndices(AlignObj &alignObj){
  AlignedIndices alignedIdx; // initialize empty struct
  TracebackType TracebackPointer;
  int ROW_IDX = alignObj.signalA_len;
  int COL_IDX = alignObj.signalB_len;
  int ROW_SIZE = alignObj.signalA_len + 1;
  int COL_SIZE = alignObj.signalB_len + 1;

  if(alignObj.FreeEndGaps){
    // For overlap alignment find maximum score in matrix "M" along last column and long row.
    float maxScore = -std::numeric_limits<float>::infinity();
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
      if(alignObj.M[i*COL_SIZE+COL_SIZE-1] >= maxScore){
        ROW_IDX = i;
        COL_IDX = COL_SIZE-1;
        maxScore = alignObj.M[i*COL_SIZE+COL_SIZE-1];
      }
    }
    for (int j = 0; j < COL_SIZE; j++){
      if(alignObj.M[(ROW_SIZE-1)*COL_SIZE+j] >= maxScore){
        ROW_IDX = ROW_SIZE-1;
        COL_IDX = j;
        maxScore = alignObj.M[(ROW_SIZE-1)*COL_SIZE+j];
      }
    }
    // Use the indices of maximum score as starting point for traceback procedure.
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
    alignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
    // Align indices higher than max-score-index to NA
    if(ROW_IDX != alignObj.signalA_len){
      // Maximum score is obtained in last column. Align all row indices below max-score-index to NA.
      for (int i = alignObj.signalA_len; i>ROW_IDX; i--){
        alignedIdx.indexA_aligned.push_back(i);
        alignedIdx.indexB_aligned.push_back(NA); // Insert NA in signalB.
        alignedIdx.score.push_back(maxScore); // Insert maxScore instead of score from the matrix M.
        alignObj.Path[i*COL_SIZE+COL_IDX] = true;
      }
    }
    else if (COL_IDX != alignObj.signalB_len){
      // Maximum score is obtained in last row. Align all column indices right to max-score-index to NA.
      for (int j = alignObj.signalB_len; j>COL_IDX; j--){
        alignedIdx.indexA_aligned.push_back(NA); // Insert NA in signalA.
        alignedIdx.indexB_aligned.push_back(j);
        alignedIdx.score.push_back(maxScore); // Insert maxScore instead of score from the matrix M.
        alignObj.Path[ROW_IDX*COL_SIZE+j] = true;
      }
    }
  }
  else{
    // In Global alignment, traceback starts at the bottom-right corner.
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
    alignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
  }
  alignObj.score_forw = alignObj.M_forw[ROW_IDX*COL_SIZE+COL_IDX];

  while(TracebackPointer != SS){
    // SS: STOP when top-left corner of the matrix is reached
    // D: Diagonal, T: Top, L: Left
    // TODO Have statement for all conditions.
    switch(TracebackPointer){
    case DM: {
      // Go diagonal (Up-Left) in the matrix M.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      alignedIdx.score.push_back(alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      alignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      ROW_IDX = ROW_IDX - 1;
      COL_IDX = COL_IDX - 1;
      break;
    }
    case TM:
    {
      // Go up in the matrix M.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(NA);
      alignedIdx.score.push_back(alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      alignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      ROW_IDX = ROW_IDX - 1;
      break;
    }
    case LM:
    {
      // Go left in the matrix M.
      alignedIdx.indexA_aligned.push_back(NA);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      alignedIdx.score.push_back(alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      alignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      COL_IDX = COL_IDX - 1;
      break;
    }
    }
    // Read traceback for the next iteration.
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
  }
  // push_back adds values at the end of vector, therefore, reverse the vector.
  std::reverse(std::begin(alignedIdx.indexA_aligned), std::end(alignedIdx.indexA_aligned));
  std::reverse(std::begin(alignedIdx.indexB_aligned), std::end(alignedIdx.indexB_aligned));
  std::reverse(std::begin(alignedIdx.score), std::end(alignedIdx.score));
  // Copy aligned indices to alignObj.
  alignObj.indexA_aligned = alignedIdx.indexA_aligned;
  alignObj.indexB_aligned = alignedIdx.indexB_aligned;
  alignObj.score = alignedIdx.score;
}
} // namespace DIAlign
