#include "alignment.h"

AlignObj doAlignment(NumericMatrix s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment){
  AlignObj alignObj(signalA_len+1, signalB_len+1);
  alignObj.FreeEndGaps = OverlapAlignment;
  alignObj.GapOpen = gap;
  alignObj.GapExten = gap;
  alignObj.signalA_len = signalA_len;
  alignObj.signalB_len = signalB_len;

  NumericMatrix M;
  M = initializeMatrix(0, signalA_len+1, signalB_len+1);
  std::vector<TracebackType> Traceback;
  Traceback.resize((signalA_len+1)*(signalB_len+1), SS);

  // Initialize first row and first column for global and overlap alignment.
  for(int i = 0; i<=signalA_len; i++){
    M(i, 0) = -i*gap;
    Traceback[i*(signalB_len+1)+0] = TM; //Top
  }

  for(int j = 0; j<=signalB_len; j++){
    M(0, j) = -j*gap;
    Traceback[0*(signalB_len+1)+j] = LM; //Left
  }
  Traceback[0*(signalB_len+1) + 0] = SS; //STOP

  // Perform dynamic programming for alignment
  float Diago, gapInA, gapInB; // signalA is along the rows, signalB is along the columns
  for(int i=1; i<=signalA_len; i++ ){
    for(int j=1; j<=signalB_len; j++){
      Diago = M(i-1, j-1) + s(i-1, j-1);
      gapInB= M(i-1, j) - gap; // Travelling from Top. signalA is along the rows, signalB is along the columns
      gapInA = M(i, j-1) - gap; // Travelling from Left. signalA is along the rows, signalB is along the columns
      if(Diago>=gapInA && Diago>=gapInB){
        Traceback[i*(signalB_len+1) + j] = DM; // D: Diagonal
        M(i, j) = Diago;
      }
      else if (gapInA>=Diago && gapInA>=gapInB){
        Traceback[i*(signalB_len+1) + j] = LM; // L: Left
        M(i, j) = gapInA;
      }
      else{
        Traceback[i*(signalB_len+1) + j] = TM; // T: Top
        M(i, j) = gapInB;
      }
    }
  }
  alignObj.Traceback = Traceback;
  // Rcpp::Rcout << M << std::endl;
  for (int i = 0; i < signalA_len+1; i++) {
    for (int j = 0; j < signalB_len+1; j++) {
      alignObj.M[i*(signalB_len+1) + j] = M(i, j); // Add an element (column) to the row
    }
  }
  // printMatrix(alignObj.M, signalA_len+1, signalB_len+1);
  // printMatrix(alignObj.Traceback, signalA_len+1, signalB_len+1);
  return alignObj;
}

void getAlignedIndices(AlignObj &alignObj){
  AlignedIndices alignedIdx;
  TracebackType TracebackPointer;
  float alignmentScore;
  int ROW_IDX = alignObj.signalA_len;
  int COL_IDX = alignObj.signalB_len;
  int ROW_SIZE = alignObj.signalA_len + 1;
  int COL_SIZE = alignObj.signalB_len + 1;

  if(alignObj.FreeEndGaps == true){
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
  }

  while(TracebackPointer != SS){
    // D: Diagonal, T: Top, L: Left
    switch(TracebackPointer){
    case DM: {
      alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
      alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
      alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      ROW_IDX = ROW_IDX - 1;
      COL_IDX = COL_IDX - 1;
      break;
    }
    case TM:
    {
      alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
      alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
      alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      ROW_IDX = ROW_IDX - 1;
      break;
    }
    case LM:
    {
      alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
      alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
      alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      COL_IDX = COL_IDX - 1;
      break;
    }
    }
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
  }
  // alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
  // alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
  // alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
  // TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
  // Rcpp::Rcout << "Traceback is starting" << std::endl;
  // Rcpp::Rcout << TracebackPointer << std::endl;
  // Traceback path and align row indices to column indices.
  // alignedIdx.score.erase(alignedIdx.score.begin());
  // alignedIdx.indexA_aligned.erase(alignedIdx.indexA_aligned.begin());
  // alignedIdx.indexB_aligned.erase(alignedIdx.indexB_aligned.begin());
  alignObj.indexA_aligned = alignedIdx.indexA_aligned;
  alignObj.indexB_aligned = alignedIdx.indexB_aligned;
  alignObj.score = alignedIdx.score;
  return;
}
