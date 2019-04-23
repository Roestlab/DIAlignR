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
  if(alignObj.FreeEndGaps){
    for(int i = 0; i<=signalA_len; i++){
      M(i, 0) = 0;
      Traceback[i*(signalB_len+1)+0] = TM; //Top
    }
    for(int j = 0; j<=signalB_len; j++){
      M(0, j) = 0;
      Traceback[0*(signalB_len+1)+j] = LM; //Left
    }
    Traceback[0*(signalB_len+1) + 0] = SS; //STOP
  }
  else{
    for(int i = 0; i<=signalA_len; i++){
      M(i, 0) = -i*gap;
      Traceback[i*(signalB_len+1)+0] = TM; //Top
    }
    for(int j = 0; j<=signalB_len; j++){
      M(0, j) = -j*gap;
      Traceback[0*(signalB_len+1)+j] = LM; //Left
    }
    Traceback[0*(signalB_len+1) + 0] = SS; //STOP
  }

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

  if(alignObj.FreeEndGaps){
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
    TracebackPointer = alignObj.Traceback[ROW_IDX*COL_SIZE+COL_IDX];
    if(ROW_IDX != alignObj.signalA_len){
      for (int i = alignObj.signalA_len; i>ROW_IDX; i--){
        alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), i);
        alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
        alignedIdx.score.insert(alignedIdx.score.begin(), maxScore);
      }
    }
    else if (COL_IDX != alignObj.signalB_len){
      for (int j = alignObj.signalB_len; j>COL_IDX; j--){
        alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
        alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), j);
        alignedIdx.score.insert(alignedIdx.score.begin(), maxScore);
      }
    }
  }
  else{
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
  alignObj.indexA_aligned = alignedIdx.indexA_aligned;
  alignObj.indexB_aligned = alignedIdx.indexB_aligned;
  alignObj.score = alignedIdx.score;
  return;
}
