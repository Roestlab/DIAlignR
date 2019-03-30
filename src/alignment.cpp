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
    Traceback[i*(signalA_len+1)+0] = TM; //Top
  }
  for(int j = 0; j<=signalB_len; j++){
    M(0, j) = -j*gap;
    Traceback[0*(signalA_len+1)+j] = LM; //Left
  }
  Traceback[0*(signalA_len+1) + 0] = SS; //STOP

  // Perform dynamic programming for alignment
  float Diago, gapInA, gapInB;
  for(int i=1; i<=signalA_len; i++ ){
    for(int j=1; j<=signalB_len; j++){
      Diago = M(i-1, j-1) + s(i-1, j-1);
      gapInA = M(i-1, j) - gap;
      gapInB = M(i, j-1) - gap;
      if(Diago>=gapInA && Diago>=gapInB){
        Traceback[i*(signalA_len+1) + j] = DM; // D: Diagonal
        M(i, j) = Diago;
      }
      else if (gapInA>=Diago && gapInA>=gapInB){
        Traceback[i*(signalA_len+1) + j] = LM; // L: Left
        M(i, j) = gapInA;
      }
      else{
        Traceback[i*(signalA_len+1) + j] = TM; // T: Top
        M(i, j) = gapInB;
      }
    }
  }

  for (int i = 0; i < signalA_len+1; i++) {
    for (int j = 0; j < signalB_len+1; j++) {
      alignObj.M.push_back(M(i, j)); // Add an element (column) to the row
      alignObj.Traceback.push_back(Traceback[i*(signalA_len+1)+j]); // Add an element (column) to the row
    }
  }

  return alignObj;
}
