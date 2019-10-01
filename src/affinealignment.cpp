#include "affinealignment.h"
#include <exception>
#include <stdexcept>

// #include "simpleFcn.h"
// Do not inclue cpp file otherwise compiler will build the Obj through two different path.

namespace {
  void validate(DIAlign::AffineAlignObj& affineAlignObj, const DIAlign::SimMatrix& s, double go, double ge) {
    if(go < 0.0){
      throw std::invalid_argument("Gap opening penalty should be non-negative");
    }
    if(ge < 0.0){
      throw std::invalid_argument("Gap extension penalty should be non-negative");
    }
    if(affineAlignObj.signalA_len != s.n_row || affineAlignObj.signalB_len != s.n_col){
      throw std::invalid_argument("AffineAlignObj should have number of rows and columns +1 each than that of similarity matrix s.");
    }
    if(affineAlignObj.signalA_len <= 1 || affineAlignObj.signalB_len <= 1){
      throw std::invalid_argument("AffineAlignObj must have more than unit size.");
    }
    if(s.n_row <= 0 || s.n_col <= 0){
      throw std::invalid_argument("similarity matrix s must have atleast unit size.");
    }
  }
}

namespace DIAlign
{
// It performs affine alignment on similarity matrix and fills three matrices M, A and B, and corresponding traceback matrices.
void doAffineAlignment(AffineAlignObj& affineAlignObj, const SimMatrix& s, double go, double ge, bool OverlapAlignment){
  validate(affineAlignObj, s, go, ge);
  int signalA_len = s.n_row;
  int signalB_len = s.n_col;
  affineAlignObj.FreeEndGaps = OverlapAlignment;
  affineAlignObj.GapOpen = go;
  affineAlignObj.GapExten = ge;

  int Traceback_M_index = 0; // First block of Traceback vector corresponds to M matrix.
  int Traceback_A_index = 1; // Second block of Traceback vector corresponds to A matrix.
  int Traceback_B_index = 2; // Third block of Traceback vector corresponds to B matrix.

  int* oPathsM;
  oPathsM = new int[(signalA_len+1)*(signalB_len+1)];
  std::memset(oPathsM, 0, (signalA_len+1)*(signalB_len+1)*sizeof(int) );
  int* oPathsA;
  oPathsA = new int[(signalA_len+1)*(signalB_len+1)];
  std::memset(oPathsA, 0, (signalA_len+1)*(signalB_len+1)*sizeof(int) );
  int* oPathsB;
  oPathsB = new int[(signalA_len+1)*(signalB_len+1)];
  std::memset(oPathsB, 0, (signalA_len+1)*(signalB_len+1)*sizeof(int) );

  // Initialize first row and first column for affine alignment.
  double Inf = std::numeric_limits<double>::infinity();
  for(int i = 0; i<=signalA_len; i++){
    // Aligning ith character of signal A with 0th character of signal B without a gap. Not possible, hence, First column of M is initialized with -Inf.
    affineAlignObj.M[i*(signalB_len+1)+0] = -Inf;
    // Score of the best alignment between ith character of A and 0th character of B that introduces a gap in A. Not possible, hence, First column of B is initialized with -Inf.
    affineAlignObj.B[i*(signalB_len+1)+0] = -Inf;
    // It is impossible for traceback to reach the cells modified above, however, STOP alignment if it happens.
    affineAlignObj.Traceback[Traceback_M_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = SS; //STOP
    affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1)) + i*(signalB_len+1)+0] = SS; //STOP
    oPathsM[i*(signalB_len+1)+0] = 0;
    oPathsB[i*(signalB_len+1)+0] = 0;
    affineAlignObj.optionalPaths[i*(signalB_len+1)+0] = 1;
    }
  for(int j = 0; j<=signalB_len; j++){
    // Aligning 0th character of signal A with jth character of signal B without a gap. Not possible, hence, First row of M is initialized with -Inf.
    affineAlignObj.M[0*(signalB_len+1)+j] = -Inf;
    // Score of the best alignment between 0th character of A and jth character of B that results a gap in B. Not possible, hence, First row of A is initialized with -Inf.
    affineAlignObj.A[0*(signalB_len+1)+j] = -Inf;
    // It is impossible for traceback to reach the cells modified above, however, STOP alignment if it happens.
    affineAlignObj.Traceback[Traceback_M_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = SS; //STOP
    affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = SS; //STOP
    oPathsM[0*(signalB_len+1)+j] = 0;
    oPathsA[0*(signalB_len+1)+j] = 0;
    affineAlignObj.optionalPaths[0*(signalB_len+1)+j] = 1;
    }
  affineAlignObj.M[0*(signalB_len+1)+0] = 0; // Match state (0,0) should have zero to begin the alignment.
  oPathsM[0*(signalB_len+1)+0] = 1;

  // Fill up remaining cells of first row and first column for global and overlap alignment.
  if(affineAlignObj.FreeEndGaps == true){
    // For overlap alignment, there is no gap penalty for alignment of ith character of A to 0th characters of B that results a gap in B.
    // Hence, aligning ith character of A to gaps in B without any penalty.
    // Since, consecutive elements of A are aligned to consecutive gaps in B. Therefore, we remain in A matrix for such alignment.
    // Hence Traceback_matrix_A should have TA for these cells, except for Traceback_matrix_A(1,0). This cell indicates first gap and is reached from M(0,0).
    for(int i = 1; i<=signalA_len; i++){
      affineAlignObj.A[i*(signalB_len+1) + 0] = 0;
      affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = TA; //TOP A
      oPathsA[i*(signalB_len+1) + 0] = 1;
      }
    affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ 1*(signalB_len+1)+0] = TM; //TOP M
    // For overlap alignment, there is no gap penalty for alignment of zero characters of A to jth characters of B that results a gap in A.
    // Hence, aligning jth character of B to gaps in A without any penalty.
    // Since, consecutive elements of B are aligned to consecutive gaps in A. Therefore, we remain in B matrix for such alignment.
    // Hence Traceback_matrix_B should have LB for these cells, except for Traceback_matrix_B(0,1). This cell indicates first gap and is reached from M(0,0).
    for(int j = 1; j<=signalB_len; j++){
      affineAlignObj.B[0*(signalB_len+1)+j] = 0;
      affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = LB; //LEFT B
      oPathsB[0*(signalB_len+1)+j] = 1;
      }
    affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+1] = LM; //LEFT M
    }
  else {
    // In global alignment, penalty for alignment of ith character of A to 0th characters of B that results a gap in B =
    // GapOpen + (i-1)*GapExten
    // Since, consecutive elements of A are aligned to consecutive gaps in B. Therefore, we remain in A matrix for such alignment.
    // Hence Traceback_matrix_A should have TA for these cells.
    for(int i = 1; i<=signalA_len; i++){
      affineAlignObj.A[i*(signalB_len+1) + 0] = -(i-1)*ge - go;
      affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = TA; //TOP A
      oPathsA[i*(signalB_len+1) + 0] = 1;
      }
    affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ 1*(signalB_len+1)+0] = TM; //TOP M
    // In global alignment, penalty for the alignment of zero characters of A to jth characters of B that results a gap in A =
    // GapOpen + (j-1)*GapExten
    // Since, consecutive elements of B are aligned to consecutive gaps in A. Therefore, we remain in B matrix for such alignment.
    // Hence Traceback_matrix_B should have LB for these cells.
    for(int j = 1; j<=signalB_len; j++){
      affineAlignObj.B[0*(signalB_len+1)+j] = -(j-1)*ge - go;
      affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = LB; //LEFT B
      oPathsB[0*(signalB_len+1)+j] = 1;
      }
    affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+1] = LM; //LEFT M
    }

  // Perform dynamic programming to fill matrix M, A and B for affine alignment
  double Diago, InsertInA, InsertInB;
  for(int i=1; i<=signalA_len; i++){
    for(int j=1; j<=signalB_len; j++){
      // Rcpp::Rcout << s.data[(i-1)*s.n_col + j-1] << std::endl;
      double sI_1J_1 = s.data[(i-1)*s.n_col + j-1]; // signal Ai is aligned to signal Bj. Hence, it will force match state or diagonal alignment.
      Diago = affineAlignObj.M[(i-1)*(signalB_len+1)+j-1] + sI_1J_1; // M(i-1, j-1) means Ai-1 is aligned to Bj-1.
      InsertInA = affineAlignObj.A[(i-1)*(signalB_len+1)+j-1] + sI_1J_1; // A(i-1, j-1) means Ai-1 is aligned to a gap in B.
      InsertInB = affineAlignObj.B[(i-1)*(signalB_len+1)+j-1] + sI_1J_1; // B(i-1, j-1) means Bj-1 is aligned to a gap in A.
      int optimalPathCntr = 0;
      // Calculate recursively for matched state or diagonal alignment
      if(InsertInA>=Diago && InsertInA>=InsertInB){
        // given InsertInA in the last alignment and a diagobal alignment in current step, the traceback will be DA.
        affineAlignObj.Traceback[Traceback_M_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DA; // DA: Diagonal TrA
        affineAlignObj.M[i*(signalB_len+1)+j] = InsertInA;
        optimalPathCntr += oPathsA[(i-1)*(signalB_len+1) + j-1];
      }
      if(InsertInB>=Diago && InsertInB>=InsertInA){
        // given InsertInB in the last alignment and a diagobal alignment in current step, the traceback will be DB.
        affineAlignObj.Traceback[Traceback_M_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DB; // DB: Diagonal TrB
        affineAlignObj.M[i*(signalB_len+1)+j] = InsertInB;
        optimalPathCntr += oPathsB[(i-1)*(signalB_len+1) + j-1];
      }
      if(Diago>=InsertInA && Diago>=InsertInB){
        // Given that signal Ai-1 is aligned to signal Bj-1, signal Ai is aligned to signal Bj. Hence Traceback_matrix_M = DM.
        affineAlignObj.Traceback[Traceback_M_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DM; // DM: Diagonal TrM
        affineAlignObj.M[i*(signalB_len+1)+j] = Diago;
        optimalPathCntr += oPathsM[(i-1)*(signalB_len+1) + j-1];
      }
      oPathsM[i*(signalB_len+1) + j] = optimalPathCntr;

      optimalPathCntr = 0;
      // Calculate recursively for insert in signalA. signalA is along rows, hence entries would be either TM, TA or TB.
      double AfromM = affineAlignObj.M[(i-1)*(signalB_len+1)+j] - go; // Signal Ai-1 is aligned to Bj. Ai is aligned to a gap. So gap opening penalty is subtracted.
      double AfromA = affineAlignObj.A[(i-1)*(signalB_len+1)+j] - ge; // Signal Ai-1 is already aligned to a gap. Ai is aligned to a gap. So gap extension penalty is subtracted.
      double AfromB = affineAlignObj.B[(i-1)*(signalB_len+1)+j] - go; // Signal Ai-1 (gap) is aligned to signal Bj. Ai is aligned to a gap in B. So a gap in B is introduced, thus, gap opening penalty is subtracted.
      if(AfromA >= AfromM && AfromA >= AfromB){
        // Given signal Ai is aligned to gap and signal Ai-1 is already aligned to gap, The way to traceback is TA (Top from A to A).
        affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TA; // TA: Top TrA
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromA;
        optimalPathCntr += oPathsA[(i-1)*(signalB_len+1) + j];
        }
      if(AfromB >= AfromM && AfromB >= AfromA){
        // Given signal Ai is aligned to gap and signal Ai-1 (gap) is aligned to Bj, The way to traceback is TB (Top from A to B).
        affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TB; // TB: Top TrB
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromB;
        optimalPathCntr += oPathsB[(i-1)*(signalB_len+1) + j];
      }
      if(AfromM >= AfromA && AfromM >= AfromB){
        // Given signal Ai is aligned to gap and signal Ai-1 is aligned to Bj, The way to traceback is TM (Top from A to M).
        affineAlignObj.Traceback[Traceback_A_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TM; // TM: Top TrM
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromM;
        optimalPathCntr += oPathsM[(i-1)*(signalB_len+1) + j];
      }
      oPathsA[i*(signalB_len+1) + j] = optimalPathCntr;

      optimalPathCntr = 0;
      // Calculate recursively for insert in signalB. signalB is along rows, hence entries would be either LM, LA or LB.
      double BfromM = affineAlignObj.M[i*(signalB_len+1)+j-1] - go; // Signal Ai is aligned to Bj-1. Because Bj is aligned to a gap, so gap opening penalty is subtracted.
      double BfromA = affineAlignObj.A[i*(signalB_len+1)+j-1] - go; // Signal Ai is aligned to signal Bj-1(gap). Because Bj is aligned to a gap, so a gap in A is introduced, thus, gap opening penalty is subtracted.
      double BfromB = affineAlignObj.B[i*(signalB_len+1)+j-1] - ge; // Signal Bj-1 is already aligned to a gap. Because  Bj is aligned to a gap also, so gap extension penalty is subtracted.
      if(BfromA >= BfromM && BfromA >= BfromB){
        // Given signal Bj is aligned to a gap and signal Bj-1 (gap) is aligned to signal Ai, The way to traceback is LA (Left from B to A).
        affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LA; // LA: Left TrA
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromA;
        optimalPathCntr += oPathsA[i*(signalB_len+1) + j-1];
        }
      if(BfromB >= BfromM && BfromB >= BfromA){
        // Given signal Bj is aligned to a gap and signal Bj-1 is already aligned to gap, The way to traceback is LB (Left from B to B).
        affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LB; // LB: Left TrB
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromB;
        optimalPathCntr += oPathsB[i*(signalB_len+1) + j-1];
      }
      if(BfromM >= BfromA && BfromM >= BfromB){
        // Given signal Bj is aligned to gap and signal Ai is aligned to Bj-1, The way to traceback is LM (Left from B to M).
        affineAlignObj.Traceback[Traceback_B_index*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LM; // LM: Left TrM
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromM;
        optimalPathCntr += oPathsM[i*(signalB_len+1) + j-1];
      }
      oPathsB[i*(signalB_len+1) + j] = optimalPathCntr;

      optimalPathCntr = 0;
      if(affineAlignObj.M[i*(signalB_len+1)+j]>=affineAlignObj.A[i*(signalB_len+1)+j] &&
         affineAlignObj.M[i*(signalB_len+1)+j]>=affineAlignObj.B[i*(signalB_len+1)+j]){
        optimalPathCntr += oPathsM[i*(signalB_len+1) + j];
      }
      if(affineAlignObj.A[i*(signalB_len+1)+j]>=affineAlignObj.M[i*(signalB_len+1)+j] &&
         affineAlignObj.A[i*(signalB_len+1)+j]>=affineAlignObj.B[i*(signalB_len+1)+j]){
        optimalPathCntr += oPathsA[i*(signalB_len+1) + j];
      }
      if(affineAlignObj.B[i*(signalB_len+1)+j]>=affineAlignObj.A[i*(signalB_len+1)+j] &&
         affineAlignObj.B[i*(signalB_len+1)+j]>=affineAlignObj.M[i*(signalB_len+1)+j]){
        optimalPathCntr += oPathsB[i*(signalB_len+1) + j];
      }
      affineAlignObj.optionalPaths[i*(signalB_len+1) + j] = optimalPathCntr;
      }
    }
  delete[] oPathsM;
  delete[] oPathsA;
  delete[] oPathsB;
}

void getAffineAlignedIndices(AffineAlignObj &affineAlignObj, int bandwidth){
  AlignedIndices alignedIdx; // initialize empty struct.
  TracebackType TracebackPointer;
  tbJump MatName; // Matrix name M = 0, A = 1 or B = 2
  double affineAlignmentScore;
  int ROW_IDX = affineAlignObj.signalA_len;
  int COL_IDX = affineAlignObj.signalB_len;
  int ROW_SIZE = (affineAlignObj.signalA_len)+1;
  int COL_SIZE = (affineAlignObj.signalB_len)+1;

  if(affineAlignObj.FreeEndGaps == true){
    // Overlap Alignment
    // Maximum score and corresponding indices along the last column and last row is searched across all three matrices.
    // Matrix name and maximum score indices are passed by reference.
    affineAlignmentScore = getOlapAffineAlignStartIndices(affineAlignObj.M, affineAlignObj.A, affineAlignObj.B, ROW_SIZE, COL_SIZE, ROW_IDX, COL_IDX, MatName);
    if(ROW_IDX != affineAlignObj.signalA_len){
      // Maximum score is obtained in last column. Align all row indices below max-score-index to NA.
      for (int i = affineAlignObj.signalA_len; i>ROW_IDX; i--){
        alignedIdx.indexA_aligned.push_back(i);
        alignedIdx.indexB_aligned.push_back(NA); // Insert NA in signalB.
        alignedIdx.score.push_back(affineAlignmentScore); // Insert maxScore instead of score from the matrix M.
        //affineAlignObj.Path[i*COL_SIZE+COL_IDX] = true;
        //fillSimPath(affineAlignObj.simPath, bandwidth, i, COL_IDX, ROW_SIZE, COL_SIZE);
        }
      }
    else if (COL_IDX != affineAlignObj.signalB_len){
      // Maximum score is obtained in last row. Align all column indices right to max-score-index to NA.
      for (int j = affineAlignObj.signalB_len; j>COL_IDX; j--){
        alignedIdx.indexA_aligned.push_back(NA); // Insert NA in signalA.
        alignedIdx.indexB_aligned.push_back(j);
        alignedIdx.score.push_back(affineAlignmentScore); // Insert maxScore instead of score from the matrix M.
        //affineAlignObj.Path[ROW_IDX*COL_SIZE+j] = true;
        //fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, j, ROW_SIZE, COL_SIZE);
        }
      }
    }
  else {
    // Global Alignment, traceback starts at the bottom-right corner.
    // Search for the matrix which has the highest score at bottom-right corner.
    double Mscore = affineAlignObj.M[ROW_IDX*COL_SIZE+COL_IDX];
    double Ascore = affineAlignObj.A[ROW_IDX*COL_SIZE+COL_IDX];
    double Bscore = affineAlignObj.B[ROW_IDX*COL_SIZE+COL_IDX];
    if (Mscore >= Ascore && Mscore >= Bscore){
      affineAlignmentScore = Mscore;
      MatName = M;
      }
    else if(Ascore >= Mscore && Ascore>= Bscore) {
      affineAlignmentScore = Ascore;
      MatName = A;
      }
    else {
      affineAlignmentScore = Bscore;
      MatName = B;
      }
    }

  alignedIdx.score.push_back(affineAlignmentScore);
  TracebackPointer = affineAlignObj.Traceback[MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX];
  affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
  fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
  // Traceback path and align row indices to column indices.

  while(TracebackPointer != SS){
    // SS: STOP when top-left corner of the matrix is reached
    // D: Diagonal, T: Top, L: Left
    // Rcpp::Rcout << TracebackPointer << std::endl;
    switch(TracebackPointer){
    // In the code below, we are appending future scores. Because, once we go to the M,A or B matrix.
    // we will not be able to tell which matrix we are currently in.
    case DM:
      {// Go diagonal (Up-Left) to the matrix M.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      ROW_IDX = ROW_IDX-1;
      COL_IDX = COL_IDX-1;
      MatName = M;
      alignedIdx.score.push_back(affineAlignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      break;}

    case DA:
      {// Go diagonal (Up-Left) to the matrix A.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      ROW_IDX = ROW_IDX-1;
      COL_IDX = COL_IDX-1;
      MatName = A;
      alignedIdx.score.push_back(affineAlignObj.A[ROW_IDX*COL_SIZE+COL_IDX]);
      affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      break;}

    case DB:
      {// Go diagonal (Up-Left) to the matrix B.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      ROW_IDX = ROW_IDX-1;
      COL_IDX = COL_IDX-1;
      MatName = B;
      alignedIdx.score.push_back(affineAlignObj.B[ROW_IDX*COL_SIZE+COL_IDX]);
      affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
      fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      break;}

    case TM:
      {// Go up to the matrix M.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(NA);
      ROW_IDX = ROW_IDX-1;
      MatName = M;
      alignedIdx.score.push_back(affineAlignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      if(COL_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;}

    case TA:
      {// Go up to the matrix A.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(NA);
      ROW_IDX = ROW_IDX-1;
      MatName = A;
      alignedIdx.score.push_back(affineAlignObj.A[ROW_IDX*COL_SIZE+COL_IDX]);
      if(COL_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;}

    case TB:
      {// Go up to the matrix B.
      alignedIdx.indexA_aligned.push_back(ROW_IDX);
      alignedIdx.indexB_aligned.push_back(NA);
      ROW_IDX = ROW_IDX-1;
      MatName = B;
      alignedIdx.score.push_back(affineAlignObj.B[ROW_IDX*COL_SIZE+COL_IDX]);
      if(COL_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;}

    case LM:
      {// Go left to the matrix M.
      alignedIdx.indexA_aligned.push_back(NA);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      COL_IDX = COL_IDX-1;
      MatName = M;
      alignedIdx.score.push_back(affineAlignObj.M[ROW_IDX*COL_SIZE+COL_IDX]);
      if(ROW_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;
      }

    case LA:
      {// Go left to the matrix A.
      alignedIdx.indexA_aligned.push_back(NA);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      COL_IDX = COL_IDX-1;
      MatName = A;
      alignedIdx.score.push_back(affineAlignObj.A[ROW_IDX*COL_SIZE+COL_IDX]);
      if(ROW_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;}

    case LB:
      {// Go left to the matrix B.
      alignedIdx.indexA_aligned.push_back(NA);
      alignedIdx.indexB_aligned.push_back(COL_IDX);
      COL_IDX = COL_IDX-1;
      MatName = B;
      alignedIdx.score.push_back(affineAlignObj.B[ROW_IDX*COL_SIZE+COL_IDX]);
      if(ROW_IDX != 0){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      else if(!affineAlignObj.FreeEndGaps){
        affineAlignObj.nGaps += 1;
        affineAlignObj.Path[ROW_IDX*COL_SIZE+COL_IDX] = true;
        fillSimPath(affineAlignObj.simPath, bandwidth, ROW_IDX, COL_IDX, ROW_SIZE, COL_SIZE);
      }
      break;}

    }
    // Read traceback for the next iteration.
    TracebackPointer = affineAlignObj.Traceback[MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX];
  }
  // push_back adds values at the end of vector, therefore, reverse the vector.
  std::reverse(std::begin(alignedIdx.indexA_aligned), std::end(alignedIdx.indexA_aligned));
  std::reverse(std::begin(alignedIdx.indexB_aligned), std::end(alignedIdx.indexB_aligned));
  std::reverse(std::begin(alignedIdx.score), std::end(alignedIdx.score));
  // remove the first index, since the score-traceback is ahead of aligned indices.
  alignedIdx.score.erase(alignedIdx.score.begin());
  // Copy aligned indices to alignObj.
  affineAlignObj.indexA_aligned = alignedIdx.indexA_aligned;
  affineAlignObj.indexB_aligned = alignedIdx.indexB_aligned;
  affineAlignObj.score = alignedIdx.score;
  //return;
}

// It finds start indices and matrix for tracebackin in case of overlap affine alignment.
double getOlapAffineAlignStartIndices(double* MatrixM, double* MatrixA, double* MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName){
  double maxScore = -std::numeric_limits<double>::infinity(); //- Inf
  int MaxRowIndex, MaxColIndex;
  for(int i = 0; i < ROW_SIZE; i++){
    if(MatrixM[i*COL_SIZE+COL_SIZE-1] >= maxScore){
      // Search the maximum score along the last column of M matrix
      MaxRowIndex = i;
      MaxColIndex = COL_SIZE-1;
      MatrixName = M;
      maxScore = MatrixM[i*COL_SIZE+COL_SIZE-1];
    }
    else if(MatrixA[i*COL_SIZE+COL_SIZE-1] >= maxScore){
      // Search the maximum score along the last column of A matrix
      MaxRowIndex = i;
      MaxColIndex = COL_SIZE-1;
      MatrixName = A;
      maxScore = MatrixA[i*COL_SIZE+COL_SIZE-1];
      }
    else if(MatrixB[i*COL_SIZE+COL_SIZE-1] >= maxScore){
      // Search the maximum score along the last column of B matrix
      MaxRowIndex = i;
      MaxColIndex = COL_SIZE-1;
      MatrixName = B;
      maxScore = MatrixB[i*COL_SIZE+COL_SIZE-1];
      }
    }
  for (int j = 0; j < COL_SIZE; j++){
    if(MatrixM[(ROW_SIZE-1)*COL_SIZE+j] >= maxScore){
      // Search the maximum score along the last row of M matrix
      MaxRowIndex = ROW_SIZE-1;
      MaxColIndex = j;
      MatrixName = M;
      maxScore = MatrixM[(ROW_SIZE-1)*COL_SIZE+j];
      }
    else if(MatrixA[(ROW_SIZE-1)*COL_SIZE+j] >= maxScore){
      // Search the maximum score along the last row of A matrix
      MaxRowIndex = ROW_SIZE-1;
      MaxColIndex = j;
      MatrixName = A;
      maxScore = MatrixA[(ROW_SIZE-1)*COL_SIZE+j];
      }
    else if(MatrixB[(ROW_SIZE-1)*COL_SIZE+j] >= maxScore){
      // Search the maximum score along the last row of B matrix
      MaxRowIndex = ROW_SIZE-1;
      MaxColIndex = j;
      MatrixName = B;
      maxScore = MatrixB[(ROW_SIZE-1)*COL_SIZE+j];
      }
    }
  // Copy max-score indices to pass-by-reference variables
  OlapStartRow = MaxRowIndex;
  OlapStartCol = MaxColIndex;
  return maxScore;
}

void fillSimPath(bool* simPath, int bandwidth, int ROW_IDX, int COL_IDX, int ROW_SIZE, int COL_SIZE){
  for (int i = ROW_IDX-bandwidth; i<=ROW_IDX+bandwidth; i++){
    if(i>=0 && i<ROW_SIZE){
      simPath[i*COL_SIZE+COL_IDX] = true;
    }
  }
  for (int j = COL_IDX-bandwidth; j<=COL_IDX+bandwidth; j++){
    if(j>=0 && j<COL_SIZE){
      simPath[ROW_IDX*COL_SIZE+j] = true;
    }
  }
}

double getForwardSim(const SimMatrix& s, bool* simPath){
  double forwardSim = 0;
  int COL_SIZE = s.n_col+1;
  for(int i=0; i<s.n_row; i++){
    for(int j=0; j<s.n_col; j++){
      if(simPath[(i+1)*COL_SIZE+(j+1)]){
        forwardSim += s.data[i*s.n_col + j];
      }
    }
  }
  return forwardSim;
}

/***
 * https://www.coursera.org/lecture/bioinformatics-pku/alignment-with-affine-gap-penalty-and-calculation-of-time-complexity-of-the-0ya7X
 * https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
 * Biological Sequence Analysis (Chapter 2) by Durbin, Eddy, Krogh, and Mitchison
 *
 ***/

} // namespace DIAlign
