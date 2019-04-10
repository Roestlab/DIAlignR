#include "affinealignment.h"
// Do not inclue cpp file because compiler will build the Obj through two different path.


AffineAlignObj doAffineAlignment(NumericMatrix s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment){
  AffineAlignObj affineAlignObj(signalA_len+1, signalB_len+1);
  affineAlignObj.FreeEndGaps = OverlapAlignment;
  affineAlignObj.GapOpen = go;
  affineAlignObj.GapExten = ge;

  // Initialize first row and first column for global and overlap alignment.
  float Inf = std::numeric_limits<float>::infinity();
  for(int i = 0; i<=signalA_len; i++){
    affineAlignObj.M[i*(signalB_len+1)+0] = -Inf;
    affineAlignObj.B[i*(signalB_len+1)+0] = -Inf;
    affineAlignObj.Traceback[0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = SS; //STOP
    affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1)) + i*(signalB_len+1)+0] = SS; //STOP
    }

  for(int j = 0; j<=signalB_len; j++){
    affineAlignObj.M[0*(signalB_len+1)+j] = -Inf;
    affineAlignObj.A[0*(signalB_len+1)+j] = -Inf;
    affineAlignObj.Traceback[0*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = SS; //STOP
    affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = SS; //STOP
    }
  affineAlignObj.M[0*(signalB_len+1)+0] = 0;

  if(affineAlignObj.FreeEndGaps == true){
    for(int i = 1; i<=signalA_len; i++){
      affineAlignObj.A[i*(signalB_len+1) + 0] = 0;
      affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = TA; //TOP A
      }
    for(int j = 1; j<=signalB_len; j++){
      affineAlignObj.B[0*(signalB_len+1)+j] = 0;
      affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = LB; //LEFT B
      }
    }
  else {
      for(int i = 1; i<=signalA_len; i++){
        affineAlignObj.A[i*(signalB_len+1) + 0] = -(i-1)*ge - go;
        affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0] = TA; //TOP A
        }
      for(int j = 1; j<=signalB_len; j++){
        affineAlignObj.B[0*(signalB_len+1)+j] = -(j-1)*ge - go;
        affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j] = LB; //LEFT B
        }
      }

  // Perform dynamic programming for affine alignment
  float Diago, gapInA, gapInB;
  for(int i=1; i<=signalA_len; i++){
    for(int j=1; j<=signalB_len; j++){
      float sI_1J_1 = s(i-1, j-1);
      Diago = affineAlignObj.M[(i-1)*(signalB_len+1)+j-1] + sI_1J_1;
      gapInA = affineAlignObj.A[(i-1)*(signalB_len+1)+j-1] + sI_1J_1;
      gapInB = affineAlignObj.B[(i-1)*(signalB_len+1)+j-1] + sI_1J_1;

      // Calculate recursively for matched alignment
      if(Diago>=gapInA && Diago>=gapInB){
        affineAlignObj.Traceback[0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DM; // DM: Diagonal TrM
        affineAlignObj.M[i*(signalB_len+1)+j] = Diago;
        }
      else if (gapInA>=Diago && gapInA>=gapInB){
        affineAlignObj.Traceback[0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DA; // DA: Diagonal TrA
        affineAlignObj.M[i*(signalB_len+1)+j] = gapInA;
        }
      else{
        affineAlignObj.Traceback[0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = DB; // DB: Diagonal TrB
        affineAlignObj.M[i*(signalB_len+1)+j] = gapInB;
        }

      // Calculate recursively for gap in signalB
      float AfromM = affineAlignObj.M[(i-1)*(signalB_len+1)+j] - go;
      float AfromA = affineAlignObj.A[(i-1)*(signalB_len+1)+j] - ge;
      float AfromB = affineAlignObj.B[(i-1)*(signalB_len+1)+j] - go;
      if(AfromM >= AfromA && AfromM >= AfromB){
        affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TM; // TM: Top TrM
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromM;
        }
      else if (AfromA >= AfromM && AfromA >= AfromB){
        affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TA; // TA: Top TrA
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromA;
        }
      else{
        affineAlignObj.Traceback[1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = TB; // TB: Top TrB
        affineAlignObj.A[i*(signalB_len+1)+j] = AfromB;
        }

      // Calculate recursively for gap in signalA
      float BfromM = affineAlignObj.M[i*(signalB_len+1)+j-1] - go;
      float BfromA = affineAlignObj.A[i*(signalB_len+1)+j-1] - go;
      float BfromB = affineAlignObj.B[i*(signalB_len+1)+j-1] - ge;
      if(BfromM >= BfromA && BfromM >= BfromB){
        affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LM; // LM: Left TrM
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromM;
        }
      else if (BfromA >= BfromM && BfromA >= BfromB){
        affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LA; // LA: Left TrA
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromA;
        }
      else{
        affineAlignObj.Traceback[2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j] = LB; // LB: Left TrB
        affineAlignObj.B[i*(signalB_len+1)+j] = BfromB;
        }
      }
    }
  printMatrix(affineAlignObj.Traceback, signalA_len+1, signalB_len+1);
  // printMatrix(affineAlignObj.Traceback[1*(signalA_len+1)*(signalB_len+1)], signalA_len+1, signalB_len+1);
  // printMatrix(affineAlignObj.Traceback[2*(signalA_len+1)*(signalB_len+1)], signalA_len+1, signalB_len+1);
  return affineAlignObj;
  }

/***
AlignedIndices getAffineAlignedIndices(AffineAlignObj &affineAlignObj){
    AlignedIndices alignedIdx;
    char TracebackPointer;
    tbJump MatName;
    float affineAlignmentScore;
    int ROW_IDX = affineAlignObj.signalA_len;
    int COL_IDX = affineAlignObj.signalB_len;
    int ROW_SIZE = (affineAlignObj.signalA_len)+1;
    int COL_SIZE = (affineAlignObj.signalB_len)+1;

    if(affineAlignObj.FreeEndGaps == true){
        // Overlap Alignment
        affineAlignmentScore = getOlapAffineAlignStartIndices((float *)affineAlignObj.M, (float *)affineAlignObj.A, (float *)affineAlignObj.B, ROW_SIZE, COL_SIZE, ROW_IDX, COL_IDX, MatName);
        if(ROW_IDX != affineAlignObj.signalA_len){
            for (int i = affineAlignObj.signalA_len; i>ROW_IDX; i--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), i);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
                alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
            }
        } else if (COL_IDX != affineAlignObj.signalB_len){
            for (int j = affineAlignObj.signalB_len; j>COL_IDX; j--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), j);
                alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
            }
        }
    } else {
        // Global Alignment
        float Mscore = *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX);
        float Ascore = *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX);
        float Bscore = *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX);
        if (Mscore >= Ascore && Mscore >= Bscore){
            affineAlignmentScore = Mscore;
            MatName = M;
        } else if(Ascore >= Mscore && Ascore>= Bscore) {
            affineAlignmentScore = Ascore;
            MatName = A;
        } else{
            affineAlignmentScore = Bscore;
            MatName = B;
        }
    }

//    std::cout << "************************" << std::endl;
//    for (std::vector<float>::iterator it = alignedIdx.score.begin(); it != alignedIdx.score.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    for (std::vector<int>::iterator it = alignedIdx.indexA_aligned.begin(); it != alignedIdx.indexA_aligned.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    for (std::vector<int>::iterator it = alignedIdx.indexB_aligned.begin(); it != alignedIdx.indexB_aligned.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    std::cout << "************************" << std::endl;

    alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
    alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
    alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
    TracebackPointer = *(affineAlignObj.Traceback+MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX);
    // Traceback path and align row indices to column indices.
    while(TracebackPointer != SS){
        // D: Diagonal, T: Top, L: Left
        switch(TracebackPointer){
        // In the code below, we are appending future values. Because, once we go to the M,A or B matrix.
        // we will not be able to tell which matrix we are currently in.
        case DM: {
            ROW_IDX = ROW_IDX-1;
            MatName = M;
            COL_IDX = COL_IDX-1;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
            break;
                }
        case DA:
        {
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case DB:
        {
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case TM:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = M;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case TA:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case TB:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case LM:
        {
            COL_IDX = COL_IDX-1;
            MatName = M;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case LA:
        {
            COL_IDX = COL_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
        case LB:
        {
            COL_IDX = COL_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
            break;
        }
            }

        TracebackPointer = *(affineAlignObj.Traceback+MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX);
        }

//    std::cout << "************************" << std::endl;
//    for (std::vector<float>::iterator it = alignedIdx.score.begin(); it != alignedIdx.score.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    for (std::vector<int>::iterator it = alignedIdx.indexA_aligned.begin(); it != alignedIdx.indexA_aligned.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    for (std::vector<int>::iterator it = alignedIdx.indexB_aligned.begin(); it != alignedIdx.indexB_aligned.end(); it++)
//        std::cout << *it << " ";
//    std::cout << std::endl;
//    std::cout << "************************" << std::endl;

    alignedIdx.score.erase(alignedIdx.score.begin());
    alignedIdx.indexA_aligned.erase(alignedIdx.indexA_aligned.begin());
    alignedIdx.indexB_aligned.erase(alignedIdx.indexB_aligned.begin());
    return alignedIdx;
}

template<class T>
T getOlapAffineAlignStartIndices(T *MatrixM, T *MatrixA, T *MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName){
    T affineAlignmentScore;
    float maxScore = -std::numeric_limits<float>::infinity();
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
        if(*(MatrixM+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = M;
            maxScore = *(MatrixM+i*COL_SIZE+COL_SIZE-1);
        } else if(*(MatrixA+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = A;
            maxScore = *(MatrixA+i*COL_SIZE+COL_SIZE-1);
        } else if(*(MatrixB+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = B;
            maxScore = *(MatrixB+i*COL_SIZE+COL_SIZE-1);
        }
    }
    for (int j = 0; j < COL_SIZE; j++){
        if(*(MatrixM+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = M;
            maxScore = *(MatrixM+(ROW_SIZE-1)*COL_SIZE+j);
        } else if(*(MatrixA+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = A;
            maxScore = *(MatrixA+(ROW_SIZE-1)*COL_SIZE+j);
        } else if(*(MatrixB+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = B;
            maxScore = *(MatrixB+(ROW_SIZE-1)*COL_SIZE+j);
        }
    }
    OlapStartRow = MaxRowIndex;
    OlapStartCol = MaxColIndex;
    affineAlignmentScore = maxScore;
    return affineAlignmentScore;
}
***/
