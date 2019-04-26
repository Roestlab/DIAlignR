#ifndef AFFINEALIGNOBJ_H
#define AFFINEALIGNOBJ_H

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};
enum tbJump {M = 0, A = 1, B = 2};

// This function overloads << to display TracebackType.
std::ostream& operator<<(std::ostream& out, const TracebackType value);

// This function converts TracebackType Enum to characters.
std::vector<char> EnumToChar(std::vector<TracebackType> v);

struct AffineAlignObj
{
  std::vector<float> M; // Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.
  std::vector<float> A; // Insert in sequence A, residue in A is aligned to gap in B. A(i,j) is the best score given that Ai is aligned to a gap in B.
  std::vector<float> B; // Insert in sequence B, residue in B is aligned to gap in A. B(i,j) is the best score given that Bj is aligned to a gap in A.
  std::vector<TracebackType> Traceback;
  int signalA_len; // Number of data-points in signal A
  int signalB_len; // Number of data-points in signal B
  float GapOpen; // Penalty for Gap opening
  float GapExten; // Penalty for Gap extension
  // For single gap: Penalty = GapOpen
  // For two consecutive gaps: Penalty = GapOpen + GapExten
  // For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten
  bool FreeEndGaps; // True for Overlap alignment
  std::vector<int> indexA_aligned; // Aligned signalA indices after affine alignment
  std::vector<int> indexB_aligned; // Aligned signalB indices after affine alignment
  std::vector<float> score;  // Score along the aligned path

  // Not a default constructor
  AffineAlignObj(int ROW_SIZE, int COL_SIZE)
  {
    M.resize(ROW_SIZE * COL_SIZE, 0);
    A.resize(ROW_SIZE * COL_SIZE, 0);
    B.resize(ROW_SIZE * COL_SIZE, 0);
    Traceback.resize(3 * ROW_SIZE * COL_SIZE, SS);
    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;
    GapOpen = 0.0;
    GapExten = 0.0;
    FreeEndGaps = true;
  }

  // Rule 1 Copy constructor and Rule 2 Copy assignment operator are not needed.
  // C++ auto-create them and will do the right-thing compared to manual code.

  // Rule 3 Not a default destructor
  ~AffineAlignObj()
  {  }
};

#endif // AFFINEALIGNOBJ_H
