#ifndef AFFINEALIGNOBJ_H
#define AFFINEALIGNOBJ_H

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};

std::ostream& operator<<(std::ostream& out, const TracebackType value);

enum tbJump {M = 0, A = 1, B = 2};

std::vector<char> EnumToChar(std::vector<TracebackType> v);

struct AffineAlignObj
{
  std::vector<float> M;
  std::vector<float> A;
  std::vector<float> B;
  std::vector<TracebackType> Traceback;
  int signalA_len; // stack allocation
  int signalB_len;
  float GapOpen;
  float GapExten;
  bool FreeEndGaps;
  std::vector<int> indexA_aligned;
  std::vector<int> indexB_aligned;
  std::vector<float> score;

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
