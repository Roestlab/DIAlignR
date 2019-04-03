#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <stdio.h>
#include <vector>

#include "simpleFcn.h"
#include "affinealignobj.h"
#define NA 0

struct AlignedIndices{
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<float> score;
};

struct AlignObj
{
    std::vector<TracebackType> Traceback;
    std::vector<float> M;
    int signalA_len;
    int signalB_len;
    float GapOpen;
    float GapExten;
    bool FreeEndGaps;

    // Not a default constructor
    AlignObj(int ROW_SIZE, int COL_SIZE)
    {
      M.resize(ROW_SIZE * COL_SIZE, 0);
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
    ~AlignObj()
    {  }
};

AlignObj doAlignment(NumericMatrix s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment);

AlignedIndices getAlignedIndices(AlignObj &alignObj);

#endif // ALIGNMENT_H
