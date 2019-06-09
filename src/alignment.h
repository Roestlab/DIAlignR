#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <limits>

#include "simpleFcn.h"
#include "affinealignobj.h"
#define NA 0

struct AlignedIndices{
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<double> score;
};

struct AlignObj
{
    std::vector<TracebackType> Traceback;
    std::vector<double> M;
    int signalA_len;
    int signalB_len;
    double GapOpen;
    double GapExten;
    bool FreeEndGaps;
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<double> score;

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

// This function performs dynamic programming and calculates "M" and "Traceback". Traceback matrix keeps record of the path as we fill matrix M.
AlignObj doAlignment(SimMatrix s, double gap, bool OverlapAlignment);

// This tracebacks along the highest scoring path, preparing list of scores and aligned indices.
void getAlignedIndices(AlignObj &alignObj);

#endif // ALIGNMENT_H
