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

    // Rule 1 Copy constructor
    AlignObj(const AlignObj &other)
    {
      signalA_len = other.signalA_len;
      signalB_len = other.signalB_len;
      GapOpen = other.GapOpen;
      GapExten = other.GapExten;
      FreeEndGaps = other.FreeEndGaps;
    }

    // Rule 2 Copy assignment operator
    AlignObj& operator=(const AlignObj& other)
    {
      if(this == &other) return *this; // handling of self assignment.
      signalA_len = other.signalA_len;
      signalB_len = other.signalB_len;
      GapOpen = other.GapOpen;
      GapExten = other.GapExten;
      FreeEndGaps = other.FreeEndGaps;
      return *this;
    }

    // Rule 3 Not a default destructor
    ~AlignObj()
    {  }
};

AlignObj doAlignment(NumericMatrix s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment);

#endif // ALIGNMENT_H
