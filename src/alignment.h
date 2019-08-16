#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <limits>

#include "utils.h"
#include "similarityMatrix.h"
#include "affinealignobj.h"

namespace DIAlign
{

const int NA = 0;

struct AlignedIndices{
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<double> score;
};

struct AlignObj
{
    std::vector<double> s_data;
    std::vector<TracebackType> Traceback;
    std::vector<double> M;
    std::vector<double> M_forw;
    std::vector<bool> Path; // Path matrix would represent alignment path through similarity matrix as binary-hot encoding.
    std::vector<bool> simPath;
    std::vector<int> OptionalPaths; // Highlight the number of all optimal paths.
    int signalA_len;
    int signalB_len;
    double GapOpen;
    double GapExten;
    bool FreeEndGaps;
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<double> score;
    double score_forw;
    double simScore_forw; // Summation of similarity score along the alignment path-band.
    double alterAlignScore; // Alignment score of 2nd best peak alignment.
    int nGaps;

    // Not a default constructor
    AlignObj(int ROW_SIZE, int COL_SIZE)
    {
      M.resize(ROW_SIZE * COL_SIZE, 0);
      M_forw.resize(ROW_SIZE * COL_SIZE, 0);
      Traceback.resize(ROW_SIZE * COL_SIZE, SS);
      Path.resize(ROW_SIZE * COL_SIZE, false);
      simPath.resize(ROW_SIZE * COL_SIZE, false);
      OptionalPaths.resize(ROW_SIZE * COL_SIZE, 0);
      signalA_len = ROW_SIZE-1;
      signalB_len = COL_SIZE-1;
      s_data.resize(signalA_len * signalB_len, 0.0);
      GapOpen = 0.0;
      GapExten = 0.0;
      FreeEndGaps = true;
      score_forw = 0.0;
      simScore_forw = 0.0;
      alterAlignScore = 0.0;
      nGaps = 0;
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

void fillsimPath(std::vector<bool> &simPath, int bandwidth, int ROW_IDX, int COL_IDX, int ROW_SIZE, int COL_SIZE);

}

#endif // ALIGNMENT_H
