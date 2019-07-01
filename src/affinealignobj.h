#ifndef AFFINEALIGNOBJ_H
#define AFFINEALIGNOBJ_H

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

namespace DIAlign
{
enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};
enum tbJump {M = 0, A = 1, B = 2};

// This function overloads << to display TracebackType.
std::ostream& operator<<(std::ostream& out, const TracebackType value);

// This function converts TracebackType Enum to characters.
std::vector<char> EnumToChar(std::vector<TracebackType> v);

struct AffineAlignObj
{
private:
  int signalA_capacity; // Capacity of matrix for signal A (rows)
  int signalB_capacity; // Capacity of matrix for signal B (columns)

public:
  double* M; // Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.
  double* A; // Insert in sequence A, residue in A is aligned to gap in B. A(i,j) is the best score given that Ai is aligned to a gap in B.
  double* B; // Insert in sequence B, residue in B is aligned to gap in A. B(i,j) is the best score given that Bj is aligned to a gap in A.
  TracebackType* Traceback;
  bool* Path; // Path matrix would represent alignment path through similarity matrix as binary-hot encoding.

  int signalA_len; // Number of data-points in signal A
  int signalB_len; // Number of data-points in signal B
  double GapOpen; // Penalty for Gap opening
  double GapExten; // Penalty for Gap extension
  // For single gap: Penalty = GapOpen
  // For two consecutive gaps: Penalty = GapOpen + GapExten
  // For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten
  bool FreeEndGaps; // True for Overlap alignment
  std::vector<int> indexA_aligned; // Aligned signalA indices after affine alignment
  std::vector<int> indexB_aligned; // Aligned signalB indices after affine alignment
  std::vector<double> score;  // Score along the aligned path

  // Not a default constructor
  AffineAlignObj() {}
  AffineAlignObj(int ROW_SIZE, int COL_SIZE, bool clearMemory = true)
  {
    M = new double[ROW_SIZE * COL_SIZE];
    A = new double[ROW_SIZE * COL_SIZE];
    B = new double[ROW_SIZE * COL_SIZE];
    Traceback = new TracebackType[3* ROW_SIZE * COL_SIZE];
    Path = new bool[ROW_SIZE * COL_SIZE];

    if (clearMemory)
    {
      std::memset(M, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(A, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(B, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(Traceback, SS, 3 * ROW_SIZE * COL_SIZE * sizeof(TracebackType));
      std::memset(Path, 0, ROW_SIZE * COL_SIZE * sizeof(bool));
    }

    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;
    GapOpen = 0.0;
    GapExten = 0.0;
    FreeEndGaps = true;

    signalA_capacity = ROW_SIZE-1;
    signalB_capacity = COL_SIZE-1;
  }

  void reset(int ROW_SIZE, int COL_SIZE)
  {
    if (ROW_SIZE -1 > signalA_capacity || COL_SIZE -1 > signalB_capacity)
    {
      std::cout << "Error: cannot reset an object beyond capacity" << std::endl;
      std::cout << ROW_SIZE << " vs " << signalA_capacity << std::endl;
      throw 1;
    }

    std::memset(M, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(A, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(B, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(Traceback, SS, 3 * ROW_SIZE * COL_SIZE * sizeof(TracebackType));
    std::memset(Path, 0, ROW_SIZE * COL_SIZE * sizeof(bool));

    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;

    GapOpen = 0.0;
    GapExten = 0.0;
    FreeEndGaps = true;

    indexA_aligned.clear();
    indexB_aligned.clear();
    score.clear();
  }

  // Rule 2 Copy assignment operator
  AffineAlignObj& operator=(const AffineAlignObj& rhs)
  {
    delete[] M;
    delete[] A;
    delete[] B;
    delete[] Traceback;
    delete[] Path;

    //std::cout << " this " << this << std::endl;
    signalA_len = rhs.signalA_len;
    signalA_capacity = rhs.signalA_capacity;
    signalB_len = rhs.signalB_len;
    signalB_capacity = rhs.signalB_capacity;

    GapOpen = rhs.GapOpen;
    GapExten = rhs.GapExten;
    FreeEndGaps = rhs.FreeEndGaps;
    indexA_aligned = rhs.indexA_aligned;
    indexB_aligned = rhs.indexB_aligned;
    score = rhs.score;

    int ROW_SIZE = rhs.signalA_len + 1;
    int COL_SIZE = rhs.signalB_len + 1;

    M = new double[ROW_SIZE * COL_SIZE];
    A = new double[ROW_SIZE * COL_SIZE];
    B = new double[ROW_SIZE * COL_SIZE];
    Traceback = new TracebackType[3* ROW_SIZE * COL_SIZE];
    Path = new bool[ROW_SIZE * COL_SIZE];

    std::memcpy(M, rhs.M, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(A, rhs.A, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(B, rhs.B, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(Traceback, rhs.Traceback, 3 *ROW_SIZE * COL_SIZE * sizeof(TracebackType));
    std::memcpy(Path, rhs.Path, ROW_SIZE * COL_SIZE * sizeof(bool));
  }

  // Rule 1 Copy constructor
  AffineAlignObj(const AffineAlignObj& rhs)
  {
    delete[] M;
    delete[] A;
    delete[] B;
    delete[] Traceback;
    delete[] Path;

    signalA_len = rhs.signalA_len;
    signalA_capacity = rhs.signalA_capacity;
    signalB_len = rhs.signalB_len;
    signalB_capacity = rhs.signalB_capacity;

    GapOpen = rhs.GapOpen;
    GapExten = rhs.GapExten;
    FreeEndGaps = rhs.FreeEndGaps;
    indexA_aligned = rhs.indexA_aligned;
    indexB_aligned = rhs.indexB_aligned;
    score = rhs.score;

    int ROW_SIZE = rhs.signalA_len + 1;
    int COL_SIZE = rhs.signalB_len + 1;

    M = new double[ROW_SIZE * COL_SIZE];
    A = new double[ROW_SIZE * COL_SIZE];
    B = new double[ROW_SIZE * COL_SIZE];
    Traceback = new TracebackType[3* ROW_SIZE * COL_SIZE];
    Path = new bool[ROW_SIZE * COL_SIZE];

    std::memcpy(M, rhs.M, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(A, rhs.A, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(B, rhs.B, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(Traceback, rhs.Traceback, 3 *ROW_SIZE * COL_SIZE * sizeof(TracebackType));
    std::memcpy(Path, rhs.Path, ROW_SIZE * COL_SIZE * sizeof(bool));
  }

  // Rule 3 Not a default destructor
  ~AffineAlignObj()
  {
    delete[] M;
    delete[] A;
    delete[] B;
    delete[] Traceback;
    delete[] Path;
  }

};
} // namespace DIAlign

#endif // AFFINEALIGNOBJ_H

// brackets ([]) specify the index of an element of the array. In fact these brackets are a dereferencing operator known as offset operator.
// We used bracked to work with STL containers. In case of arrays/ pointers, they work in the same way.
// double* M;
// M[i] == *(M+i);
