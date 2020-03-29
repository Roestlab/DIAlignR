#ifndef AFFINEALIGNOBJ_H
#define AFFINEALIGNOBJ_H

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

/**
 * @namespace DIAlign
 * @brief Generic namespace for all classes and functions of DIAlign
 *
 */
namespace DIAlign
{

/**
 * @brief Generic namespace for Traceback that includes source matrix name and the arrow direction.
 *
 */
namespace Traceback
{
/// All combinations of cumulative score-matrix name and arrow directions.
enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};

/// cumulative score-matrix names
enum tbJump {M = 0, A = 1, B = 2};

/// This function overloads << to display TracebackType.
std::ostream& operator<<(std::ostream& out, const TracebackType value);

/// This function converts TracebackType Enum to characters.
std::vector<char> EnumToChar(std::vector<TracebackType> v);
}

/**
 * @brief An affine alignment object.
 *
 * This object contains similarity matrix, three matrices M, A and B storing cumulative-scores for dynamic programming.
 * Traceback matrices store source matrix name and direction as matrices are filled with dynamic programming.
 * Path matrix encode alignment path that results in the highest cumulative score.
 * The aligned indices are also stored for signal A and signal B.
 */
struct AffineAlignObj
{
private:
  int signalA_capacity; ///< Capacity of matrix for signal A (rows).
  int signalB_capacity; ///< Capacity of matrix for signal B (columns).

public:
  double* s_data; ///< similarity score matrix.
  double* M; ///< Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.
  double* A; ///< Insert in sequence A, residue in A is aligned to gap in B. A(i,j) is the best score given that Ai is aligned to a gap in B.
  double* B; ///< Insert in sequence B, residue in B is aligned to gap in A. B(i,j) is the best score given that Bj is aligned to a gap in A.
  Traceback::TracebackType* Traceback; ///< Traceback matrices store source matrix name and direction as matrices are filled with dynamic programming.
  bool* Path; ///< Path matrix would represent alignment path through similarity matrix as binary-hot encoding.
  bool* simPath; ///< Not needed, will be removed.
  // s_data, M, A and B should be private. Now there is a possibility of memory-leak.
  // TODO Make above variables private.
  int signalA_len; ///< Number of data-points in signal A.
  int signalB_len; ///< Number of data-points in signal B.
  double GapOpen; ///< Penalty for Gap opening. For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten.
  double GapExten; ///< Penalty for Gap extension. For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten.
  bool FreeEndGaps; ///< True for Overlap alignment.
  std::vector<int> indexA_aligned; ///< Aligned signalA indices after affine alignment.
  std::vector<int> indexB_aligned; ///< Aligned signalB indices after affine alignment.
  std::vector<double> score;  ///< Cumulative score along the aligned path.
  int nGaps; ///< Total number of gaps in the alignment path.

  /**
   * @brief Constructor for AffineAlignObj.
   *
   * Allocates memory for s_data, M, A, B, Traceback, Path matrices. Initialize them with zero if clearMemory is set true.
   * @param ROW_SIZE Number of rows in matrix M.
   * @param COL_SIZE Number of columns in matrix M.
   * @param clearMemory If true, matrices are initialized with zero.
   */
  // Not a default constructor
  AffineAlignObj(int ROW_SIZE, int COL_SIZE, bool clearMemory = true)
  {
    allocateMemory_(ROW_SIZE, COL_SIZE);

    // clearMemory means having default zero values.
    // We could use a for-loop but memset is faster for contiguous location in memory.
    // It makes byte value = unsigned(int_0)
    if (clearMemory)
    {
      std::memset(s_data, 0, (ROW_SIZE -1) * (COL_SIZE-1) * sizeof(double));
      std::memset(M, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(A, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(B, 0, ROW_SIZE * COL_SIZE * sizeof(double));
      std::memset(Traceback, Traceback::SS, 3 * ROW_SIZE * COL_SIZE * sizeof(Traceback::TracebackType));
      std::memset(Path, 0, ROW_SIZE * COL_SIZE * sizeof(bool));
      std::memset(simPath, 0, ROW_SIZE * COL_SIZE * sizeof(bool));
    }

    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;
    GapOpen = 0.0;
    GapExten = 0.0;
    FreeEndGaps = true;
    nGaps = 0;

    signalA_capacity = ROW_SIZE-1;
    signalB_capacity = COL_SIZE-1;
  }

  /// Reset object to initial state (without allocating new memory)
  void reset(int ROW_SIZE, int COL_SIZE)
  {
    if (ROW_SIZE -1 > signalA_capacity || COL_SIZE -1 > signalB_capacity)
    {
      std::cout << "Error: cannot reset an object beyond capacity" << std::endl;
      std::cout << ROW_SIZE << " vs " << signalA_capacity << std::endl;
      throw 1;
    }

    // resetting all values to zero.
    // We could use a for-loop but memset is faster for contiguous location in memory.
    std::memset(s_data, 0, (ROW_SIZE -1) * (COL_SIZE-1) * sizeof(double));
    std::memset(M, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(A, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(B, 0, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memset(Traceback, Traceback::SS, 3 * ROW_SIZE * COL_SIZE * sizeof(Traceback::TracebackType));
    std::memset(Path, 0, ROW_SIZE * COL_SIZE * sizeof(bool));
    std::memset(simPath, 0, ROW_SIZE * COL_SIZE * sizeof(bool));

    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;
    GapOpen = 0.0;
    GapExten = 0.0;
    FreeEndGaps = true;
    indexA_aligned.clear();
    indexB_aligned.clear();
    score.clear();
    nGaps = 0;
  }

  /**
   * @brief Overloading copy assignment operator.
   */
  AffineAlignObj& operator=(const AffineAlignObj& rhs)
  {
    freeMemory_();
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
    nGaps = rhs.nGaps;

    int ROW_SIZE = rhs.signalA_len + 1;
    int COL_SIZE = rhs.signalB_len + 1;

    allocateMemory_(ROW_SIZE, COL_SIZE);
    copyData_(rhs, ROW_SIZE, COL_SIZE);
    return *this;
  }

  /**
   * @brief Copy constructor.
   */
  AffineAlignObj(const AffineAlignObj& rhs)
  {
    freeMemory_();

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
    nGaps = rhs.nGaps;

    int ROW_SIZE = rhs.signalA_len + 1;
    int COL_SIZE = rhs.signalB_len + 1;

    allocateMemory_(ROW_SIZE, COL_SIZE);
    copyData_(rhs, ROW_SIZE, COL_SIZE);
  }

  /// Destructor: frees memory.
  ~AffineAlignObj()
  {
    freeMemory_();
  }

private:

  // Should be deleted?
  AffineAlignObj() {}

  /// Frees memory for s_data, M, A, B, Traceback, Path matrices.
  void freeMemory_()
  {
    delete[] s_data;
    delete[] M;
    delete[] A;
    delete[] B;
    delete[] Traceback;
    delete[] Path;
    delete[] simPath;
  }

  /// Allocates memory for s_data, M, A, B, Traceback, Path matrices.
  void allocateMemory_(int ROW_SIZE, int COL_SIZE)
  {
    // new allocate memory in heap. Here we just keep the memory address in our object's member.
    // Memory will remain valid outside of this constructor's scope.
    // Therefore, we need to explicitly free it.
    s_data = new double[(ROW_SIZE -1) * (COL_SIZE-1)];
    M = new double[ROW_SIZE * COL_SIZE];
    A = new double[ROW_SIZE * COL_SIZE];
    B = new double[ROW_SIZE * COL_SIZE];
    Traceback = new Traceback::TracebackType[3* ROW_SIZE * COL_SIZE];
    Path = new bool[ROW_SIZE * COL_SIZE];
    simPath = new bool[ROW_SIZE * COL_SIZE];
  }

  /// Memory copy s_data, M, A, B, Traceback, Path matrices.
  void copyData_(const AffineAlignObj& rhs, int ROW_SIZE, int COL_SIZE)
  {
    std::memcpy(s_data, rhs.s_data, (ROW_SIZE -1) * (COL_SIZE-1) * sizeof(double));
    std::memcpy(M, rhs.M, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(A, rhs.A, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(B, rhs.B, ROW_SIZE * COL_SIZE * sizeof(double));
    std::memcpy(Traceback, rhs.Traceback, 3 *ROW_SIZE * COL_SIZE * sizeof(Traceback::TracebackType));
    std::memcpy(Path, rhs.Path, ROW_SIZE * COL_SIZE * sizeof(bool));
    std::memcpy(simPath, rhs.simPath, ROW_SIZE * COL_SIZE * sizeof(bool));
  }
};
} // namespace DIAlign

#endif // AFFINEALIGNOBJ_H

// brackets ([]) specify the index of an element of the array. In fact these brackets are a dereferencing operator known as offset operator.
// We used bracked to work with STL containers. In case of arrays/ pointers, they work in the same way.
// double* M;
// M[i] == *(M+i);
