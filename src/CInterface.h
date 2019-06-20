#ifndef DIALIGN_CINTERFACE_H
#define DIALIGN_CINTERFACE_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/// enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};
struct  CAffineAlignObj
{
  /// double* M; // Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.
  /// double* A; // Insert in sequence A, residue in A is aligned to gap in B. A(i,j) is the best score given that Ai is aligned to a gap in B.
  /// double* B; // Insert in sequence B, residue in B is aligned to gap in A. B(i,j) is the best score given that Bj is aligned to a gap in A.
  /// TracebackType* Traceback;
  /// bool* Path; // Path matrix would represent alignment path through similarity matrix as binary-hot encoding.
  /// int n_matrix; // length of the matriices (all have the same length except ...
  /// int signalA_len; // Number of data-points in signal A
  /// int signalB_len; // Number of data-points in signal B
  /// double GapOpen; // Penalty for Gap opening
  /// double GapExten; // Penalty for Gap extension
  /// // For single gap: Penalty = GapOpen
  /// // For two consecutive gaps: Penalty = GapOpen + GapExten
  /// // For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten
  /// bool FreeEndGaps; // True for Overlap alignment
  int* indexA_aligned; // Aligned signalA indices after affine alignment
  int* indexB_aligned; // Aligned signalB indices after affine alignment
  int n_indexA_aligned; // Aligned signalA indices after affine alignment
  int n_indexB_aligned; // Aligned signalB indices after affine alignment
  int n_score;  // length of the score array
  double* score;  // Score along the aligned path
} ;

  struct CAffineAlignObj alignChromatogramsC(
                             int nchrom,
                             double* r1, int ndata1,
                             double* r2, int ndata2,
                             const char* alignType,
                             double* tA, int ntA,
                             double* tB, int ntB,
                             const char* normalization, const char* simType,
                             double B1p, double B2p, int noBeef,
                             double goFactor, double geFactor,
                             double cosAngleThresh, bool OverlapAlignment,
                             double dotProdThresh, double gapQuantile,
                             bool hardConstrain, double samples4gradient);

#ifdef __cplusplus
}
#endif

#endif
