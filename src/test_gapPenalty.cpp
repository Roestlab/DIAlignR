#include <cmath> // require for std::abs
#include <assert.h>
#include "gapPenalty.h"
#include "utils.h" //To propagate #define DIALIGN_USE_Rcpp

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;

void test_getGapPenalty(){
  SimMatrix s;
  s.data = {
    -2.19074364,  1.05454755,  5.26353599, -2.89112702,  0.25924473,  0.89726085,
    2.62386419, -0.21909407,  6.45342181,  0.08363896,  1.75295225,  3.44525833,
    -0.67808607, -2.61900693,  5.84668688, -6.43320725,  3.13581374,  0.60742015,
    3.53848608,  1.79679546,  6.77245762, -3.09977746,  5.26891460,  6.36395493,
    0.51481333, -6.85511916,  1.93171191, -1.28967451,  2.87660981,  1.36891013,
    2.71681581,  1.45688120, 3.72849306, -0.35247316, -1.83002582, -1.28698150,
    -4.67793934, -2.20775344, -1.17718574, -0.23953770, -0.65075869, -5.37730953,
    -2.02511518,  6.21064240,  2.36748179,  6.47276131, -0.41645117,  0.22746729,
    -0.05248435, -3.09630330, -2.01486144,  6.69890407, -1.18674116,  4.32714654,
    -2.64271788, -5.39763472, -0.46891328,  3.30758758,  3.91768941,  5.51485630,
    -4.86472662,  6.59372756, -1.60943300,  0.97449429
  }; //set.seed(2); rnorm(64, 0.5, 3)
  s.n_col = 8;
  s.n_row = 8;

  double gPM70 = getGapPenalty(s, 0.7, "dotProductMasked");
  double gP80 = getGapPenalty(s, 0.8, "dotProduct");
  double gCA75 = getGapPenalty(s, 0.75, "cosineAngle");
  double gC2A70 = getGapPenalty(s, 0.7, "cosine2Angle");
  double gED20 = getGapPenalty(s, 0.2, "euclideanDist");
  double gCOV90 = getGapPenalty(s, 0.90, "covariance");
  double gCOR95 = getGapPenalty(s, 0.95, "correlation");
  double gNONE50 = getGapPenalty(s, 0.5, "NONE");

  ASSERT(std::abs(gPM70 -  2.732795) < 1e-6);
  ASSERT(std::abs(gP80 -  3.804172) < 1e-6);
  ASSERT(std::abs(gCA75 -  0.95) < 1e-6);
  ASSERT(std::abs(gC2A70 - 0.95) < 1e-6);
  ASSERT(std::abs(gED20 -  -2.091367) < 1e-6);
  ASSERT(std::abs(gCOV90 - 6.101456) < 1e-6);
  ASSERT(std::abs(gCOR95 - 6.46986) < 1e-6);
  // TODO How to check this case?
  //ASSERT(std::abs(gNONE50 - 0.0) < 1e-6);

  std::fill(s.data.begin(), s.data.end(), 0.0);
  gPM70 = getGapPenalty(s, 0.7, "dotProductMasked");
  gP80 = getGapPenalty(s, 0.8, "dotProduct");
  gCA75 = getGapPenalty(s, 0.75, "cosineAngle");
  gC2A70 = getGapPenalty(s, 0.7, "cosine2Angle");
  gED20 = getGapPenalty(s, 0.2, "euclideanDist");
  gCOV90 = getGapPenalty(s, 0.90, "covariance");
  gCOR95 = getGapPenalty(s, 0.95, "correlation");
  gNONE50 = getGapPenalty(s, 0.5, "NONE");

  ASSERT(std::abs(gPM70 - 0.0) < 1e-6);
  ASSERT(std::abs(gP80 - 0.0) < 1e-6);
  ASSERT(std::abs(gCA75 - 0.95) < 1e-6);
  ASSERT(std::abs(gC2A70 - 0.95) < 1e-6);
  ASSERT(std::abs(gED20 - 0.0) < 1e-6);
  ASSERT(std::abs(gCOV90 - 0.0) < 1e-6);
  ASSERT(std::abs(gCOR95 - 0.0) < 1e-6);
  // TODO How to check this case?
  //ASSERT(std::abs(gNONE50 - 0.0) < 1e-6);

  // TODO make gapPenalty robust against high entries of zero.
  // Thinking of taking maximum/20 = 5% noise cutoff.
  s.data = {
    0.0,  0.0,  0.0, 0.0,  -0.00473,  0.00085,
    0.0, -0.009407,  0.0,  0.00363896,  0.0,  0.0,
    -0.00808607, -0.0000693,  0.00068688, 0.0,  0.0,  0.005,
    0.0,  0.0, 0.0, -0.09977746,  0.0,  6.36395493,
    0.51481333, 6.85511916,  1.93171191, 1.28967451,  2.87660981,  1.36891013,
    0.0,  0.0, 0.0, -0.0, -0.0, 0.0,
    0.0, 0.0, 0.0, -0.0003770, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0, -0.00045117,  0.00046729,
    0.00048435, 0.00630330, -0.0,  0.0, 0.0, 0.0,
    0.0, 0.0, -0.000328,  0.0, 0.0,  0.0,
    -0.00072662, 0.0, 0.0, 0.0
  };

  // std::cout.precision(10);

  gPM70 = getGapPenalty(s, 0.7, "dotProductMasked");
  gP80 = getGapPenalty(s, 0.8, "dotProduct");
  gCA75 = getGapPenalty(s, 0.75, "cosineAngle");
  gC2A70 = getGapPenalty(s, 0.7, "cosine2Angle");
  gED20 = getGapPenalty(s, 0.2, "euclideanDist");
  gCOV90 = getGapPenalty(s, 0.90, "covariance");
  gCOR95 = getGapPenalty(s, 0.95, "correlation");
  gNONE50 = getGapPenalty(s, 0.5, "NONE");

  ASSERT(std::abs(gPM70 - 0.0) < 1e-6);
  ASSERT(std::abs(gP80 - 0.000474114) < 1e-6);
  ASSERT(std::abs(gCA75 - 0.95) < 1e-6);
  ASSERT(std::abs(gC2A70 - 0.95) < 1e-6);
  ASSERT(std::abs(gED20 - 0.0) < 1e-6);
  ASSERT(std::abs(gCOV90 - 0.3622603) < 1e-6);
  ASSERT(std::abs(gCOR95 - 1.8472916) < 1e-6);
  // TODO How to check this case?
  //ASSERT(std::abs(gNONE50 - 0.0) < 1e-6);
}

#ifdef DIALIGN_USE_Rcpp
int main_gapPenalty(){
#else
int main(){
#endif
  test_getGapPenalty();
  std::cout << "test gapPenalty successful" << std::endl;
  return 0;
}
