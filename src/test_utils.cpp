#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "utils.h"

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;

void test_getQuantile(){
  std::vector<double> vec = {
    -2.19074364,  1.05454755,  5.26353599, -2.89112702,  0.25924473,  0.89726085,
    2.62386419, -0.21909407,  6.45342181,  0.08363896,  1.75295225,  3.44525833,
    -0.67808607, -2.61900693,  5.84668688, -6.43320725,  3.13581374,  0.60742015,
    3.53848608,  1.79679546,  6.77245762, -3.09977746,  5.26891460,  6.36395493,
    0.51481333, -6.85511916,  1.93171191, -1.28967451,  2.87660981,  1.36891013,
    2.71681581,  1.45688120, 3.72849306, -0.35247316, -1.83002582, -1.28698150,
    -4.67793934, -2.20775344, -1.17718574, -0.23953770, -0.65075869, -5.37730953,
    -2.02511518,  6.21064240,  2.36748179,  6.47276131, -0.41645117,  0.22746729,
    -0.05248435, -3.09630330, -2.01486144,  6.69890407, -1.18674116,  4.32714654,
    -2.64271788, -5.39763472, -0.46891328,  3.30758758,  3.91768941,  5.51485630
  }; //set.seed(2); rnorm(60, 0.5, 3)

  double q10 = getQuantile(vec, 0.1);
  double q50 = getQuantile(vec, 0.5);
  double q75 = getQuantile(vec, 0.75);
  double q95 = getQuantile(vec, 0.95);

  ASSERT( std::abs(q10 -  -3.096651) < 1e-6);
  ASSERT( std::abs(q50 -  0.387029) < 1e-6);
  ASSERT( std::abs(q75 -  3.342005) < 1e-6);
  ASSERT(std::abs(q95 - 6.454389) < 1e-6);

  std::fill(vec.begin(), vec.end(), 0.0);
  q10 = getQuantile(vec, 0.1);
  q50 = getQuantile(vec, 0.5);
  q75 = getQuantile(vec, 0.75);
  q95 = getQuantile(vec, 0.95);

  ASSERT( std::abs(q10 -  0.0) < 1e-6);
  ASSERT( std::abs(q50 -  0.0) < 1e-6);
  ASSERT( std::abs(q75 -  0.0) < 1e-6);
  ASSERT(std::abs(q95 - 0.0) < 1e-6);
}

#ifdef USE_Rcpp
int main_utils(){
#else
int main(){
#endif
  test_getQuantile();
  std::cout << "test utils successful" << std::endl;
  return 0;
}
