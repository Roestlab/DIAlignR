#include <cmath> // require for std::abs
#include <assert.h>
#include "constrainMat.h"
#include "utils.h" //To propagate #define USE_Rcpp

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;

void test_calcNoBeefMask(){
  double A1 = 3353.2, A2 = 3363.5;
  double B1 = 3325.9, B2 = 3339.5;
  double B1p = 3324.7, B2p = 3336.119;
  SimMatrix MASK;
  MASK.data.resize(5*4,0);
  MASK.n_col = 5;
  MASK.n_row = 4;

  //........................  CASE 1 ........................................
  int noBeef = 2;
  calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, false);
  // std::cout << "MASK is : " << std::endl;
  // 0 0 0 0.707 1.414
  // 0 0 0 0 0.707
  // 0 0 0 0 0
  // 0.707 0 0 0 0

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {0, 0, 0, 0.7071068, 1.4142136}; cmp_arr.push_back(tmp);
  tmp = {0, 0, 0, 0, 0.7071068}; cmp_arr.push_back(tmp);
  tmp = {0, 0, 0, 0, 0}; cmp_arr.push_back(tmp);
  tmp = {0.7071068, 0, 0, 0, 0}; cmp_arr.push_back(tmp);

  // std::cout.precision(10);
  for (int i = 0; i < MASK.n_row; i++)
    for (int j = 0; j < MASK.n_col; j++)
      ASSERT(std::abs(MASK.data[i*MASK.n_col+j] - cmp_arr[i][j]) < 1e-06);

  //........................  CASE 2 ........................................
  std::fill(MASK.data.begin(), MASK.data.end(), 0.0);
  noBeef = 1;
  calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, true);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0, 0, 1, 1, 1}; cmp_arr2.push_back(tmp);
  tmp = {0, 0, 0, 1, 1}; cmp_arr2.push_back(tmp);
  tmp = {1, 0, 0, 0, 1}; cmp_arr2.push_back(tmp);
  tmp = {1, 1, 0, 0, 0}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < MASK.n_row; i++)
    for (int j = 0; j < MASK.n_col; j++)
      ASSERT(std::abs(MASK.data[i*MASK.n_col+j] - cmp_arr2[i][j]) < 1e-07);

  //........................  CASE 3 ........................................
  std::fill(MASK.data.begin(), MASK.data.end(), 0.0);
  noBeef = 0;
  calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, true);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {0, 1, 1, 1, 1}; cmp_arr3.push_back(tmp);
  tmp = {1, 0, 1, 1, 1}; cmp_arr3.push_back(tmp);
  tmp = {1, 1, 0, 1, 1}; cmp_arr3.push_back(tmp);
  tmp = {1, 1, 1, 0, 1}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < MASK.n_row; i++)
    for (int j = 0; j < MASK.n_col; j++)
      ASSERT(std::abs(MASK.data[i*MASK.n_col+j] - cmp_arr3[i][j]) < 1e-07);

}

void test_constrainSimilarity(){
  SimMatrix s, MASK;
  //........................  CASE1 ........................................
  s.data = {-2, -2, 10, -2, 10,
            10, -2, -2, -2, -2,
            -2, 10, -2, -2, -2,
            -2, -2, -2, 10, -2};
  MASK.data = {0.000, 0.000, 0.707, 1.414, 2.121,
               0.000, 0.000, 0.000, 0.707, 1.414,
               0.707, 0.000, 0.000, 0.000, 0.000,
               1.414, 0.707, 0.000, 0.000, 0.000};
  s.n_col = 5;
  s.n_row = 4;
  MASK.n_col = 5;
  MASK.n_row = 4;
  double constrainVal = 0.0;

  constrainSimilarity(s, MASK, constrainVal);
  // std::cout << "Constrained similarity matrix is : " << std::endl;
  // -2 -2 10 -2 10
  // 10 -2 -2 -2 -2
  // -2 10 -2 -2 -2
  // -2 -2 -2 10 -2

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {-2, -2, 10, -2, 10}; cmp_arr.push_back(tmp);
  tmp = {10, -2, -2, -2, -2}; cmp_arr.push_back(tmp);
  tmp = {-2, 10, -2, -2, -2}; cmp_arr.push_back(tmp);
  tmp = {-2, -2, -2, 10, -2}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++)
    for (int j = 0; j < s.n_col; j++)
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-07);

  //........................  CASE2 ........................................
  s.data = {2, 2, 10, 2, 10,
            10, 2, 2, 2, 2,
            2, 10, 2, 2, 2,
            2, 2, 2, 10, 2};
  constrainVal = -10.0;

  constrainSimilarity(s, MASK, constrainVal);
  // std::cout << "Constrained similarity matrix is : " << std::endl;
  // 2 2 2.93 -12.14 -11.21
  // 10 2 2 -5.07 -12.14
  // -5.07 10 2 2 2
  // -12.14 -5.07 2 10 2

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {2, 2, 2.93, -12.14, -11.21}; cmp_arr2.push_back(tmp);
  tmp = {10, 2, 2, -5.07, -12.14}; cmp_arr2.push_back(tmp);
  tmp = {-5.07, 10, 2, 2, 2}; cmp_arr2.push_back(tmp);
  tmp = {-12.14, -5.07, 2, 10, 2}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s.n_row; i++)
    for (int j = 0; j < s.n_col; j++)
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr2[i][j]) < 1e-07);

  //........................  CASE3 ........................................
  s.data = {2, 2, 10, 2, 10,
            10, 2, 2, 2, 2,
            2, 10, 2, 2, 2,
            2, 2, 2, 10, 2};
  constrainVal = 10.0;
  std::fill(MASK.data.begin(), MASK.data.end(), 0.0);
  constrainSimilarity(s, MASK, constrainVal);
  // std::cout << "Constrained similarity matrix is : " << std::endl;
  // 2 2 10 2 10
  // 10 2 2 2 2
  // 2 10 2 2 2
  // 2 2 2 10 2
  std::vector< std::vector< double > > cmp_arr3;
  tmp = {2, 2, 10, 2, 10}; cmp_arr3.push_back(tmp);
  tmp = {10, 2, 2, 2, 2}; cmp_arr3.push_back(tmp);
  tmp = {2, 10, 2, 2, 2}; cmp_arr3.push_back(tmp);
  tmp = {2, 2, 2, 10, 2}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr3[i][j]) < 1e-07);
    }
  }
  // ASSERT(cmp_arr3[0][0] == 2);
  // ASSERT(cmp_arr3[1][0] == 10);
  // ASSERT(cmp_arr3[2][0] == 2);
  // ASSERT(cmp_arr3[3][0] == 2);

  // ASSERT(cmp_arr3[0][1] == 2);
  // ASSERT(cmp_arr3[1][1] == 2);
  // ASSERT(cmp_arr3[2][1] == 10);
  // ASSERT(cmp_arr3[3][1] == 2);

  // ASSERT(cmp_arr3[0][2] == 10);
  // ASSERT(cmp_arr3[1][2] == 2);
  // ASSERT(cmp_arr3[2][2] == 2);
  // ASSERT(cmp_arr3[3][2] == 2);

  // ASSERT(cmp_arr3[0][3] == 2);
  // ASSERT(cmp_arr3[1][3] == 2);
  // ASSERT(cmp_arr3[2][3] == 2);
  // ASSERT(cmp_arr3[3][3] == 10);

  // ASSERT(cmp_arr3[0][4] == 10);
  // ASSERT(cmp_arr3[1][4] == 2);
  // ASSERT(cmp_arr3[2][4] == 2);
  // ASSERT(cmp_arr3[3][4] == 2);
}

#ifdef USE_Rcpp
int main_constrainMat(){
#else
int main(){
#endif
  test_calcNoBeefMask();
  test_constrainSimilarity();
  std::cout << "test constrainMat successful" << std::endl;
  return 0;
}
