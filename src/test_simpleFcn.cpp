#include <iostream>
#include <vector>
#include <assert.h>
#include "simpleFcn.h"

#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

void test_getseqSim(){
  double Match=10.0, MisMatch=-2.0;
  std::string seq1 = "GCAT";
  std::string seq2 = "CAGTG";
  int seq1Len = seq1.size();
  int seq2Len = seq2.size();
  SimMatrix s = getseqSim(seq1, seq2, Match, MisMatch);

  // std::cout << "Similarity matrix is : " << std::endl;
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

  for (int i = 0; i < seq1Len; i++)
    for (int j = 0; j < seq2Len; j++)
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-07);

  // ASSERT(s[0][0] == -2);
  // ASSERT(s[1][0] == 10);
  // ASSERT(s[2][0] == -2);
  // ASSERT(s[3][0] == -2);

  // ASSERT(s[0][1] == -2);
  // ASSERT(s[1][1] == -2);
  // ASSERT(s[2][1] == 10);
  // ASSERT(s[3][1] == -2);

  // ASSERT(s[0][2] == 10);
  // ASSERT(s[1][2] == -2);
  // ASSERT(s[2][2] == -2);
  // ASSERT(s[3][2] == -2);

  // ASSERT(s[0][3] == -2);
  // ASSERT(s[1][3] == -2);
  // ASSERT(s[2][3] == -2);
  // ASSERT(s[3][3] == 10);

  // ASSERT(s[0][4] == 10);
  // ASSERT(s[1][4] == -2);
  // ASSERT(s[2][4] == -2);
  // ASSERT(s[3][4] == -2);

  Match=0.0, MisMatch=-0.0;
  s = getseqSim(seq1, seq2, Match, MisMatch);
  // std::cout << "Similarity matrix is : " << std::endl;
  // 0 0 0 0 0
  // 0 0 0 0 0
  // 0 0 0 0 0
  // 0 0 0 0 0

  std::vector< std::vector< double > > cmp_arr1;
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < seq1Len; i++)
    for (int j = 0; j < seq2Len; j++)
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-07);

  // ASSERT(s[0][0] == 0);
  // ASSERT(s[1][0] == 0);
  // ASSERT(s[2][0] == 0);
  // ASSERT(s[3][0] == 0);

  // ASSERT(s[0][1] == 0);
  // ASSERT(s[1][1] == 0);
  // ASSERT(s[2][1] == 0);
  // ASSERT(s[3][1] == 0);

  // ASSERT(s[0][2] == 0);
  // ASSERT(s[1][2] == 0);
  // ASSERT(s[2][2] == 0);
  // ASSERT(s[3][2] == 0);

  // ASSERT(s[0][3] == 0);
  // ASSERT(s[1][3] == 0);
  // ASSERT(s[2][3] == 0);
  // ASSERT(s[3][3] == 0);

  // ASSERT(s[0][4] == 0);
  // ASSERT(s[1][4] == 0);
  // ASSERT(s[2][4] == 0);
  // ASSERT(s[3][4] == 0);

  //TODO: How to test for such cases?
  seq1 = "";
  seq2 = "GCAT";
  s = getseqSim(seq1, seq2, Match, MisMatch);
  // std::cout << "Similarity matrix is : " << std::endl;
  // 0 0 0 0 0
}

#ifdef USE_Rcpp
int main_simpleFcn(){
#else
int main(){
#endif
  test_getseqSim();
  std::cout << "test simpleFcn successful" << std::endl;
  return 0;
}
