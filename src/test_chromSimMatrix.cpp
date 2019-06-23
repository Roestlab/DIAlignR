#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "chromSimMatrix.h"

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

void test_meanVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, -2, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, -2, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, -2, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 10, -2}; vov.push_back(tmp);

  double mean = meanVecOfVec(vov);
  ASSERT(std::abs(mean - 1.0) < 1e-6);
}

void test_eucLenVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, -2, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, -2, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, -2, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 10, -2}; vov.push_back(tmp);

  double len = eucLenVecOfVec(vov);
  ASSERT(std::abs(len - 23.66432) < 1e-6);
}

void test_perSampleEucLenVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 0, -2}; vov.push_back(tmp);

  std::vector< double > vec = perSampleEucLenVecOfVec(vov);
  tmp = {10.583005, 10.583005, 10.583005, 0.0, 10.583005};
  for (int i = 0; i < vec.size(); i++){
    ASSERT(std::abs(vec[i] - tmp[i]) < 1e-6);
  }
}

void test_perSampleSqrSumVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 9, -2}; vov.push_back(tmp);

  std::vector< double > vec = perSampleSqrSumVecOfVec(vov);
  tmp = {112.0, 112.0, 112.0, 81.0, 112.0};
  for (int i = 0; i < vec.size(); i++){
    ASSERT(std::abs(vec[i] - tmp[i]) < 1e-6);
  }
}

void test_perSampleMeanVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 9, -2}; vov.push_back(tmp);

  std::vector< double > vec = perSampleMeanVecOfVec(vov);
  tmp = {1.0, 1.0, 1.0, 2.25, 1.0};
  for (int i = 0; i < vec.size(); i++){
    ASSERT(std::abs(vec[i] - tmp[i]) < 1e-6);
  }
}

void test_perSampleSumVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, -2, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, -2, 9, -2}; vov.push_back(tmp);

  std::vector< double > vec = perSampleSumVecOfVec(vov);
  tmp = {4.0, 4.0, 4.0, 9.0, 4.0};
  for (int i = 0; i < vec.size(); i++){
    ASSERT(std::abs(vec[i] - tmp[i]) < 1e-6);
  }
}

void test_distToSim(){
  SimMatrix s;
  //........................  CASE1 ........................................
  s.data = {-2, 2, 10, -2, 10,
            10, 2.5, -2, -2, -2,
            -2, 10, -2, 0, -2,
            -2, -2, -1, 10, -2};
  s.n_col = 5;
  s.n_row = 4;
  double offset = 1, Numerator = 1;
  distToSim(s, offset, Numerator);

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {-1.000001, 0.33333322, 0.09090908, -1.000001, 0.09090908}; cmp_arr.push_back(tmp);
  tmp = {0.09090908, 0.2857142, -1.000001, -1.000001, -1.000001}; cmp_arr.push_back(tmp);
  tmp = {-1.000001, 0.09090908, -1.000001, 0.999999, -1.000001}; cmp_arr.push_back(tmp);
  tmp = {-1.000001, -1.000001, 1e06, 0.09090908, -1.000001}; cmp_arr.push_back(tmp);

  std::cout.precision(10);
  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-04);
      // Threshold is 1e-04 instead of 1e-06 because of subtraction of 1e06 - 1e06 at s(4,4)
    }
  }

  //........................  CASE2 ........................................
  s.data = {-2, 2, 10, -2, 10,
            10, 2.5, -2, -2, -2,
            -2, 10, -2, 0, -2,
            -2, -2, -1, 10, -2};
  offset = 1, Numerator = 0;
  distToSim(s, offset, Numerator);

  std::vector< std::vector< double > > cmp_arr1;
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-06);
    }
  }
}

void test_clamp(){

}

void test_meanNormalizeVecOfVec(){

}

void test_L2NormalizeVecOfVec(){

}

void test_divideVecOfVec(){

}

void test_ElemWiseSumOuterProd(){

}

void test_ElemWiseSumOuterProdMeanSub(){

}

void test_ElemWiseSumOuterEucl(){

}

void test_ElemWiseOuterCosine(){

}

void test_SumOuterProd(){

}

void test_SumOuterCov(){

}

void test_SumOuterCorr(){

}

void test_SumOuterEucl(){

}

void test_SumOuterCosine(){

}

void test_getSimilarityMatrix(){

}

#ifdef USE_Rcpp
int main_chromSimMatrix(){
#else
int main(){
#endif
  test_meanVecOfVec();
  test_eucLenVecOfVec();
  test_perSampleEucLenVecOfVec();
  test_perSampleSqrSumVecOfVec();
  test_perSampleMeanVecOfVec();
  test_perSampleSumVecOfVec();
  test_distToSim();
  test_clamp();
  test_meanNormalizeVecOfVec();
  test_L2NormalizeVecOfVec();
  test_divideVecOfVec();
  test_ElemWiseSumOuterProd();
  test_ElemWiseSumOuterProdMeanSub();
  test_ElemWiseSumOuterEucl();
  test_ElemWiseOuterCosine();
  test_SumOuterProd();
  test_SumOuterCov();
  test_SumOuterCorr();
  test_SumOuterEucl();
  test_SumOuterCosine();
  test_getSimilarityMatrix();
  std::cout << "test chromSimMatrix successful" << std::endl;
  return 0;
}
