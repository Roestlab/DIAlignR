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

  //std::cout.precision(10);
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
  std::vector< double > vec;
  vec = {-2.0, 2.5, 0.9, -0.50, 0.0};
  double minVal = -1.0;
  double maxVal = 1.0;
  clamp(vec, minVal, maxVal);

  std::vector< double > cmp = {-1.0, 1.0, 0.9, -0.50, 0.0};
  for(int i = 0; i < 5; i++){
    ASSERT(std::abs(vec[i] - cmp[i]) < 1e-06);
  }
}

void test_meanNormalizeVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, 3, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, 1, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, 1, 0, -2}; vov.push_back(tmp);
  std::vector< std::vector< double > > vov_new = meanNormalizeVecOfVec(vov);
  // mean = 1.35

  std::vector< std::vector< double > > cmp_arr;
  tmp = {-1.481481, -1.481481, 7.407407, 0.0, 7.407407}; cmp_arr.push_back(tmp);
  tmp = {7.407407, -1.481481, 2.2222222, 0.0, -1.481481}; cmp_arr.push_back(tmp);
  tmp = {-1.481481, 7.407407, 0.7407407, 0.0, -1.481481}; cmp_arr.push_back(tmp);
  tmp = {-1.481481, -1.481481, 0.7407407, 0.0, -1.481481}; cmp_arr.push_back(tmp);

  int n_frag = vov_new.size(); // 4
  for (int i = 0; i < n_frag; i++){
    for (int j = 0; j < vov_new[i].size(); j++){
      ASSERT(std::abs(vov_new[i][j] - cmp_arr[i][j]) < 1e-06);
    }
  }
}

void test_L2NormalizeVecOfVec(){
  std::vector< std::vector< double > > vov;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, 3, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, 1, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, 1, 0, -2}; vov.push_back(tmp);
  std::vector< std::vector< double > > vov_new = L2NormalizeVecOfVec(vov);
  // mean = 21.14237

  std::vector< std::vector< double > > cmp_arr;
  tmp = {-0.09459677, -0.09459677, 0.47298387, 0.0, 0.47298387}; cmp_arr.push_back(tmp);
  tmp = {0.47298387, -0.09459677, 0.14189516, 0.0, -0.09459677}; cmp_arr.push_back(tmp);
  tmp = {-0.09459677, 0.47298387, 0.04729839, 0.0, -0.09459677}; cmp_arr.push_back(tmp);
  tmp = {-0.09459677, -0.09459677, 0.04729839, 0.0, -0.09459677}; cmp_arr.push_back(tmp);

  int n_frag = vov_new.size(); // 4
  for (int i = 0; i < n_frag; i++){
    for (int j = 0; j < vov_new[i].size(); j++){
      ASSERT(std::abs(vov_new[i][j] - cmp_arr[i][j]) < 1e-06);
    }
  }
}

void test_divideVecOfVec(){
  std::vector< std::vector< double > > vov;
  double num = 21.14237;
  std::vector< double > tmp;
  tmp = {-2.0, -2, 10, 0, 10}; vov.push_back(tmp);
  tmp = {10, -2, 3, 0, -2}; vov.push_back(tmp);
  tmp = {-2, 10, 1, 0, -2}; vov.push_back(tmp);
  tmp = {-2, -2, 1, 0, -2}; vov.push_back(tmp);
  std::vector< std::vector< double > > vov_new = divideVecOfVec(vov, num);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {-0.09459677, -0.09459677, 0.47298387, 0.0, 0.47298387}; cmp_arr.push_back(tmp);
  tmp = {0.47298387, -0.09459677, 0.14189516, 0.0, -0.09459677}; cmp_arr.push_back(tmp);
  tmp = {-0.09459677, 0.47298387, 0.04729839, 0.0, -0.09459677}; cmp_arr.push_back(tmp);
  tmp = {-0.09459677, -0.09459677, 0.04729839, 0.0, -0.09459677}; cmp_arr.push_back(tmp);

  int n_frag = vov_new.size(); // 4
  for (int i = 0; i < n_frag; i++){
    for (int j = 0; j < vov_new[i].size(); j++){
      ASSERT(std::abs(vov_new[i][j] - cmp_arr[i][j]) < 1e-06);
    }
  }
}

void test_ElemWiseSumOuterProd(){
  std::vector< double > d1 = {1.0, 0.0, 0.5, -0.2};
  std::vector< double > d2 = {11.0, -6.0, 0.0};
  SimMatrix s;
  s.n_row = 4;
  s.n_col = 3;

  //........................  CASE 1 ........................................
  s.data.resize(3*4, 0.0);
  ElemWiseSumOuterProd(d1, d2, s);

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {11.0, -6.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {5.5, -3.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {-2.2, 1.2, 0.0}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s.data.resize(3*4, 10); // Only new elements will be initialized with value = 10.0
  ElemWiseSumOuterProd(d1, d2, s);
  std::vector< std::vector< double > > cmp_arr1;
  tmp = {22.0, -12.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {11, -6.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {-4.4, 2.4, 0.0}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-06);
    }
  }
}

void test_ElemWiseSumOuterProdMeanSub(){
  std::vector< double > d1 = {1.0, 0.0, 0.5, -0.2};
  std::vector< double > d2 = {11.0, -6.0, 0.0};
  std::vector< double > m1 = {0.3, 0.0, 0.0, 0.3};
  std::vector< double > m2 = {5.0, -1.0, 0.2};
  SimMatrix s;
  s.n_row = 4;
  s.n_col = 3;

  //........................  CASE 1 ........................................
  s.data.resize(3*4, 0.0);
  ElemWiseSumOuterProdMeanSub(d1, d2, s, m1, m2);

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {4.2, -3.5, -0.14}; cmp_arr.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {3.0, -2.5, -0.1}; cmp_arr.push_back(tmp);
  tmp = {-3.0, 2.5, 0.1}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s.data.resize(3*4, 10); // Only new elements will be initialized with value = 10.0
  ElemWiseSumOuterProdMeanSub(d1, d2, s, m1, m2);
  std::vector< std::vector< double > > cmp_arr1;
  tmp = {8.4, -7.0, -0.28}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {6, -5.0, -0.2}; cmp_arr1.push_back(tmp);
  tmp = {-6.0, 5.0, 0.2}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-06);
    }
  }
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
