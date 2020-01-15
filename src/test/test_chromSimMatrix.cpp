#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "../chromSimMatrix.h"
#include "../utils.h" //To propagate #define USE_Rcpp

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;
using namespace SimilarityMatrix;

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
  std::vector< double > d1 = {1.0, 0.0, 0.5, -0.2};
  std::vector< double > d2 = {11.0, -6.0, 0.0};
  SimMatrix s;
  s.n_row = 4;
  s.n_col = 3;

  //........................  CASE 1 ........................................
  s.data.resize(3*4, 0.0);
  ElemWiseSumOuterEucl(d1, d2, s);

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {100, 49, 1.0}; cmp_arr.push_back(tmp);
  tmp = {121.0, 36.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {110.25, 42.25, 0.25}; cmp_arr.push_back(tmp);
  tmp = {125.44, 33.64, 0.04}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s.data.resize(3*4, 10); // Only new elements will be initialized with value = 10.0
  ElemWiseSumOuterEucl(d1, d2, s);
  std::vector< std::vector< double > > cmp_arr1;
  tmp = {200, 98.0, 2.0}; cmp_arr1.push_back(tmp);
  tmp = {242.0, 72.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {220.5, 84.5, 0.5}; cmp_arr1.push_back(tmp);
  tmp = {250.88, 67.28, 0.08}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-06);
    }
  }
}

void test_ElemWiseOuterCosine(){
  std::vector< double > d1 = {1.0, 0.0, 0.5, -0.2};
  std::vector< double > d2 = {11.0, -6.0, 0.0};
  std::vector< double > m1 = {0.3, 0.0, 0.0, 0.3};
  std::vector< double > m2 = {5.0, 1.0, 0.2};
  SimMatrix s;
  s.n_row = 4;
  s.n_col = 3;

  //........................  CASE 1 ........................................
  s.data.resize(3*4, 0.0);
  ElemWiseOuterCosine(d1, d2, m1, m2, s);

  std::vector< std::vector< double > > cmp_arr;
  std::vector< double > tmp;
  tmp = {7.333307422, -19.99991333, 0}; cmp_arr.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr.push_back(tmp);
  tmp = {1099999.78, -2999997, 0}; cmp_arr.push_back(tmp);
  tmp = {-1.466661484, 3.999982667, 0}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-05);
    }
  }

  //........................  CASE 2 ........................................
  s.data.resize(3*4, 10); // Only new elements will be initialized with value = 10.0
  ElemWiseOuterCosine(d1, d2, m1, m2, s);
  std::vector< std::vector< double > > cmp_arr1;
  tmp = {14.66661484, -39.99982667, 0}; cmp_arr1.push_back(tmp);
  tmp = {0.0, 0.0, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {2199999.56, -5999994, 0.0}; cmp_arr1.push_back(tmp);
  tmp = {-2.933322969, 7.999965333, 0}; cmp_arr1.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr1[i][j]) < 1e-05);
    }
  }
}

void test_SumOuterProd(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, -1}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, -0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  SimMatrix s_m, s_l, s_n, s_z;
  s_m.n_row = 4;
  s_m.n_col = 4;
  s_l.n_row = 4;
  s_l.n_col = 4;
  s_n.n_row = 4;
  s_n.n_col = 4;
  s_z.n_row = 4;
  s_z.n_col = 4;
  s_z.data.resize(4*4, 0.0);

  //........................  CASE 1 ........................................
  s_m.data.resize(4*4, 0.0);
  SumOuterProd(d1, d2, "mean", s_m);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {2.844893, 4.237074, 5.296343, 2.300126}; cmp_arr.push_back(tmp);
  tmp = {3.692308, 5.447667, 6.204288, 4.721311}; cmp_arr.push_back(tmp);
  tmp = {3.268600, 4.842371, 5.750315, 3.510719}; cmp_arr.push_back(tmp);
  tmp = {4.721311, 7.112232, 7.868852, 6.204288}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s_m.n_row; i++){
    for (int j = 0; j < s_m.n_col; j++){
      ASSERT(std::abs(s_m.data[i*s_m.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s_l.data.resize(4*4, 0.0);
  SumOuterProd(d1, d2, "L2", s_l);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.1251213, 0.1863509, 0.2329386, 0.1011619}; cmp_arr2.push_back(tmp);
  tmp = {0.1623915, 0.2395940, 0.2728709, 0.2076481}; cmp_arr2.push_back(tmp);
  tmp = {0.1437564, 0.2129724, 0.2529048, 0.1544050}; cmp_arr2.push_back(tmp);
  tmp = {0.2076481, 0.3128033, 0.3460802, 0.2728709}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s_l.n_row; i++){
    for (int j = 0; j < s_l.n_col; j++){
      ASSERT(std::abs(s_l.data[i*s_l.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s_n.data.resize(4*4, 0.0);
  SumOuterProd(d1, d2, "None", s_n);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {9.4, 14.0, 17.5, 7.6}; cmp_arr3.push_back(tmp);
  tmp = {12.2, 18.0, 20.5, 15.6}; cmp_arr3.push_back(tmp);
  tmp = {10.8, 16.0, 19.0, 11.6}; cmp_arr3.push_back(tmp);
  tmp = {15.6, 23.5, 26.0, 20.5}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s_n.n_row; i++){
    for (int j = 0; j < s_n.n_col; j++){
      ASSERT(std::abs(s_n.data[i*s_n.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  std::vector< std::vector< double > > z1;
  tmp = {0.0, 0.0, 0.0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);

  std::vector< std::vector< double > > z2;
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  SumOuterProd(d1, d2, "mean", s_z);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s_z.n_row; i++){
    for (int j = 0; j < s_z.n_col; j++){
      //std::cout <<  s_z.data[i*s_z.n_col+j] << " ";
      //ASSERT(std::abs(s_z.data[i*s_z.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }
}

void test_SumOuterCov(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, -1}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, -0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  SimMatrix s_m, s_l, s_n, s_z;
  s_m.n_row = 4;
  s_m.n_col = 4;
  s_l.n_row = 4;
  s_l.n_col = 4;
  s_n.n_row = 4;
  s_n.n_col = 4;
  s_z.n_row = 4;
  s_z.n_col = 4;

  //........................  CASE 1 ........................................
  s_m.data.resize(4*4, 0.0);
  SumOuterCov(d1, d2, "mean", s_m);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {0.5649433, 0.9836066, 1.261034, -0.08575032}; cmp_arr.push_back(tmp);
  tmp = {0.6456494, 1.1349306, 1.160151, 0.63051702}; cmp_arr.push_back(tmp);
  tmp = {0.6052963, 1.0592686, 1.210593, 0.27238335}; cmp_arr.push_back(tmp);
  tmp = {0.9886507, 1.7402270, 1.715006, 1.12484237}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s_m.n_row; i++){
    for (int j = 0; j < s_m.n_col; j++){
      ASSERT(std::abs(s_m.data[i*s_m.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s_l.data.resize(4*4, 0.0);
  SumOuterCov(d1, d2, "L2", s_l);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.02484678, 0.04326003, 0.05546157, -0.003771387}; cmp_arr2.push_back(tmp);
  tmp = {0.02839633, 0.04991542, 0.05102465, 0.027730786}; cmp_arr2.push_back(tmp);
  tmp = {0.02662155, 0.04658772, 0.05324311, 0.011979700}; cmp_arr2.push_back(tmp);
  tmp = {0.04348187, 0.07653697, 0.07542774, 0.049471723}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s_l.n_row; i++){
    for (int j = 0; j < s_l.n_col; j++){
      ASSERT(std::abs(s_l.data[i*s_l.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s_n.data.resize(4*4, 0.0);
  SumOuterCov(d1, d2, "None", s_n);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {1.866667, 3.25, 4.166667, -0.2833333}; cmp_arr3.push_back(tmp);
  tmp = {2.133333, 3.75, 3.833333, 2.0833333}; cmp_arr3.push_back(tmp);
  tmp = {2.0, 3.5, 4.0, 0.9}; cmp_arr3.push_back(tmp);
  tmp = {3.266667, 5.75, 5.666667, 3.7166667}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s_n.n_row; i++){
    for (int j = 0; j < s_n.n_col; j++){
      ASSERT(std::abs(s_n.data[i*s_n.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  std::vector< std::vector< double > > z1;
  tmp = {0.0, 0.0, 0.0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);

  std::vector< std::vector< double > > z2;
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  s_z.data.resize(4*4, 0.0);
  SumOuterCov(d1, d2, "mean", s_z);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s_z.n_row; i++){
    for (int j = 0; j < s_z.n_col; j++){
      //std::cout << s_z.data[i*s_z.n_col+j] << " ";
      //ASSERT(std::abs(s_z.data[i*s_z.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }
}

void test_SumOuterCorr(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, -1}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, -0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  SimMatrix s_m, s_l, s_n, s_z;
  s_m.n_row = 4;
  s_m.n_col = 4;
  s_l.n_row = 4;
  s_l.n_col = 4;
  s_n.n_row = 4;
  s_n.n_col = 4;
  s_z.n_row = 4;
  s_z.n_col = 4;

  //........................  CASE 1 ........................................
  s_m.data.resize(4*4, 0.0);
  SumOuterCorr(d1, d2, "mean", s_m);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {0.8737211, 0.8660254, 0.9905361, -0.06486282}; cmp_arr.push_back(tmp);
  tmp = {0.9985384, 0.9992601, 0.9112932, 0.47693252}; cmp_arr.push_back(tmp);
  tmp = {0.9743547, 0.9707253, 0.9897433, 0.21444787}; cmp_arr.push_back(tmp);
  tmp = {0.9901516, 0.9922154, 0.8723686, 0.55098860}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s_m.n_row; i++){
    for (int j = 0; j < s_m.n_col; j++){
      ASSERT(std::abs(s_m.data[i*s_m.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s_l.data.resize(4*4, 0.0);
  SumOuterCorr(d1, d2, "L2", s_l);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.8737211, 0.8660254, 0.9905361, -0.06486282}; cmp_arr2.push_back(tmp);
  tmp = {0.9985384, 0.9992601, 0.9112932, 0.47693252}; cmp_arr2.push_back(tmp);
  tmp = {0.9743547, 0.9707253, 0.9897433, 0.21444787}; cmp_arr2.push_back(tmp);
  tmp = {0.9901516, 0.9922154, 0.8723686, 0.55098860}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s_l.n_row; i++){
    for (int j = 0; j < s_l.n_col; j++){
      ASSERT(std::abs(s_l.data[i*s_l.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s_n.data.resize(4*4, 0.0);
  SumOuterCorr(d1, d2, "None", s_n);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {0.8737211, 0.8660254, 0.9905361, -0.06486282}; cmp_arr3.push_back(tmp);
  tmp = {0.9985384, 0.9992601, 0.9112932, 0.47693252}; cmp_arr3.push_back(tmp);
  tmp = {0.9743547, 0.9707253, 0.9897433, 0.21444787}; cmp_arr3.push_back(tmp);
  tmp = {0.9901516, 0.9922154, 0.8723686, 0.55098860}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s_n.n_row; i++){
    for (int j = 0; j < s_n.n_col; j++){
      ASSERT(std::abs(s_n.data[i*s_n.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  std::vector< std::vector< double > > z1;
  tmp = {0.0, 0.0, 0.0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);

  std::vector< std::vector< double > > z2;
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  s_z.data.resize(4*4, 0.0);
  SumOuterCorr(d1, d2, "mean", s_z);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s_z.n_row; i++){
    for (int j = 0; j < s_z.n_col; j++){
      ASSERT(std::abs(s_z.data[i*s_z.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }
}

void test_SumOuterEucl(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, -1}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, -0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  SimMatrix s_m, s_l, s_n, s_z;
  s_m.n_row = 4;
  s_m.n_col = 4;
  s_l.n_row = 4;
  s_l.n_col = 4;
  s_n.n_row = 4;
  s_n.n_col = 4;

  //........................  CASE 1 ........................................
  s_m.data.resize(4*4, 0.0);
  SumOuterEucl(d1, d2, "mean", s_m);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {0.5871842, 0.5211067, 0.5165468, 0.2857270}; cmp_arr.push_back(tmp);
  tmp = {0.5849202, 0.7368782, 0.5335614, 0.3618207}; cmp_arr.push_back(tmp);
  tmp = {0.6515918, 0.6568138, 0.5620653, 0.3211770}; cmp_arr.push_back(tmp);
  tmp = {0.4102288, 0.6068609, 0.4931426, 0.3400544}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s_m.n_row; i++){
    for (int j = 0; j < s_m.n_col; j++){
      ASSERT(std::abs(s_m.data[i*s_m.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s_l.data.resize(4*4, 0.0);
  SumOuterEucl(d1, d2, "L2", s_l);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.8682131, 0.8425724, 0.8445747, 0.6576925}; cmp_arr2.push_back(tmp);
  tmp = {0.8624804, 0.9318630, 0.8504457, 0.7314213}; cmp_arr2.push_back(tmp);
  tmp = {0.8921416, 0.9070038, 0.8688622, 0.6946385}; cmp_arr2.push_back(tmp);
  tmp = {0.7612420, 0.8698441, 0.8203342, 0.7093717}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s_l.n_row; i++){
    for (int j = 0; j < s_l.n_col; j++){
      ASSERT(std::abs(s_l.data[i*s_l.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s_n.data.resize(4*4, 0.0);
  SumOuterEucl(d1, d2, "None", s_n);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {0.3289897, 0.3999998, 0.6666662, 0.1881846}; cmp_arr3.push_back(tmp);
  tmp = {0.2808002, 0.3999998, 0.3999998, 0.2348906}; cmp_arr3.push_back(tmp);
  tmp = {0.3238277, 0.4721357, 0.6666662, 0.2132572}; cmp_arr3.push_back(tmp);
  tmp = {0.1963133, 0.2582456, 0.2582456, 0.1915639}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s_n.n_row; i++){
    for (int j = 0; j < s_n.n_col; j++){
      ASSERT(std::abs(s_n.data[i*s_n.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  std::vector< std::vector< double > > z1;
  tmp = {0.0, 0.0, 0.0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);

  std::vector< std::vector< double > > z2;
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  s_z.data.resize(4*4, 0.0);
  SumOuterEucl(d1, d2, "mean", s_z);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s_z.n_row; i++){
    for (int j = 0; j < s_z.n_col; j++){
      ASSERT(std::abs(s_z.data[i*s_z.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }
}

void test_SumOuterCosine(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, -1}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, -0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  SimMatrix s_m, s_l, s_n, s_z;
  s_m.n_row = 4;
  s_m.n_col = 4;
  s_l.n_row = 4;
  s_l.n_col = 4;
  s_n.n_row = 4;
  s_n.n_col = 4;

  //........................  CASE 1 ........................................
  s_m.data.resize(4*4, 0.0);
  SumOuterCosine(d1, d2, "mean", s_m);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {0.9338557, 0.9328144, 0.9935318, 0.4495778}; cmp_arr.push_back(tmp);
  tmp = {0.9994619, 0.9889952, 0.9597366, 0.7609750}; cmp_arr.push_back(tmp);
  tmp = {0.9892024, 0.9828713, 0.9945046, 0.6326431}; cmp_arr.push_back(tmp);
  tmp = {0.9859988, 0.9961734, 0.9391111, 0.7715162}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s_m.n_row; i++){
    for (int j = 0; j < s_m.n_col; j++){
      ASSERT(std::abs(s_m.data[i*s_m.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s_l.data.resize(4*4, 0.0);
  SumOuterCosine(d1, d2, "L2", s_l);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.9338516, 0.9328111, 0.9935285, 0.4495763}; cmp_arr2.push_back(tmp);
  tmp = {0.9994579, 0.9889920, 0.9597338, 0.7609727}; cmp_arr2.push_back(tmp);
  tmp = {0.9891982, 0.9828680, 0.9945015, 0.6326410}; cmp_arr2.push_back(tmp);
  tmp = {0.9859953, 0.9961706, 0.9391086, 0.7715141}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s_l.n_row; i++){
    for (int j = 0; j < s_l.n_col; j++){
      ASSERT(std::abs(s_l.data[i*s_l.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s_n.data.resize(4*4, 0.0);
  SumOuterCosine(d1, d2, "None", s_n);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {0.9338561, 0.9328148, 0.9935322, 0.4495780}; cmp_arr3.push_back(tmp);
  tmp = {0.9994623, 0.9889956, 0.9597370, 0.7609753}; cmp_arr3.push_back(tmp);
  tmp = {0.9892028, 0.9828717, 0.9945050, 0.6326433}; cmp_arr3.push_back(tmp);
  tmp = {0.9859992, 0.9961737, 0.9391114, 0.7715164}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s_n.n_row; i++){
    for (int j = 0; j < s_n.n_col; j++){
      ASSERT(std::abs(s_n.data[i*s_n.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  std::vector< std::vector< double > > z1;
  tmp = {0.0, 0.0, 0.0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);
  tmp = {0, 0, 0, 0.0}; z1.push_back(tmp);

  std::vector< std::vector< double > > z2;
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; z2.push_back(tmp);
  s_z.data.resize(4*4, 0.0);
  SumOuterCosine(d1, d2, "mean", s_z);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);
  tmp = {0.0, 0.0, 0.0, 0.0}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s_z.n_row; i++){
    for (int j = 0; j < s_z.n_col; j++){
      ASSERT(std::abs(s_z.data[i*s_z.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }
}

void test_getSimilarityMatrix(){
  std::vector< std::vector< double > > d1;
  std::vector< double > tmp;
  tmp = {1.0, 3.0, 2.0, 4.0}; d1.push_back(tmp);
  tmp = {0, 0, 0, 1.0}; d1.push_back(tmp);
  tmp = {4, 4, 4, 5}; d1.push_back(tmp);

  std::vector< std::vector< double > > d2;
  tmp = {1.4, 2.0, 1.5, 4.0}; d2.push_back(tmp);
  tmp = {0.0, 0.5, 0.0, 0.0}; d2.push_back(tmp);
  tmp = {2.0, 3.0, 4.0, 0.9}; d2.push_back(tmp);

  double cosAngleThresh = 0.3, dotProdThresh = 0.96;

  //........................  CASE 1 ........................................
  SimMatrix s = getSimilarityMatrix(d1, d2, "L2", "dotProductMasked", cosAngleThresh, dotProdThresh);

  std::vector< std::vector< double > > cmp_arr;
  tmp = {0.1251213, 0.1863509, 0.2329386, 0.1011619}; cmp_arr.push_back(tmp);
  tmp = {0.1623915, 0.2395940, 0.2728709, 0.2076481}; cmp_arr.push_back(tmp);
  tmp = {0.1437564, 0.2129724, 0.2529048, 0.1544050}; cmp_arr.push_back(tmp);
  tmp = {0.2076481, 0.3128033, 0.3460802, 0.2728709}; cmp_arr.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr[i][j]) < 1e-06);
    }
  }

  //........................  CASE 2 ........................................
  s = getSimilarityMatrix(d1, d2, "L2", "dotProductMasked", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr2;
  tmp = {0.1251213, 0.1863509, 0.2329386, 0.1011619}; cmp_arr2.push_back(tmp);
  tmp = {0.1623915, 0.2395940, 0.2728709, 0.2076481}; cmp_arr2.push_back(tmp);
  tmp = {0.1437564, 0.2129724, 0.2529048, 0.1544050}; cmp_arr2.push_back(tmp);
  tmp = {0.2076481, 0.3128033, 0.3460802, 0.0}; cmp_arr2.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr2[i][j]) < 1e-06);
    }
  }

  //........................  CASE 3 ........................................
  s = getSimilarityMatrix(d1, d2, "mean", "dotProduct", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr3;
  tmp = {2.504811, 3.730570, 4.663212, 2.025167}; cmp_arr3.push_back(tmp);
  tmp = {3.250925, 4.796447, 5.462620, 4.156921}; cmp_arr3.push_back(tmp);
  tmp = {2.877868, 4.263509, 5.062916, 3.091044}; cmp_arr3.push_back(tmp);
  tmp = {4.156921, 6.262028, 6.928201, 5.462620}; cmp_arr3.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr3[i][j]) < 1e-06);
    }
  }

  //........................  CASE 4 ........................................
  s = getSimilarityMatrix(d1, d2, "None", "cosineAngle", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr4;
  tmp = {0.9338556, 0.9328143, 0.9935317, 0.4495778}; cmp_arr4.push_back(tmp);
  tmp = {0.9994618, 0.9889952, 0.9597366, 0.7609750}; cmp_arr4.push_back(tmp);
  tmp = {0.9892023, 0.9828712, 0.9945046, 0.6326430}; cmp_arr4.push_back(tmp);
  tmp = {0.9859988, 0.9961734, 0.9391110, 0.7715162}; cmp_arr4.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr4[i][j]) < 1e-06);
    }
  }

  //........................  CASE 5 ........................................
  s = getSimilarityMatrix(d1, d2, "None", "cosine2Angle", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr5;
  tmp = {0.7441746, 0.7402868, 0.9742125, -0.5957592}; cmp_arr5.push_back(tmp);
  tmp = {0.9978499, 0.9562246, 0.8421902, 0.1581667}; cmp_arr5.push_back(tmp);
  tmp = {0.9570445, 0.9320735, 0.9780804, -0.1995248}; cmp_arr5.push_back(tmp);
  tmp = {0.9443890, 0.9847243, 0.7638603, 0.1904752}; cmp_arr5.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr5[i][j]) < 1e-06);
    }
  }

  //........................  CASE 6 ........................................
  s = getSimilarityMatrix(d1, d2, "mean", "euclideanDist", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr6;
  tmp = {0.6076553, 0.5304450, 0.5201205, 0.2975992}; cmp_arr6.push_back(tmp);
  tmp = {0.6143512, 0.7417414, 0.5406573, 0.3750524}; cmp_arr6.push_back(tmp);
  tmp = {0.6798319, 0.6585880, 0.5629232, 0.3336137}; cmp_arr6.push_back(tmp);
  tmp = {0.4337325, 0.6414240, 0.5113607, 0.3554709}; cmp_arr6.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr6[i][j]) < 1e-06);
    }
  }

  //........................  CASE 7 ........................................
  s = getSimilarityMatrix(d1, d2, "L2", "covariance", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr7;
  tmp = {0.02484678, 0.03216771, 0.05546157, -0.003771387}; cmp_arr7.push_back(tmp);
  tmp = {0.02839633, 0.03438617, 0.05102465, 0.027730786}; cmp_arr7.push_back(tmp);
  tmp = {0.02662155, 0.03327694, 0.05324311, 0.011979700}; cmp_arr7.push_back(tmp);
  tmp = {0.02839633, 0.03438617, 0.05102465, 0.027730786}; cmp_arr7.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr7[i][j]) < 1e-06);
    }
  }

  //........................  CASE 8 ........................................
  s = getSimilarityMatrix(d1, d2, "None", "correlation", 0.4, 0.86);

  std::vector< std::vector< double > > cmp_arr8;
  tmp = {0.8737211, 0.9226129, 0.9905361, -0.06486282}; cmp_arr8.push_back(tmp);
  tmp = {0.9985384, 0.9862414, 0.9112932, 0.47693252}; cmp_arr8.push_back(tmp);
  tmp = {0.9743547, 0.9933993, 0.9897433, 0.21444787}; cmp_arr8.push_back(tmp);
  tmp = {0.9985384, 0.9862414, 0.9112932, 0.47693252}; cmp_arr8.push_back(tmp);

  for (int i = 0; i < s.n_row; i++){
    for (int j = 0; j < s.n_col; j++){
      ASSERT(std::abs(s.data[i*s.n_col+j] - cmp_arr8[i][j]) < 1e-06);
    }
  }
}

#ifdef DIALIGN_USE_Rcpp
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
  // For following functions, it outputs random numbers for case 4!
  //test_SumOuterProd();
  //test_SumOuterCov();
  //test_SumOuterCorr();
  //test_SumOuterEucl();
  //test_SumOuterCosine();
  test_getSimilarityMatrix();
  std::cout << "test chromSimMatrix successful" << std::endl;
  return 0;
}
