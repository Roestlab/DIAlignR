#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "../miscell.h"
#include "../utils.h" //To propagate #define USE_Rcpp
#include "../spline.h"
//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;
using namespace tk;

// Anonymous namespace: Only valid for this file.
namespace {
std::vector<double> getT1(){
  std::vector<double> t1 = {4978.4,4981.8,4985.2,4988.6,4992,4995.4,4998.9,5002.3,5005.7,5009.1,
                            5012.5,5015.9,5019.3,5022.8,5026.2,5029.6,5033,5036.4,5039.8,5043.2,
                            5046.7,5050.1,5053.5,5056.9,5060.3,5063.7,5067.1,5070.5,5074,5077.4,
                            5080.8,5084.2,5087.6,5091,5094.4,5097.9,5101.3,5104.7,5108.1,5111.5,
                            5114.9,5118.3,5121.8,5125.2,5128.6,5132,5135.4,5138.8,5142.2,5145.7,
                            5149.1,5150.8,5152.5,5154.2,5155.9,5157.6,5159.3,5161,5162.7,5164.4,
                            5166.1,5167.85,5169.6,5171.3,5173,5174.7,5176.4,5178.1,5179.8,5181.5,
                            5183.2,5184.9,5186.6,5188.3,5190,5191.7,5193.4,5195.15,5196.9,5198.6,
                            5200.3,5202,5203.7,5205.4,5207.1,5208.8,5210.5,5212.2,5213.9,5215.6,
                            5217.3,5219.05,5220.8,5222.5,5224.2,5225.9,5227.6,5231,5234.4,5237.8,
                            5241.2,5244.7,5248.1,5251.5,5254.9,5258.3,5261.7,5265.1,5268.6,5272,
                            5275.4,5278.8,5282.2,5285.6,5289,5292.4,5295.9,5299.3,5302.7,5306.1,
                            5309.5,5312.9,5316.3,5319.8,5323.2,5326.6,5330,5333.4,5336.8,5340.2,
                            5343.7,5347.1,5350.5,5353.9,5357.3,5359,5360.7,5364.1,5367.6,5371,5374.4,
                            5377.8,5381.2,5384.6,5388,5391.4,5394.9,5398.3,5401.7,5405.1,5408.5,5411.9,
                            5415.3,5418.8,5422.2,5425.6,5429,5432.4,5435.8,5439.2,5442.7,5446.1,5449.5,
                            5452.9,5456.3,5459.7,5463.1,5466.5,5470,5473.4,5476.8,5480.2,5483.6,5487,
                            5490.4,5493.9,5497.3,5500.7,5504.1,5507.5,5510.9,5514.3,5517.8,5521.2,5524.6,
                            5528,5531.4,5534.8,5538.2,5539.95,5541.7,5545.1,5548.5,5550.2,5551.9,5555.3,
                            5558.7,5562.1,5563.8,5565.5,5569,5572.4,5574.1,5575.8,-1};
  return t1;
}

std::vector<double> getT2(){
  std::vector<double> t2 = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4988.6,4990.35,4992.1,4995.5,4998.9,5002.3,
                            5005.7,5009.1,5012.5,5015.9,5017.65,5019.4,5021.1,5022.8,5024.5,5026.2,5029.6,5031.3,
                            5033,5034.7,5036.4,5038.1,5039.8,5041.55,5043.3,5046.7,5050.1,5053.5,5056.9,5060.3,5063.7,5067.2,
                            5070.6,5074,5077.4,5080.8,5084.2,5087.6,5091,5094.5,5097.9,5101.3,5104.7,5108.1,5111.5,5114.9,
                            5118.4,5121.8,5125.2,5128.6,5132,5135.4,5138.8,5142.3,5145.7,5149.1,5152.5,5155.9,5159.3,5162.7,
                            5166.2,5169.6,5173,5176.4,5179.8,5183.2,5186.6,5190,5193.5,5196.9,5200.3,5203.7,5207.1,5210.5,
                            5213.9,5217.4,5220.8,5224.2,5227.6,5231,5234.4,5237.8,5241.3,5244.7,5248.1,5251.5,5254.9,5258.3,
                            5261.7,5265.2,5268.6,5272,5275.4,5278.8,5282.2,5285.6,5289.1,5290.8,5292.5,5294.2,5295.9,5299.3,
                            5302.7,5306.1,5309.5,5312.9,5316.4,5319.8,5323.2,5326.6,5330,5333.4,5336.8,5340.3,5343.7,5347.1,
                            5350.5,5353.9,5357.3,5360.7,5364.2,5367.6,5371,5374.4,5377.8,5381.2,5384.6,5388.1,5391.5,5394.9,
                            5398.3,5401.7,5405.1,5408.5,5412,5415.4,5418.8,5422.2,5425.6,5429,5432.4,5435.8,5439.3,5442.7,
                            5446.1,5449.5,5452.9,5456.3,5459.7,5463.2,5466.6,5470,5473.4,5476.8,5480.2,5483.6,5485.35,5487.1,
                            5490.5,5493.9,5497.3,5500.7,5504.1,5507.5,5511,5514.4,5517.8,5521.2,5524.6,5528,5531.4,5534.9,5538.3,
                            5541.7,5545.1,5548.5,5551.9,5555.3,5558.7,5562.2,5565.6,5569,5572.4,5575.8,5579.2,5582.6,5586.1};
  return t2;
}
}

void test_xicIntersect(){
    std::vector< std::vector< double > > time;
    std::vector< std::vector< double > > intensity;
    std::vector< std::vector< double > > timeExp;
    std::vector< std::vector< double > > intensityExp;
    std::vector<double> tmp;

    // Case 1
    tmp = {0, 3.4, 6.8, 9.4, 12.5}; time.push_back(tmp);
    tmp = {0, 3.4, 6.8, 9.4, 12.5}; time.push_back(tmp);
    tmp = {0, 3.4, 6.8, 9.4, 12.5}; time.push_back(tmp);
    tmp = {10, 14, 28, 34, 125}; intensity.push_back(tmp);
    tmp = {20, 34, 68, 94, 22.5}; intensity.push_back(tmp);
    tmp = {90, 31.4, 61.8, 92.4, 12.5}; intensity.push_back(tmp);

    tmp = {0, 3.4, 6.8, 9.4, 12.5}; timeExp.push_back(tmp);
    tmp = {0, 3.4, 6.8, 9.4, 12.5}; timeExp.push_back(tmp);
    tmp = {0, 3.4, 6.8, 9.4, 12.5}; timeExp.push_back(tmp);
    tmp = {10, 14, 28, 34, 125}; intensityExp.push_back(tmp);
    tmp = {20, 34, 68, 94, 22.5}; intensityExp.push_back(tmp);
    tmp = {90, 31.4, 61.8, 92.4, 12.5}; intensityExp.push_back(tmp);

    xicIntersect(time, intensity);
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 5; j++){
        ASSERT(std::abs(time[i][j] - timeExp[i][j]) < 1e-6);
        ASSERT(std::abs(intensity[i][j] - intensityExp[i][j]) < 1e-6);
      }
    }

    // Case 2
    time.clear();
    intensity.clear();
    timeExp.clear();
    intensityExp.clear();
    tmp = {1.0, 3.4, 6.8, 9.4, 12.5}; time.push_back(tmp);
    tmp = {1.0, 2.4, 4.8, 7.4, 12.1}; time.push_back(tmp);
    tmp = {1.5, 3.4, 6.8, 9.8, 12.4, 14.5}; time.push_back(tmp);
    tmp = {10, 14, 28, 34, 125}; intensity.push_back(tmp);
    tmp = {20, 34, 68, 94, 22.5}; intensity.push_back(tmp);
    tmp = {90, 31.4, 61.8, 89.6, 92.4, 12.5}; intensity.push_back(tmp);

    tmp = {1.0, 3.4, 6.8, 9.4, 12.5}; timeExp.push_back(tmp);
    tmp = {1.0, 2.4, 4.8, 7.4, 12.1}; timeExp.push_back(tmp);
    tmp = {1.5, 3.4, 6.8, 9.8, 12.4}; timeExp.push_back(tmp);
    tmp = {10, 14, 28, 34, 125}; intensityExp.push_back(tmp);
    tmp = {20, 34, 68, 94, 22.5}; intensityExp.push_back(tmp);
    tmp = {90, 31.4, 61.8, 89.6, 92.4}; intensityExp.push_back(tmp);

    xicIntersect(time, intensity);
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 5; j++){
        ASSERT(std::abs(time[i][j] - timeExp[i][j]) < 1e-6);
        ASSERT(std::abs(intensity[i][j] - intensityExp[i][j]) < 1e-6);
      }
    }
}

void test_interpolateZero(){
  std::vector<double> x = {-1, -1,2,3,-1, -1, 5, 9, -1, 10, -1 , -1};
  interpolateZero(x);
  std::vector<double> y = {-1,-1,2,3,3.667, 4.333, 5, 9, 9.5, 10, -1 , -1};
  for (int i = 0; i < x.size(); i++){
    ASSERT(std::abs(x[i] - y[i]) < 1e-2);
  }

}

void test_getKeep(){
  std::vector<int> x = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,51,53,55,57,59,61,63,65,67,69,
                           71,73,75,77,79,81,83,85,87,89,91,93,95,135,189,193,198,202,204};
  std::vector<int> keep = getKeep(205, x);
  std::vector<int> y = {18,19,20,21,22,23,24,25,26,27};
  for (int i = 0; i < 11; i++){
    ASSERT(keep[i] == y[i]);
  }
  ASSERT(keep.back() == 203);
}

void test_getFlank(){
  std::vector<double> t1 = getT1();
  std::vector<double> t2 = getT2();
  std::vector<int> flank = getFlank(t1, t2);
  std::vector<int> cmp_arr = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,204};
  for (int i = 0; i < cmp_arr.size(); i++){
    ASSERT(flank[i] - cmp_arr[i]);
  }
}

void test_getFlankN(){
  std::vector<int> flank = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,204};
  std::vector<double> t2 = getT2();
  std::vector<int> flank2 = getFlankN(t2, flank);
  std::vector<int> cmp_arr = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
  for (int i = 0; i < cmp_arr.size(); i++){
    ASSERT(flank2[i] == cmp_arr[i]);
  }

  std::vector<double> t1 = getT1();
  std::vector<int> flank1 = getFlankN(t1, flank);
  cmp_arr = {204};
  for (int i = 0; i < cmp_arr.size(); i++){
    ASSERT(flank1[i] == cmp_arr[i]);
  }
}

void test_addFlankToLeft(){
  std::vector<double> t = {3003.4, 3006.8, 3010.2, 3013.6, 3017.0, 3020.4, 3023.8,
                           3027.2, 3030.6, 3034.0, 3037.4, 3040.8, 3044.2, 3047.6};
  std::vector<std::vector<double>> inten{
    { 0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
      4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923}};
  inten.push_back(t);
  std::vector<double> tA(14, -1);
  std::vector<std::vector<double>> intenN{ { 1.2, 3.4, 5.6}};
  std::vector<double> tN = {3013.4, 3016.0, 3020.0};
  std::vector<int> flank{0,1,12,13};

  addFlankToLeft(t, tN, tA, inten, intenN, flank);

  std::vector<int> cmp_fk = {12, 13};
  std::vector<double> cmp_tN = {3006.6, 3010, 3013.4, 3016, 3020};
  std::vector<double> cmp_tA = {3006.6, 3010, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  std::vector< std::vector< double > > cmp_arr;
  cmp_arr.push_back({0.20506, 0.885007, 1.2, 3.4, 5.6});
  //........................  CASE 1 ........................................
  for (int i = 0; i < intenN.size(); i++)
    for (int j = 0; j < intenN[i].size(); j++)
      ASSERT(std::abs(intenN[i][j] - cmp_arr[i][j]) < 1e-06);

  for (int j = 0; j < tN.size(); j++)
    ASSERT(std::abs(tN[j] - cmp_tN[j]) < 1e-06);

  for (int j = 0; j < tA.size(); j++)
    ASSERT(std::abs(tA[j] - cmp_tA[j]) < 1e-06);

  for (int j = 0; j < flank.size(); j++)
    ASSERT( flank[j] == cmp_fk[j]);
}

void test_addFlankToRight(){
  std::vector<double> t = {3003.4, 3006.8, 3010.2, 3013.6, 3017.0, 3020.4, 3023.8,
                           3027.2, 3030.6, 3034.0, 3037.4, 3040.8, 3044.2, 3047.6};
  std::vector<std::vector<double>> inten{
    { 0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
      4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923}};
  inten.push_back(t);
  std::vector<double> tA(14, -1);
  std::vector<std::vector<double>> intenN{ { 1.2, 3.4, 5.6}};
  std::vector<double> tN = {3013.4, 3016.0, 3020.0};
  std::vector<int> flank{0,1,12,13};

  addFlankToRight(t, tN, tA, inten, intenN, flank);

  std::vector<int> cmp_fk = {0, 1};
  std::vector<double> cmp_tN = {3013.4, 3016, 3020, 3023.4, 3026.8};
  std::vector<double> cmp_tA = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3023.4, 3026.8};
  std::vector< std::vector< double > > cmp_arr;
  cmp_arr.push_back({1.2, 3.4, 5.6, 0.200958, 0.142092});
  //........................  CASE 1 ........................................
  for (int i = 0; i < intenN.size(); i++)
    for (int j = 0; j < intenN[i].size(); j++)
      ASSERT(std::abs(intenN[i][j] - cmp_arr[i][j]) < 1e-06);

  for (int j = 0; j < tN.size(); j++)
    ASSERT(std::abs(tN[j] - cmp_tN[j]) < 1e-06);

  for (int j = 0; j < tA.size(); j++)
    ASSERT(std::abs(tA[j] - cmp_tA[j]) < 1e-06);

  for (int j = 0; j < flank.size(); j++)
    ASSERT( flank[j] == cmp_fk[j]);
}

#ifdef DIALIGN_USE_Rcpp
int main_miscell(){
#else
  int main(){
#endif
    test_xicIntersect();
    test_interpolateZero();
    //test_getKeep();
    //test_getFlank();
    test_getFlankN();
    test_addFlankToLeft();
    test_addFlankToRight();
    std::cout << "test miscell successful" << std::endl;
    return 0;
  }
