#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "../miscell.h"
#include "../utils.h" //To propagate #define USE_Rcpp

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;

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

#ifdef DIALIGN_USE_Rcpp
int main_xicIntersect(){
#else
  int main(){
#endif
    test_xicIntersect();
    test_interpolateZero();
    std::cout << "test miscell successful" << std::endl;
    return 0;
  }
