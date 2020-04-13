#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "../integrateArea.h"
#include "../utils.h" //To propagate #define USE_Rcpp

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;
using namespace PeakGroupIntensity;

void test_peakGroupArea(){
  const double left = 2.472833334;
  const double right = 3.022891666;

  const std::vector<double> position = {
    2.23095,2.239716667,2.248866667,2.25765,2.266416667,
    2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
    2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
    2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
    2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
    2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
    2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
    2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
    2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
    2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
    2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
    3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
    3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
    3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
    3.190133333,3.198016667,3.205916667,3.213166667
  };

  const std::vector<double> intensity = {
    1447,2139,1699,755,1258,1070,944,1258,1573,1636,
    1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
    1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
    51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
    184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
    6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
    2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
    2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
  };

  //........................  CASE 1 ........................................
  std::vector< std::vector< double > > position_arr;
  position_arr.push_back(position);

  std::vector< std::vector< double > > intensity_arr;
  intensity_arr.push_back(intensity);

  double area;
  area = peakGroupArea(position_arr, intensity_arr, left, right, "intensity_sum", "base_to_base", false);
  ASSERT(std::abs(area - 123446.661339019) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "trapezoid", "base_to_base", false);
  ASSERT(std::abs(area - 1140.392865964) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "simpson", "base_to_base", false);
  //ASSERT(std::abs(area - 123446.661339019) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "intensity_sum", "vertical_division_min", false);
  ASSERT(std::abs(area - 50217) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "trapezoid", "vertical_division_min", false);
  ASSERT(std::abs(area - 476.606316373) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "simpson", "vertical_division_min", false);
  //ASSERT(std::abs(area - 123446.661339019) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "intensity_sum", "vertical_division_max", false);
  ASSERT(std::abs(area - 190095) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "trapezoid", "vertical_division_max", false);
  ASSERT(std::abs(area - 1804.179415555) < 1e-03);

  area = peakGroupArea(position_arr, intensity_arr, left, right, "simpson", "vertical_division_max", false);
  //ASSERT(std::abs(area - 123446.661339019) < 1e-03);

  //........................  CASE 2 ........................................

  //........................  CASE 3 ........................................

}

#ifdef DIALIGN_USE_Rcpp
int main_peakGroupArea(){
#else
  int main(){
#endif
    test_peakGroupArea();
    std::cout << "test peakGroupArea successful" << std::endl;
    return 0;
  }
