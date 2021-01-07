#include "spline.h"
namespace DIAlign
{

std::vector<double> naturalSpline(const std::vector<double> & x, const std::vector<double> & y,
                                  const std::vector<double> &  xout){
  tk::spline s;
  s.set_points(x, y);
  std::vector<double> result(xout.size());
  std::transform (xout.begin(), xout.end(), result.begin(), s);
  return result;

}
}
