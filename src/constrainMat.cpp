#include "constrainMat.h"


namespace DIAlign 
{

// TODO Use ascii art
// TODO Make sure A1 < A2, B1 < B2 and B1p < B2p
// TODO Make sure same deltaTime

void calcNoBeefMask(SimMatrix& MASK, double A1, double A2, double B1, double B2, double B1p, double B2p, int noBeef, bool hardConstrain){
  // This equation is to from parallel lines: slope = (A2-A1)/(B2-B1) =
  // (y-A1)/(x- B1). Also, we are considering a very simple linear fit with
  // slope ~= 1. This will allow a fitting line to cross square only at
  // adjacent sides. If the slope comes out to be deviating too much from 1,
  // then this equation will fail and also in that case better to fit
  // non-linear boundaries for No-Beef region.
  /***
   * mmB1pmmmmmmmB2pmmm
   * A1mmmmmmmmmmmmmmmm
   * mmmmmmmmmmmmmmmmmm
   * mmmmmmmmmmmmmmmmmm
   * mmmmmmmmmmmmmmmmmm
   * A2mmmmmmmmmmmmmmmm
   */
  double deltaTime = (A2-A1)/(MASK.n_row-1);
  int idx = floor((B1p - B1)/deltaTime) + 1; // Index of tB which will correspond to B1p.
  int leftTri1[2] = {idx-noBeef, 0}; // (X,Y) coordinate of imaginary left triangle.
  int rightTri1[2] = {idx+noBeef, 0};
  idx = ceil((B2p - B1)/deltaTime); // Index of tB which will correspond to B2p.
  int leftTri2[2] = {idx-noBeef, MASK.n_row-1};
  int rightTri2[2] = {idx+noBeef, MASK.n_row-1}; // Considering our square starts at (0,0)
  double distFromLeftBndry = 0.0, distFromRightBndry = 0.0;
  double lessThanLeftTri = 0.0, greaterThanRightTri = 0.0;
  if(hardConstrain){
    for(int x = 0; x < MASK.n_col; x++){
      for(int y = 0; y < MASK.n_row; y++){
        lessThanLeftTri = y-leftTri1[1] - (leftTri2[1]-leftTri1[1])*(x-leftTri1[0])/(leftTri2[0]-leftTri1[0]);
        greaterThanRightTri = y-rightTri1[1] - (rightTri2[1]-rightTri1[1])*(x-rightTri1[0])/(rightTri2[0]-rightTri1[0]);
        if(lessThanLeftTri<=0 & greaterThanRightTri>=0)
          MASK.data[y*MASK.n_col + x] = 0.0;
        else
          MASK.data[y*MASK.n_col + x] = 1.0;
      }
    }
  }
  else {
    for(int x = 0; x < MASK.n_col; x++){
      for(int y = 0; y < MASK.n_row; y++){
        // distance of point (X0, Y0) from a line that passes through P1 and P2.
        // distance(P1, P2, (X0, Y0)) <- -((Y2-Y1)X0 - (X2-X1)Y0 + X2*Y1 -
        // Y2*X1)/sqrt((Y2-Y1)^2 +(X2-X1)^2)
        distFromLeftBndry = -((leftTri2[1]-leftTri1[1])*x - (leftTri2[0]-leftTri1[0])*y + leftTri2[0]*leftTri1[1] -
          leftTri2[1]*leftTri1[0])/sqrt((leftTri2[1]-leftTri1[1])*(leftTri2[1]-leftTri1[1]) + (leftTri2[0]-leftTri1[0])*(leftTri2[0]-leftTri1[0]));
        distFromRightBndry = -((rightTri2[1]-rightTri1[1])*x - (rightTri2[0]-rightTri1[0])*y + rightTri2[0]*rightTri1[1] -
          rightTri2[1]*rightTri1[0])/sqrt((rightTri2[1]-rightTri1[1])*(rightTri2[1]-rightTri1[1]) + (rightTri2[0]-rightTri1[0])*(rightTri2[0]-rightTri1[0]));
        if(distFromLeftBndry<=0 & distFromRightBndry>=0)
          MASK.data[y*MASK.n_col + x] = 0.0;
        else if(distFromLeftBndry > 0.0)
          MASK.data[y*MASK.n_col + x] = distFromLeftBndry;
        else if(distFromRightBndry < 0.0)
          MASK.data[y*MASK.n_col + x] = -distFromRightBndry;
      }
    }
  }
}

// TODO Make sure s and MASK have same size
void constrainSimilarity(SimMatrix& s, const SimMatrix& MASK, double constrainVal){
  for(int i = 0; i < s.data.size(); i++)
    s.data[i] += constrainVal*MASK.data[i];
}
} // namespace DIAlign
