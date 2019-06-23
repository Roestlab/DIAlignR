#include <cmath> // require for std::abs
#include <assert.h>
#include "constrainMat.h"

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

void test_calcNoBeefMask(){

}

void test_constrainSimilarity(){

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
