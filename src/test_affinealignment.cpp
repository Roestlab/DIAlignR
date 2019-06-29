#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "affinealignment.h"

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.

using namespace DIAlign;

void test_getOlapAffineAlignStartIndices(){

}

void test_getAffineAlignedIndices(){

}

void test_doAffineAlignment(){

}

#ifdef USE_Rcpp
int main_affinealignment(){
#else
int main(){
#endif
  test_getOlapAffineAlignStartIndices();
  test_getAffineAlignedIndices();
  test_doAffineAlignment();
  std::cout << "test affinealignment successful" << std::endl;
  return 0;
}
