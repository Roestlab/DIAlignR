#include <iostream>
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "alignment.h"

void test_doAlignment(){
  AlignObj obj(3,4);
  std::cout << "test alignment successfu1l" << std::endl;
}

#ifdef USE_Rcpp
int main_alignment(){
#else
int main(){
#endif
  test_doAlignment();
  std::cout << "test alignment successful" << std::endl;
  return 0;
}
