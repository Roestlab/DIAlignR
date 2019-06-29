#include <vector>
#include <cmath> // require for std::abs
#include <assert.h>
#include "affinealignobj.h"
#include "utils.h"

//TODO update this statement so we know which line failed.
#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.
using namespace DIAlign;

void test_EnumToChar(){
  std::vector<TracebackType> tb = {DA, SS, TM, DM, TA, DB, TB, TM, LA, SS, LB, LM};
  std::vector<char> ch = EnumToChar(tb);
  std::vector<char> cmp_vec = {'2', '0', '4', '1', '5', '3', '6', '4', '8', '0', '9', '7'};
  // character literals are in single quotes only.

  for(int i = 0; i < ch.size(); i++){
    ASSERT(ch[i] == cmp_vec[i]);
  }
}

#ifdef USE_Rcpp
int main_affinealignobj(){
#else
int main(){
#endif
  test_EnumToChar();
  std::cout << "test affinealignobj successful" << std::endl;
  return 0;
}
