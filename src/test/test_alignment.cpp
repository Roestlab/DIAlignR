#include <iostream>
#include <cmath> // require for std::abs
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "../alignment.h"
#include "../utils.h" //To propagate #define USE_Rcpp

#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.
using namespace DIAlign;
using namespace Alignment;

// Anonymous namespace: Only valid for this file.
namespace {
  std::vector< std::vector< TracebackType > > assertThisTraceback(int caseNum){
    std::vector< std::vector< TracebackType > > cmp_arr_Traceback;
    std::vector<TracebackType> tmp_tb;
    switch(caseNum){
    case 1: {
      tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, LM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, TM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      break;
    }
    case 2:{
      tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, LM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, TM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, TM, TM, DM, DM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
      break;
    }
    case 3: {
      tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      break;
    }
    case 4:{
      tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
      break;
    }
    }
    return cmp_arr_Traceback;
  }
  std::vector< std::vector< double > > assertThisM(int caseNum){
    std::vector< std::vector< double > > cmp_arr_M;
    std::vector<double> tmp;
    switch(caseNum){
    case 1:{
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      tmp = {0, -2, -2, 10, -2, 10}; cmp_arr_M.push_back(tmp);
      tmp = {0, 10, -4, -4, 8, -4}; cmp_arr_M.push_back(tmp);
      tmp = {0, -2, 20, -2, -6, 6}; cmp_arr_M.push_back(tmp);
      tmp = {0, -2, -2, 18, 8, -8}; cmp_arr_M.push_back(tmp);
      break;
    }
    case 2:{
      tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M.push_back(tmp);
      tmp = {-22, -2, -24, -34, -56, -78}; cmp_arr_M.push_back(tmp);
      tmp = {-44, -12, -4, -26, -36, -58}; cmp_arr_M.push_back(tmp);
      tmp = {-66, -34, -2, -6, -28, -38}; cmp_arr_M.push_back(tmp);
      tmp = {-88, -56, -24, -4, 4, -18}; cmp_arr_M.push_back(tmp);
      break;
    }
    case 3:{
      tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M.push_back(tmp);
      tmp = {-22, 0, -22, -44, -66, -88}; cmp_arr_M.push_back(tmp);
      tmp = {-44, -22, 0, -22, -44, -66}; cmp_arr_M.push_back(tmp);
      tmp = {-66, -44, -22, 0, -22, -44}; cmp_arr_M.push_back(tmp);
      tmp = {-88, -66, -44, -22, 0, -22}; cmp_arr_M.push_back(tmp);
      break;
    }
    case 4:{
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
      break;
    }
    }
    return cmp_arr_M;
  }
  // Following Path assertion for test_doAlignment ONLY.
  std::vector< std::vector< bool > > assertThisPath(int caseNum){
    std::vector< std::vector<bool> > cmp_arr_Path;
    std::vector<bool> tmp_b;
    switch(caseNum){
    case 1:{
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 2:{
      tmp_b = {false, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, false, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, false, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, false, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, false, false}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 3:{
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 4:{
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    }

    return cmp_arr_Path;
  }
  // Following Path assertion for test_getAlignedIndices ONLY.
  std::vector< std::vector< bool > > assertThisPath2(int caseNum){
    std::vector< std::vector<bool> > cmp_arr_Path;
    std::vector<bool> tmp_b;
    switch(caseNum){
    case 1:{
      tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {1, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 1, 1, 1}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 2:{
      tmp_b = {false, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, true, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, true, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, true, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, true, true}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 3:{
      tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 1, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 1, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 1}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    case 4:{
      tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 1, 0, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 1, 0}; cmp_arr_Path.push_back(tmp_b);
      tmp_b = {0, 0, 0, 0, 0, 1}; cmp_arr_Path.push_back(tmp_b);
      break;
    }
    }

    return cmp_arr_Path;
  }
  std::vector< std::vector< int > > assertThisOptionalPaths(int caseNum){
    std::vector< std::vector<int> > cmp_arr_OptionalPaths;
    std::vector<int> tmp_i;
    switch(caseNum){
    case 1:{
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      break;
    }
    case 2:{
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 2, 1, 1, 2}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 3, 1, 2}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 4, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      break;
    }
    case 3:{
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 1, 2, 3, 4, 5}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 2, 1, 3, 6, 10}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 3, 3, 1, 4, 10}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 4, 6, 4, 1, 5}; cmp_arr_OptionalPaths.push_back(tmp_i);
      break;
    }
    case 4:{
      tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 3, 5, 7, 9, 11}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 5, 13, 25, 41, 61}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 7, 25, 63, 129, 231}; cmp_arr_OptionalPaths.push_back(tmp_i);
      tmp_i = {1, 9, 41, 129, 321, 681}; cmp_arr_OptionalPaths.push_back(tmp_i);
      break;
    }
    }
    return cmp_arr_OptionalPaths;
  }
  std::vector< std::vector< double > > assertThisM_forw(int caseNum){
    std::vector< std::vector< double > > cmp_arr_M_forw;
    std::vector<double> tmp;
    switch(caseNum){
    case 1:{
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, -46, -136, -258, -436, -646}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, -124, -532, -1376, -2832, -5076}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, -258, -1304, -4338, -10884, -23054}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, -436, -2760, -10740, -31008, -77348}; cmp_arr_M_forw.push_back(tmp);
      break;
    }
    case 2:{
      tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-22, -90, -246, -478, -810, -1218}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-44, -234, -796, -1970, -4020, -7210}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-66, -478, -1898, -5790, -14118, -29610}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-88, -810, -3948, -13974, -38928, -95058}; cmp_arr_M_forw.push_back(tmp);
      break;
    }
    case 3:{
      tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-22, -88, -242, -484, -814, -1232}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-44, -242, -792, -1958, -4004, -7194}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-66, -484, -1958, -5808, -14058, -29436}; cmp_arr_M_forw.push_back(tmp);
      tmp = {-88, -814, -4004, -14058, -39600, -95238}; cmp_arr_M_forw.push_back(tmp);
      break;
    }
    case 4:{
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
      break;
    }
    }
    return cmp_arr_M_forw;
  }
  bool assertThisFreeEndGaps(int caseNum){
    bool tmp;
    switch(caseNum){
    case 1:
    case 4: {
      tmp = true;
      break;
    }
    case 2:
    case 3: {
      tmp = false;
      break;
    }
    }
    return tmp;
  }
  double assertThisGap(int caseNum){
    double tmp;
    switch(caseNum){
    case 1:
    case 2:
    case 3:{
      tmp = 22.0;
      break;
    }
    case 4: {
      tmp = 0.0;
      break;
    }
    }
    return tmp;
  }
  std::vector<int> assertThisindexA_aligned(int caseNum){
    std::vector<int> tmp;
    switch(caseNum){
    case 1:{
      tmp = {1, 2, 3, 4, 0, 0};
      break;
    }
    case 2:{
      tmp = {1, 2, 3, 4, 0};
      break;
    }
    case 3:{
      tmp = {0, 1, 2, 3, 4};
      break;
    }
    case 4: {
      tmp = {0, 1, 2, 3, 4};
      break;
    }
    }
    return tmp;
  }
  std::vector<int> assertThisindexB_aligned(int caseNum){
    std::vector<int> tmp;
    switch(caseNum){
    case 1:{
      tmp = {0, 1, 2, 3, 4, 5};
      break;
    }
    case 2:{
      tmp = {1, 2, 3, 4, 5};
      break;
    }
    case 3:{
      tmp = {1, 2, 3, 4, 5};
      break;
    }
    case 4: {
      tmp = {1, 2, 3, 4, 5};
      break;
    }
    }
    return tmp;
  }
  std::vector<double> assertThisscore(int caseNum){
    std::vector<double> tmp;
    switch(caseNum){
    case 1:{
      tmp = {0, 10, 20, 18, 18, 18};
      break;
    }
    case 2:{
      tmp = {-2.0, -4.0, -6.0, 4.0, -18};
      break;
    }
    case 3:{
      tmp = {-22.0, -22.0, -22.0, -22.0, -22.0};
      break;
    }
    case 4: {
      tmp = {0.0, 0.0, 0.0, 0.0, 0.0};
      break;
    }
    }
    return tmp;
  }
  double assertThisscore_forw(int caseNum){
    double tmp;
    switch(caseNum){
    case 1:{
      tmp = -10740.0;
      break;
    }
    case 2:{
      tmp = -95058.0;
      break;
    }
    case 3:{
      tmp = -95238.0;
      break;
    }
    case 4: {
      tmp = 0.0;
      break;
    }
    }
    return tmp;
  }

  template<class T>
  std::vector<T> vov2v(std::vector< std::vector< T > > vov){
    std::vector<T> vec;
    for (const auto& v : vov) vec.insert(vec.end(), v.begin(), v.end());
    return vec;
  }
  void fillAlignObj(AlignObj& obj, int caseNum){
    obj.Traceback = vov2v(assertThisTraceback(caseNum));
    obj.M = vov2v(assertThisM(caseNum));
    obj.Path = vov2v(assertThisPath(caseNum));
    obj.OptionalPaths = vov2v(assertThisOptionalPaths(caseNum));
    obj.M_forw = vov2v(assertThisM_forw(caseNum));
    obj.signalA_len = 4;
    obj.signalB_len = 5;
    obj.GapOpen = assertThisGap(caseNum);
    obj.GapExten = assertThisGap(caseNum);
    obj.FreeEndGaps = assertThisFreeEndGaps(caseNum);
  }
}

void test_doAlignment_cases(AlignObj obj, int caseNum){
  // Traceback
  std::vector< std::vector< TracebackType > > cmp_arr_Traceback = assertThisTraceback(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback[i][j]);
    }
  }
  // M
  std::vector< std::vector< double > > cmp_arr_M = assertThisM(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M[i][j]) < 1e-06);
    }
  }
  // Path
  std::vector< std::vector< bool > > cmp_arr_Path = assertThisPath(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path[i][j]);
    }
  }
  // OptionalPaths
  std::vector< std::vector< int > > cmp_arr_OptionalPaths = assertThisOptionalPaths(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.OptionalPaths[i*6+j] == cmp_arr_OptionalPaths[i][j]);
    }
  }
  // M_forw
  std::vector< std::vector< double > > cmp_arr_M_forw = assertThisM_forw(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M_forw[i*6+j] - cmp_arr_M_forw[i][j]) < 1e-06);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - assertThisGap(caseNum)) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - assertThisGap(caseNum)) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == assertThisFreeEndGaps(caseNum));
  // indexA_aligned
  ASSERT(obj.indexA_aligned.size() == 0);
  // indexB_aligned
  ASSERT(obj.indexB_aligned.size() == 0);
  // score
  ASSERT(obj.score.size() == 0);
  // score_forw
  ASSERT(std::abs(obj.score_forw - 0.0) < 1e-6);
}

void test_doAlignment(){
  SimMatrix s;
  s.data = {-2, -2, 10, -2, 10,
            10, -2, -2, -2, -2,
            -2, 10, -2, -2, -2,
            -2, -2, -2, 10, -2};
  s.n_col = 5;
  s.n_row = 4;
  double gap = 22.0;

  //........................  CASE 1 ........................................
  int caseNum = 1;
  AlignObj obj = doAlignment(s, gap, true);
  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  test_doAlignment_cases(obj, caseNum);

  //........................  CASE 2 ........................................
  caseNum = 2;
  obj = doAlignment(s, gap, false);
  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  test_doAlignment_cases(obj, caseNum);

  //........................  CASE 3 ........................................
  caseNum = 3;
  s.data = {0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0};
  s.n_col = 5;
  s.n_row = 4;
  gap = 22.0;
  obj = doAlignment(s, gap, false);
  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  test_doAlignment_cases(obj, caseNum);

  //........................  CASE 4 ........................................
  caseNum = 4;
  s.data = {0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0};
  s.n_col = 5;
  s.n_row = 4;
  gap = 0.0;
  obj = doAlignment(s, gap, true);
  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  test_doAlignment_cases(obj, caseNum);
}

void test_getAlignedIndices_cases(AlignObj obj, int caseNum){
  // Traceback
  std::vector< std::vector< TracebackType > > cmp_arr_Traceback = assertThisTraceback(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback[i][j]);
    }
  }
  // M
  std::vector< std::vector< double > > cmp_arr_M = assertThisM(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M[i][j]) < 1e-06);
    }
  }
  // Path
  std::vector< std::vector< bool > > cmp_arr_Path = assertThisPath2(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path[i][j]);
    }
  }
  // OptionalPaths
  std::vector< std::vector< int > > cmp_arr_OptionalPaths = assertThisOptionalPaths(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.OptionalPaths[i*6+j] == cmp_arr_OptionalPaths[i][j]);
    }
  }
  // M_forw
  std::vector< std::vector< double > > cmp_arr_M_forw = assertThisM_forw(caseNum);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M_forw[i*6+j] - cmp_arr_M_forw[i][j]) < 1e-06);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - assertThisGap(caseNum)) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - assertThisGap(caseNum)) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == assertThisFreeEndGaps(caseNum));
  // indexA_aligned
  std::vector<int> cmp_indexA_aligned = assertThisindexA_aligned(caseNum);
  for(int i = 0; i < obj.indexA_aligned.size(); i++){
    ASSERT(obj.indexA_aligned[i] == cmp_indexA_aligned[i]);
  }
  // indexB_aligned
  std::vector<int> cmp_indexB_aligned = assertThisindexB_aligned(caseNum);
  for(int i = 0; i < obj.indexB_aligned.size(); i++){
    ASSERT(obj.indexB_aligned[i] == cmp_indexB_aligned[i]);
  }
  // score
  std::vector<double> cmp_score = assertThisscore(caseNum);
  for(int i = 0; i < obj.score.size(); i++){
    ASSERT(std::abs(obj.score[i] - cmp_score[i]) < 1e-06);
  }
  // score_forw
  ASSERT(std::abs(obj.score_forw - assertThisscore_forw(caseNum)) < 1e-6);
}

void test_getAlignedIndices(){

  //........................  CASE 1 ........................................
  int caseNum = 1;

  AlignObj obj(4,5);
  obj.s_data = {-2, -2, 10, -2, 10,
                10, -2, -2, -2, -2,
                -2, 10, -2, -2, -2,
                -2, -2, -2, 10, -2};

  fillAlignObj(obj, caseNum);
  getAlignedIndices(obj);

  std::vector< std::vector< double > > cmp_arr_s;
  std::vector<double> tmp;
  tmp = {-2, -2, 10, -2, 10}; cmp_arr_s.push_back(tmp);
  tmp = {10, -2, -2, -2, -2}; cmp_arr_s.push_back(tmp);
  tmp = {-2, 10, -2, -2, -2}; cmp_arr_s.push_back(tmp);
  tmp = {-2, -2, -2, 10, -2}; cmp_arr_s.push_back(tmp);
  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - cmp_arr_s[i][j]) < 1e-06);
    }
  }
  test_getAlignedIndices_cases(obj, caseNum);

  //........................  CASE 2 ........................................
  caseNum = 2;
  fillAlignObj(obj, caseNum);
  getAlignedIndices(obj);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - cmp_arr_s[i][j]) < 1e-06);
    }
  }
  test_getAlignedIndices_cases(obj, caseNum);

  //........................  CASE 3 ........................................
  caseNum = 3;
  fillAlignObj(obj, caseNum);
  obj.s_data = {0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0};
  getAlignedIndices(obj);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - 0.0) < 1e-06);
    }
  }
  test_getAlignedIndices_cases(obj, caseNum);

  //........................  CASE 4 ........................................
  caseNum = 4;
  fillAlignObj(obj, caseNum);
  obj.s_data = {0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0};
  getAlignedIndices(obj);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - 0.0) < 1e-06);
    }
  }
  test_getAlignedIndices_cases(obj, caseNum);
}

#ifdef DIALIGN_USE_Rcpp
int main_aligment(){
#else
int main(){
#endif
  test_doAlignment();
  test_getAlignedIndices();
  std::cout << "test alignment successful" << std::endl;
  return 0;
}
