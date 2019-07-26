#include <iostream>
#include <cmath> // require for std::abs
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "alignment.h"
#include "utils.h" //To propagate #define USE_Rcpp

#define ASSERT(condition) if(!(condition)) throw 1; // If you don't put the message, C++ will output the code.
using namespace DIAlign;

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
  AlignObj obj = doAlignment(s, gap, true);

  std::vector< std::vector< TracebackType > > cmp_arr_Traceback;
  std::vector<TracebackType> tmp_tb;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, LM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, TM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M;
  std::vector<double> tmp;
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, -2, 10, -2, 10}; cmp_arr_M.push_back(tmp);
  tmp = {0, 10, -4, -4, 8, -4}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, 20, -2, -6, 6}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, -2, 18, 8, -8}; cmp_arr_M.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path;
  std::vector<bool> tmp_b;
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);

  std::vector< std::vector<int> > cmp_arr_OptionalPaths;
  std::vector<int> tmp_i;
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths.push_back(tmp_i);

  std::vector< std::vector< double > > cmp_arr_M_forw;
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M_forw.push_back(tmp);
  tmp = {0, -46, -48, -36, -36, -36}; cmp_arr_M_forw.push_back(tmp);
  tmp = {0, -36, -40, -42, -42, -30}; cmp_arr_M_forw.push_back(tmp);
  tmp = {0, -36, -30, -34, -44, -48}; cmp_arr_M_forw.push_back(tmp);
  tmp = {0, -48, -30, -30, -24, -38}; cmp_arr_M_forw.push_back(tmp);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path[i][j]);
    }
  }
  // OptionalPaths
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.OptionalPaths[i*6+j] == cmp_arr_OptionalPaths[i][j]);
    }
  }
  // M_forw
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
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == true);
  // indexA_aligned
  ASSERT(obj.indexA_aligned.size() == 0);
  // indexB_aligned
  ASSERT(obj.indexB_aligned.size() == 0);
  // score
  ASSERT(obj.score.size() == 0);

  //........................  CASE 2 ........................................
  obj = doAlignment(s, gap, false);
  std::vector< std::vector< TracebackType > > cmp_arr_Traceback2;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, LM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, TM, DM, DM, DM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, TM, TM, DM, DM, LM}; cmp_arr_Traceback2.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M2;
  tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M2.push_back(tmp);
  tmp = {-22, -2, -24, -34, -56, -78}; cmp_arr_M2.push_back(tmp);
  tmp = {-44, -12, -4, -26, -36, -58}; cmp_arr_M2.push_back(tmp);
  tmp = {-66, -34, -2, -6, -28, -38}; cmp_arr_M2.push_back(tmp);
  tmp = {-88, -56, -24, -4, 4, -18}; cmp_arr_M2.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path2;
  tmp_b = {false, 0, 0, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, false, 0, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, false, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, 0, false, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, false, false}; cmp_arr_Path2.push_back(tmp_b);

  std::vector< std::vector<int> > cmp_arr_OptionalPaths2;
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths2.push_back(tmp_i);
  tmp_i = {1, 1, 2, 1, 1, 2}; cmp_arr_OptionalPaths2.push_back(tmp_i);
  tmp_i = {1, 1, 1, 3, 1, 2}; cmp_arr_OptionalPaths2.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 4, 1}; cmp_arr_OptionalPaths2.push_back(tmp_i);
  tmp_i = {1, 1, 1, 1, 1, 1}; cmp_arr_OptionalPaths2.push_back(tmp_i);

  std::vector< std::vector< double > > cmp_arr_M_forw2;
  tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M_forw2.push_back(tmp);
  tmp = {-22, -90, -114, -168, -234, -288}; cmp_arr_M_forw2.push_back(tmp);
  tmp = {-44, -102, -84, -108, -162, -216}; cmp_arr_M_forw2.push_back(tmp);
  tmp = {-66, -168, -84, -78, -114, -168}; cmp_arr_M_forw2.push_back(tmp);
  tmp = {-88, -234, -138, -78, -72, -108}; cmp_arr_M_forw2.push_back(tmp);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback2[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M2[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path2[i][j]);
    }
  }
  // OptionalPaths
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.OptionalPaths[i*6+j] == cmp_arr_OptionalPaths2[i][j]);
    }
  }
  // M_forw
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M_forw[i*6+j] - cmp_arr_M_forw2[i][j]) < 1e-06);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == false);
  // indexA_aligned
  ASSERT(obj.indexA_aligned.size() == 0);
  // indexB_aligned
  ASSERT(obj.indexB_aligned.size() == 0);
  // score
  ASSERT(obj.score.size() == 0);


  //........................  CASE 3 ........................................
  s.data = {0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0};
  s.n_col = 5;
  s.n_row = 4;
  gap = 22.0;
  obj = doAlignment(s, gap, false);
  std::vector< std::vector< TracebackType > > cmp_arr_Traceback3;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M3;
  tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M3.push_back(tmp);
  tmp = {-22, 0, -22, -44, -66, -88}; cmp_arr_M3.push_back(tmp);
  tmp = {-44, -22, 0, -22, -44, -66}; cmp_arr_M3.push_back(tmp);
  tmp = {-66, -44, -22, 0, -22, -44}; cmp_arr_M3.push_back(tmp);
  tmp = {-88, -66, -44, -22, 0, -22}; cmp_arr_M3.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path3;
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback3[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M3[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path3[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == false);
  // indexA_aligned
  ASSERT(obj.indexA_aligned.size() == 0);
  // indexB_aligned
  ASSERT(obj.indexB_aligned.size() == 0);
  // score
  ASSERT(obj.score.size() == 0);

  //........................  CASE 4 ........................................
  s.data = {0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0};
  s.n_col = 5;
  s.n_row = 4;
  gap = 0.0;
  obj = doAlignment(s, gap, true);
  std::vector< std::vector< TracebackType > > cmp_arr_Traceback4;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M4;
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path4;
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - s.data[i*5+j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback4[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M4[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path4[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 0.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 0.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == true);
  // indexA_aligned
  ASSERT(obj.indexA_aligned.size() == 0);
  // indexB_aligned
  ASSERT(obj.indexB_aligned.size() == 0);
  // score
  ASSERT(obj.score.size() == 0);
}

void test_getAlignedIndices(){

  AlignObj obj(4,5);
  obj.s_data = {-2, -2, 10, -2, 10,
                10, -2, -2, -2, -2,
                -2, 10, -2, -2, -2,
                -2, -2, -2, 10, -2};
  obj.Traceback = {SS, LM, LM, LM, LM, LM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, LM, DM, DM,
                   TM, DM, TM, DM, DM, DM};
  obj.M = {0, 0, 0, 0, 0, 0,
           0, -2, -2, 10, -2, 10,
           0, 10, -4, -4, 8, -4,
           0, -2, 20, -2, -6, 6,
           0, -2, -2, 18, 8, -8};
  obj.Path = {0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0};
  obj.signalA_len = 4;
  obj.signalB_len = 5;
  obj.GapOpen = 22.0;
  obj.GapExten = 22.0;
  obj.FreeEndGaps = true;

  getAlignedIndices(obj);

  std::vector< std::vector< double > > cmp_arr_s;
  std::vector<double> tmp;
  tmp = {-2, -2, 10, -2, 10}; cmp_arr_s.push_back(tmp);
  tmp = {10, -2, -2, -2, -2}; cmp_arr_s.push_back(tmp);
  tmp = {-2, 10, -2, -2, -2}; cmp_arr_s.push_back(tmp);
  tmp = {-2, -2, -2, 10, -2}; cmp_arr_s.push_back(tmp);

  std::vector< std::vector< TracebackType > > cmp_arr_Traceback;
  std::vector<TracebackType> tmp_tb;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, LM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);
  tmp_tb = {TM, DM, TM, DM, DM, DM}; cmp_arr_Traceback.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M;
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, -2, 10, -2, 10}; cmp_arr_M.push_back(tmp);
  tmp = {0, 10, -4, -4, 8, -4}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, 20, -2, -6, 6}; cmp_arr_M.push_back(tmp);
  tmp = {0, -2, -2, 18, 8, -8}; cmp_arr_M.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path;
  std::vector<bool> tmp_b;
  tmp_b = {0, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {1, 0, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path.push_back(tmp_b);
  tmp_b = {0, 0, 0, 1, 1, 1}; cmp_arr_Path.push_back(tmp_b);
  std::vector<int> cmp_indexA_aligned = {1, 2, 3, 4, 0, 0};
  std::vector<int> cmp_indexB_aligned = {0, 1, 2, 3, 4, 5};
  std::vector<double> cmp_score = {0, 10, 20, 18, 18, 18};

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - cmp_arr_s[i][j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == true);
  // indexA_aligned
  for(int i = 0; i < obj.indexA_aligned.size(); i++){
    ASSERT(obj.indexA_aligned[i] == cmp_indexA_aligned[i]);
  }
  // indexB_aligned
  for(int i = 0; i < obj.indexB_aligned.size(); i++){
    ASSERT(obj.indexB_aligned[i] == cmp_indexB_aligned[i]);
  }
  // score
  for(int i = 0; i < obj.score.size(); i++){
    ASSERT(std::abs(obj.score[i] - cmp_score[i]) < 1e-06);
  }

  //........................  CASE 2 ........................................

  obj.Traceback = {SS, LM, LM, LM, LM, LM,
                   TM, DM, DM, DM, LM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, TM, DM, DM, DM, DM,
                   TM, TM, TM, DM, DM, LM};
  obj.M = {0, -22, -44, -66, -88, -110,
           -22, -2, -24, -34, -56, -78,
           -44, -12, -4, -26, -36, -58,
           -66, -34, -2, -6, -28, -38,
           -88, -56, -24, -4, 4, -18};
  obj.Path = {0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0};
  obj.signalA_len = 4;
  obj.signalB_len = 5;
  obj.GapOpen = 22.0;
  obj.GapExten = 22.0;
  obj.FreeEndGaps = false;

  getAlignedIndices(obj);

  std::vector< std::vector< TracebackType > > cmp_arr_Traceback2;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, LM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, TM, DM, DM, DM, DM}; cmp_arr_Traceback2.push_back(tmp_tb);
  tmp_tb = {TM, TM, TM, DM, DM, LM}; cmp_arr_Traceback2.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M2;
  tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M2.push_back(tmp);
  tmp = {-22, -2, -24, -34, -56, -78}; cmp_arr_M2.push_back(tmp);
  tmp = {-44, -12, -4, -26, -36, -58}; cmp_arr_M2.push_back(tmp);
  tmp = {-66, -34, -2, -6, -28, -38}; cmp_arr_M2.push_back(tmp);
  tmp = {-88, -56, -24, -4, 4, -18}; cmp_arr_M2.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path2;
  tmp_b = {false, 0, 0, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, true, 0, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, true, 0, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, 0, true, 0, 0}; cmp_arr_Path2.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, true, true}; cmp_arr_Path2.push_back(tmp_b);
  cmp_indexA_aligned = {1, 2, 3, 4, 0};
  cmp_indexB_aligned = {1, 2, 3, 4, 5};
  cmp_score = {-2.0, -4.0, -6.0, 4.0, -18};

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - cmp_arr_s[i][j]) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback2[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M2[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path2[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == false);
  // indexA_aligned
  for(int i = 0; i < obj.indexA_aligned.size(); i++){
    ASSERT(obj.indexA_aligned[i] == cmp_indexA_aligned[i]);
  }
  // indexB_aligned
  for(int i = 0; i < obj.indexB_aligned.size(); i++){
    ASSERT(obj.indexB_aligned[i] == cmp_indexB_aligned[i]);
  }
  // score
  for(int i = 0; i < obj.score.size(); i++){
    ASSERT(std::abs(obj.score[i] - cmp_score[i]) < 1e-06);
  }

  //........................  CASE 3 ........................................
  obj.s_data = {0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0};
  obj.Traceback = {SS, LM, LM, LM, LM, LM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM};
  obj.M = {0, -22, -44, -66, -88, -110,
           -22, 0, -22, -44, -66, -88,
           -44, -22, 0, -22, -44, -66,
           -66, -44, -22, 0, -22, -44,
           -88, -66, -44, -22, 0, -22};
  obj.Path = {0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0};
  obj.signalA_len = 4;
  obj.signalB_len = 5;
  obj.GapOpen = 22.0;
  obj.GapExten = 22.0;
  obj.FreeEndGaps = false;

  getAlignedIndices(obj);

  std::vector< std::vector< TracebackType > > cmp_arr_Traceback3;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback3.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M3;
  tmp = {0, -22, -44, -66, -88, -110}; cmp_arr_M3.push_back(tmp);
  tmp = {-22, 0, -22, -44, -66, -88}; cmp_arr_M3.push_back(tmp);
  tmp = {-44, -22, 0, -22, -44, -66}; cmp_arr_M3.push_back(tmp);
  tmp = {-66, -44, -22, 0, -22, -44}; cmp_arr_M3.push_back(tmp);
  tmp = {-88, -66, -44, -22, 0, -22}; cmp_arr_M3.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path3;
  tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 1, 0, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 1, 0}; cmp_arr_Path3.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 1}; cmp_arr_Path3.push_back(tmp_b);
  cmp_indexA_aligned = {0, 1, 2, 3, 4};
  cmp_indexB_aligned = {1, 2, 3, 4, 5};
  cmp_score = {-22.0, -22.0, -22.0, -22.0, -22.0};

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - 0.0) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback3[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M3[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path3[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 22.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 22.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == false);
  // indexA_aligned
  for(int i = 0; i < obj.indexA_aligned.size(); i++){
    ASSERT(obj.indexA_aligned[i] == cmp_indexA_aligned[i]);
  }
  // indexB_aligned
  for(int i = 0; i < obj.indexB_aligned.size(); i++){
    ASSERT(obj.indexB_aligned[i] == cmp_indexB_aligned[i]);
  }
  // score
  for(int i = 0; i < obj.score.size(); i++){
    ASSERT(std::abs(obj.score[i] - cmp_score[i]) < 1e-06);
  }

  //........................  CASE 4 ........................................
  obj.s_data = {0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0};
  obj.Traceback = {SS, LM, LM, LM, LM, LM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM,
                   TM, DM, DM, DM, DM, DM};
  obj.M = {0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0};
  obj.Path = {0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0};
  obj.signalA_len = 4;
  obj.signalB_len = 5;
  obj.GapOpen = 0.0;
  obj.GapExten = 0.0;
  obj.FreeEndGaps = true;

  getAlignedIndices(obj);

  std::vector< std::vector< TracebackType > > cmp_arr_Traceback4;
  tmp_tb = {SS, LM, LM, LM, LM, LM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);
  tmp_tb = {TM, DM, DM, DM, DM, DM}; cmp_arr_Traceback4.push_back(tmp_tb);

  std::vector< std::vector< double > > cmp_arr_M4;
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);
  tmp = {0, 0, 0, 0, 0, 0}; cmp_arr_M4.push_back(tmp);

  std::vector< std::vector<bool> > cmp_arr_Path4;
  tmp_b = {0, 1, 0, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 1, 0, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 1, 0, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 1, 0}; cmp_arr_Path4.push_back(tmp_b);
  tmp_b = {0, 0, 0, 0, 0, 1}; cmp_arr_Path4.push_back(tmp_b);
  cmp_indexA_aligned = {0, 1, 2, 3, 4};
  cmp_indexB_aligned = {1, 2, 3, 4, 5};
  cmp_score = {0.0, 0.0, 0.0, 0.0, 0.0};

  // s_data
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 5; j++){
      ASSERT(std::abs(obj.s_data[i*5+j] - 0.0) < 1e-06);
    }
  }
  // Traceback
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Traceback[i*6+j] == cmp_arr_Traceback4[i][j]);
    }
  }
  // M
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(std::abs(obj.M[i*6+j] - cmp_arr_M4[i][j]) < 1e-06);
    }
  }
  // Path
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 6; j++){
      ASSERT(obj.Path[i*6+j] == cmp_arr_Path4[i][j]);
    }
  }
  // signalA_len
  ASSERT(obj.signalA_len == 4);
  // signalB_len
  ASSERT(obj.signalB_len == 5);
  // GapOpen
  ASSERT(std::abs(obj.GapOpen - 0.0) < 1e-6);
  // GapExten
  ASSERT(std::abs(obj.GapExten - 0.0) < 1e-6);
  // FreeEndGaps
  ASSERT(obj.FreeEndGaps == true);
  // indexA_aligned
  for(int i = 0; i < obj.indexA_aligned.size(); i++){
    ASSERT(obj.indexA_aligned[i] == cmp_indexA_aligned[i]);
  }
  // indexB_aligned
  for(int i = 0; i < obj.indexB_aligned.size(); i++){
    ASSERT(obj.indexB_aligned[i] == cmp_indexB_aligned[i]);
  }
  // score
  for(int i = 0; i < obj.score.size(); i++){
    ASSERT(std::abs(obj.score[i] - cmp_score[i]) < 1e-06);
  }
}

#ifdef USE_Rcpp
int main_aligment(){
#else
int main(){
#endif
  test_doAlignment();
  test_getAlignedIndices();
  std::cout << "test alignment successful" << std::endl;
  return 0;
}
