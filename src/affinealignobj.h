#ifndef AFFINEALIGNOBJ_H
#define AFFINEALIGNOBJ_H

#include <iostream>
#include <cstring>

enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};

std::ostream& operator<<(std::ostream& out, const TracebackType value);

enum tbJump {M = 0, A = 1, B = 2};

struct AffineAlignObj
{
    float* M;
    float* A;
    float* B;
    TracebackType* Traceback;
    int signalA_len; // stack allocation
    int signalB_len;
    float GapOpen;
    float GapExten;
    bool FreeEndGaps;

    // Not a default constructor
    AffineAlignObj(int ROW_SIZE, int COL_SIZE)
    {
        M = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        A = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        B = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        Traceback = new TracebackType[3*ROW_SIZE * COL_SIZE]; // heap allocation
        signalA_len = ROW_SIZE-1;
        signalB_len = COL_SIZE-1;
    }

    // Rule 1 Copy constructor
    AffineAlignObj(const AffineAlignObj &other)
    {
          signalA_len = other.signalA_len;
          signalB_len = other.signalB_len;
          GapOpen = other.GapOpen;
          GapExten = other.GapExten;
          FreeEndGaps = other.FreeEndGaps;
          M = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          A = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(A, other.A, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          B = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(B, other.B, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          Traceback = new TracebackType[3*(signalA_len+1)*(signalB_len+1)];
          std::memcpy(Traceback, other.Traceback, sizeof(TracebackType) * 3 * (signalA_len+1)*(signalB_len+1));
       }

    // Rule 2 Copy assignment operator
    AffineAlignObj& operator=(const AffineAlignObj& other)
    {
        if(this == &other) return *this; // handling of self assignment.
        delete[] M; // freeing previously used memory
        delete[] A;
        delete[] B;
        delete[] Traceback;
        signalA_len = other.signalA_len;
        signalB_len = other.signalB_len;
        GapOpen = other.GapOpen;
        GapExten = other.GapExten;
        FreeEndGaps = other.FreeEndGaps;
        M = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        A = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(A, other.A, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        B = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(B, other.B, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        Traceback = new TracebackType[3*(signalA_len+1)*(signalB_len+1)];
        std::memcpy(Traceback, other.Traceback, sizeof(TracebackType) * 3 * (signalA_len+1)*(signalB_len+1));
        return *this;
    }

    // Rule 3 Not a default destructor
    ~AffineAlignObj()
    {
        delete[] M; // since we declared with new, manually clear memory from heap.
        delete[] A;
        delete[] B;
        delete[] Traceback;
    }
};

struct AffineAlignObj1
{
  float* M;
  int signalA_len; // stack allocation
  int signalB_len;
  float GapOpen;
  float GapExten;
  bool FreeEndGaps;

  // Not a default constructor
  AffineAlignObj1(int ROW_SIZE, int COL_SIZE)
  {
    M = new float[ROW_SIZE * COL_SIZE]; // heap allocation
    signalA_len = ROW_SIZE-1;
    signalB_len = COL_SIZE-1;
  }

  // Rule 1 Copy constructor
  AffineAlignObj1(const AffineAlignObj1 &other)
  {
    signalA_len = other.signalA_len;
    signalB_len = other.signalB_len;
    GapOpen = other.GapOpen;
    GapExten = other.GapExten;
    FreeEndGaps = other.FreeEndGaps;
    M = new float[(signalA_len+1)*(signalB_len+1)];
    std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
  }

  // Rule 2 Copy assignment operator
  AffineAlignObj1& operator=(const AffineAlignObj1& other)
  {
    if(this == &other) return *this; // handling of self assignment.
    delete[] M; // freeing previously used memory
    signalA_len = other.signalA_len;
    signalB_len = other.signalB_len;
    GapOpen = other.GapOpen;
    GapExten = other.GapExten;
    FreeEndGaps = other.FreeEndGaps;
    M = new float[(signalA_len+1)*(signalB_len+1)];
    std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
    return *this;
  }

  // Rule 3 Not a default destructor
  ~AffineAlignObj1()
  {
    delete[] M; // since we declared with new, manually clear memory from heap.
  }
};

#endif // AFFINEALIGNOBJ_H
