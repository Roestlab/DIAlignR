#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "alignment.h"

#include "chromSimMatrix.h"
#include "alignment.h"
#include "gapPenalty.h"
#include "affinealignobj.h"
#include "affinealignment.h"
#include "constrainMat.h"

using namespace DIAlign;

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

struct ChromData
{
  std::vector<double > rt_data;
  std::vector< std::vector<double > > data;
};

ChromData readFile(std::string inputfile)
{
  ChromData res;

  std::ifstream file(inputfile.c_str()); //Read File in tsv format, with header-line

  //read the header
  std::string meta;
  std::string header;
  getline(file, meta, '\n'); 
  getline(file, header, '\n'); 

  if (meta.empty() || meta[0] != '#')
  {
    std::cout << "First line needs to start with character #" << std::endl;
  }
  meta = meta.substr(1);

  int nr_transitions = std::atoi(meta.c_str());

  std::string line;
  res.data.resize(nr_transitions);
  while(getline(file, line, '\n'))
  {
    auto vec = split(line, '\t');
    res.rt_data.push_back( std::atoi(vec[0].c_str()) );
    for (int k = 1; k <= nr_transitions; k++)
    {
      res.data[k-1].push_back( std::atoi(vec[k].c_str()) );
    }
  }

  // for (auto k : data)
  // {
  //   for (auto el : k)
  //     std::cout << el << " ";
  //   std::cout << std::endl;
  // }
  return res;
}

void alignChromatogramsCpp( AffineAlignObj& obj,
                            const std::vector<std::vector<double> > & r1,
                            const std::vector<std::vector<double> > & r2,
                            std::string alignType,
                            const std::vector<double>& tA, const std::vector<double>& tB,
                            const std::string & normalization, const std::string& simType,
                            double B1p = 0.0, double B2p =0.0, int noBeef = 0,
                            double goFactor = 0.125, double geFactor = 40,
                            double cosAngleThresh = 0.3, bool OverlapAlignment = true,
                            double dotProdThresh = 0.96, double gapQuantile = 0.5,
                            bool hardConstrain = false, double samples4gradient = 100.0)
{

  SimMatrix s = getSimilarityMatrix(r1, r2, normalization, simType, cosAngleThresh, dotProdThresh);
  obj.reset(s.n_row + 1, s.n_col + 1);

  double gapPenalty = getGapPenalty(s, gapQuantile, simType);
  if (alignType == "hybrid")
  {
    SimMatrix MASK;
    MASK.n_row = tA.size();
    MASK.n_col = tB.size();
    MASK.data.resize(MASK.n_row*MASK.n_col, 0.0);
    double A1 = tA[0], A2 = tA[MASK.n_row-1];
    double B1 = tB[0], B2 = tB[MASK.n_col-1];
    calcNoBeefMask(MASK, A1, A2, B1, B2, B1p, B2p, noBeef, hardConstrain);
    auto maxIt = max_element(std::begin(s.data), std::end(s.data));
    double maxVal = *maxIt;
    constrainSimilarity(s, MASK, -2.0*maxVal/samples4gradient);
  }
  doAffineAlignment(obj, s, gapPenalty*goFactor, gapPenalty*geFactor, OverlapAlignment); // Performs alignment on s matrix and returns AffineAlignObj struct
  getAffineAlignedIndices(obj); // Performs traceback and fills aligned indices in AffineAlignObj struct
}

void doAlignment()
{
  auto data1 = readFile("test_large.tsv");
  auto data2 = readFile("test_large.tsv");

  const std::vector<double> tA;
  const std::vector<double> tB;
  int noBeef = 5; // nr of matrix cells w/o penalty
  AffineAlignObj obj(256, 256); // Initialize AffineAlignObj
  alignChromatogramsCpp(obj, data1.data, data2.data, "hybrid", 
      data1.rt_data, data2.rt_data, "mean", "dotProductMasked", data2.rt_data.front(), data2.rt_data.back(), noBeef);
}

int main(){
  doAlignment();
  return 0;
}
