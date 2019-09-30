#include "integrateArea.h"
#include <numeric>

namespace DIAlign
{
double areaBwBoundaries(std::vector<std::vector<double> > vov, int leftIdx, int rightIdx){
  double area = 0.0;
  int n_frag = vov.size();
  for(int fragIon = 0; fragIon < n_frag; fragIon++){
    area += std::accumulate(vov[fragIon].begin()+leftIdx, vov[fragIon].begin()+rightIdx+1, 0);
  }
  return area;
}
}
