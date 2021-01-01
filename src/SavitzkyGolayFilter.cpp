#include <RcppEigen.h>
#include "SavitzkyGolayFilter.h"

namespace DIAlign
{
void SavitzkyGolayFilter::updateMembers(){
  for (int nl = 0; nl <= (int) (frame_size_ / 2); ++nl)
  {
    coeffs_.resize(frame_size_ * (frame_size_ / 2 + 1));
    int nr = frame_size_ - 1 - nl;

    // compute a Vandermonde matrix whose columns are powers of the vector [-nL,...,nR]
    Eigen::MatrixXd A (frame_size_, order_ + 1);
    for (int i = -nl; i <= nr; i++)
    {
      for (int j = 0; j <= static_cast<int>(order_); j++)
      {
        A(i + nl, j) = std::pow((float)i, j); // pow(int, int) is not defined
      }
    }

    // compute the singular-value decomposition of A
    Eigen::JacobiSVD<Eigen::MatrixXd> svd (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd B (order_ + 1);
    for (int i = 0; i <= order_; ++i)
    {
      B(i) = svd.matrixV()(0, i) / svd.singularValues()(i);
    }

    // compute B*transpose(U)*b, where b is the unit vector b=[1 0 ... 0]
    for (int i = 0; i < frame_size_; ++i)
    {
      coeffs_[(nl + 1) * frame_size_ - i - 1] = 0;
      for (int j = 0; j <= order_; ++j)
      {
        coeffs_[(nl + 1) * frame_size_ - i - 1] += B(j) * svd.matrixU()(i, j);
      }
    }
  }
}

std::vector<double> smoothChroms(std::vector<double> intensity, SavitzkyGolayFilter sgolay){
  PeakIntegration::MSChromatogram chromatogram;
  int len = intensity.size();
  chromatogram.resize(len);
  PeakIntegration::MSChromatogram::Iterator it = chromatogram.begin();
  for (int i=0; i<len; ++i, ++it)
  {
    it->setIntensity(intensity[i]);
  }
  sgolay.filter(chromatogram);
  std::vector<double> newChrom(len);
  it = chromatogram.begin();
  for (int i=0; i<len; ++i, ++it)
  {
    newChrom[i] = it->getIntensity();
  }

  return newChrom;
}
}
