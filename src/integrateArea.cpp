#include "integrateArea.h"

namespace DIAlign
{
namespace PeakIntegration
{
  double peakGroupArea(std::vector<std::vector<double> > position, std::vector<std::vector<double> > intensity,
                        double left, double right, const std::string integrationType, const std::string baselineType, bool fitEMG){
    double area = 0.0;

    /***
    PeakIntegrator* ptr = new PeakIntegrator();
    Param params;
    params.setIntegrationType(integrationType);
    params.setBaselineType(baselineType);
    ptr->updateMembers(params);

    PeakIntegrator::PeakArea pa;
    PeakIntegrator::PeakBackground pb;

    int n_frag = position.size();
    for(int fragIon = 0; fragIon < n_frag; fragIon++){

      MSChromatogram chromatogram;
      for (Size i = 0; i < position[fragIon].size(); ++i)
      {
        chromatogram.push_back(ChromatogramPeak(position[fragIon][i], intensity[fragIon][i]));
      }
      pa = ptr->integratePeak(chromatogram, left, right);
      pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
      area += pb.area;
    }*/
    return area;
  }
} //namespace PeakIntegration
} // namespace DIAlign
