#include "integrateArea.h"
// using namespace PeakIntegration;

namespace DIAlign
{
namespace PeakGroupIntensity
{
    std::vector<std::vector<double> > peakGroupArea(std::vector<std::vector<double> > position, std::vector<std::vector<double> > intensity,
                        double left, double right, const std::string integrationType, const std::string baselineType, bool fitEMG, bool baseline_subtraction){
    std::vector<double> area(position.size(), 0.0);
    std::vector<double> apex(position.size(), 0.0);
    double peak_integral = 0.0;
    double peak_apex_int = 0.0;

    PeakIntegration::PeakIntegrator* ptr = new PeakIntegration::PeakIntegrator();
    PeakIntegration::Param params;
    params.setIntegrationType(integrationType);
    params.setBaselineType(baselineType);
    ptr->updateMembers(params);

    PeakIntegration::PeakIntegrator::PeakArea pa;
    PeakIntegration::PeakIntegrator::PeakBackground pb;

    int n_frag = position.size();
    for(int fragIon = 0; fragIon < n_frag; fragIon++){
      PeakIntegration::MSChromatogram chromatogram;
      for (PeakIntegration::Size i = 0; i < position[fragIon].size(); ++i)
      {
        chromatogram.push_back(PeakIntegration::ChromatogramPeak(position[fragIon][i], intensity[fragIon][i]));
      }
      pa = ptr->integratePeak(chromatogram, left, right);
      pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
      if(baseline_subtraction){
        peak_integral = pa.area - pb.area;
        peak_apex_int = pa.height - pb.height;
      } else {
        peak_integral = pa.area;
        peak_apex_int = pa.height;
      }
      if (peak_integral < 0) {peak_integral = 0;}
      if (peak_apex_int < 0) {peak_apex_int = 0;}
      area[fragIon] = peak_integral;
      apex[fragIon] = peak_apex_int;
    }
    std::vector<std::vector<double> > output;
    output.push_back(area);
    output.push_back(apex);
    return output;
  }
} //namespace PeakGroupIntensity
} // namespace DIAlign
