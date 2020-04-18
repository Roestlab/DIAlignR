#ifndef INTEGRATEAREA_H
#define INTEGRATEAREA_H

#include <vector>
#include "utils.h"
#include "PeakIntegrator.h"

namespace DIAlign
{
namespace PeakGroupIntensity
{
   /**
    * @brief returns the summation of signals between leftIdx and rightIdx from vov.
    *
    */
   std::vector<std::vector<double> > peakGroupArea(std::vector<std::vector<double> > position, std::vector<std::vector<double> > intensity,
                        double left, double right, const std::string integrationType, const std::string baselineType, bool fitEMG);
} //namespace PeakGroupIntensity
} // namespace DIAlign

#endif // INTEGRATEAREA_H
