// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 20022-200.
// <Shubham Gupta PeakIntegrator.cpp>
// Copyright (C) 2020-2040 Shubham Gupta
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Shubham Gupta$
// $Authors: Shubham Gupta, Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include "PeakIntegrator.h"

namespace PeakIntegration
{
const String String::EMPTY;

String::String() :
  std::string()
{
}

String::String(const char* s) :
  std::string(s)
{
}


void Param::setIntegrationType(std::string i){
  integration_type_ = i;
}

void Param::setBaselineType(std::string b){
  baseline_type_ = b;
}

const char* Param::getIntegrationType(){
  return(integration_type_.c_str());
}

const char* Param::getBaselineType(){
  return(baseline_type_.c_str());
}

PeakIntegrator::PeakIntegrator()
{
  Param params;
  getDefaultParameters(params);
  updateMembers(params); // write defaults into Param object param_
}

PeakIntegrator::~PeakIntegrator() {}

PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, const double left, const double right) const
{
  return integratePeak_(chromatogram, left, right);
}

PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right) const
{
  return integratePeak_(chromatogram, left->getRT(), right->getRT());
}

PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, const double left, const double right, const double peak_apex_pos) const
{
  return estimateBackground_(chromatogram, left, right, peak_apex_pos);
}

PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_apex_pos) const
{
  return estimateBackground_(chromatogram, left->getRT(), right->getRT(), peak_apex_pos);
}

PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, const double left, const double right, const double peak_height, const double peak_apex_pos) const
{
  return calculatePeakShapeMetrics_(chromatogram, left, right, peak_height, peak_apex_pos);
}

PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_height, const double peak_apex_pos) const
{
  return calculatePeakShapeMetrics_(chromatogram, left->getRT(), right->getRT(), peak_height, peak_apex_pos);
}

void PeakIntegrator::getDefaultParameters(Param& params)
{
  params.setIntegrationType((std::string)INTEGRATION_TYPE_INTENSITYSUM); // "The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration."
  params.setBaselineType((std::string)BASELINE_TYPE_BASETOBASE); // "The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base)."
}

void PeakIntegrator::updateMembers(Param& param)
{
  integration_type_ = (String)param.getIntegrationType();
  baseline_type_ = (String)param.getBaselineType();
  fit_EMG_ = false;
}

}
