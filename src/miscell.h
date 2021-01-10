#ifndef XICINTERSECTS_H
#define XICINTERSECTS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include "spline.h"

namespace DIAlign
{

static bool const detect_end_na(double a, double b);

static bool const detect_start_na(double a, double b);

static bool const lessZero(double a);

void xicIntersect(std::vector<std::vector<double> > & time, std::vector<std::vector<double> > & intensity);

void interpolateZero(std::vector<double> & x);

std::vector<int> getNegIndices(const std::vector<double> & A);

std::vector<std::vector<double>> imputeChromatogram(const std::vector<std::vector<double>> & A,
                                       const std::vector<double> & t, const std::vector<int> & index);

std::vector<int> getFlank(const std::vector<double> & t1, const std::vector<double> & t2);
std::vector<int> getSkip(const std::vector<int> & index, const std::vector<int> & flank);
std::vector<int> getFlankN(const std::vector<double> & t, const std::vector<int> & flank);
std::vector<int> getKeep(int length, std::vector<int> skip);
void addFlankToLeft(const std::vector<double> & t, std::vector<double> & tN, std::vector<double> & tA,
                    const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                    std::vector<int> & flank);
void addFlankToRight(const std::vector<double> & t, std::vector<double> & tN, std::vector<double> & tA,
                     const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                     std::vector<int> & flank);

void mergeTime(std::vector<double> & t1, const std::vector<double> & t2, std::string mergeStrategy);

void mergeIntensity(std::vector<std::vector<double>> & A,
                    std::vector<std::vector<double>> & B, double w);


} // namespace DIAlign
#endif // XICINTERSECTS_H
