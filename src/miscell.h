#ifndef XICINTERSECTS_H
#define XICINTERSECTS_H

#include <vector>

namespace DIAlign
{

static bool const detect_end_na(double a, double b);

static bool const detect_start_na(double a, double b);

void xicIntersect(std::vector<std::vector<double> > & time, std::vector<std::vector<double> > & intensity);

void interpolateZero(std::vector<double> & x);
} // namespace DIAlign
#endif // XICINTERSECTS_H
