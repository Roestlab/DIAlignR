#include <algorithm>
#include <stdexcept>
#include "miscell.h"

namespace DIAlign
{
void xicIntersect(std::vector<std::vector<double> > &time,
                  std::vector<std::vector<double> > &intensity){
  int len = intensity.size();
  std::vector<std::vector<int> > t;
  for (unsigned int i = 0; i < len; i++){
    std::vector<int> k(time[i].begin(), time[i].end());
    t.push_back(k);
  }
  int strt = 0, end = 1e7;
  for(unsigned int i = 0; i < len; i++){
    if(t[i][0] > strt) {
      strt = t[i][0];
    }
    if(t[i].back() < end) {
      end = t[i].back();
    }
  }

  // Get sub chromatogram
  int strti, endi;
  for (unsigned int i = 0; i < len; i++){
    auto lower = std::find(t[i].begin(), t[i].end(), strt);
    strti = std::distance(t[i].begin(), lower);
    auto upper = std::find(t[i].begin(), t[i].end(), end);
    endi = std::distance(t[i].begin(), upper);
    intensity[i].erase(intensity[i].begin(),intensity[i].begin()+strti);
    intensity[i].erase(intensity[i].begin()+endi-strti+1,intensity[i].end());
    time[i].erase(time[i].begin(), time[i].begin()+strti);
    time[i].erase(time[i].begin()+endi-strti+1, time[i].end());
  }

  // Check if fragment-ions are of same length.
  std::vector<int> len2(len, 0);
  for (unsigned int i = 0; i < len; i++) len2[i] = intensity[i].size();
  if(std::adjacent_find(len2.begin(), len2.end(), std::not_equal_to<int>()) != len2.end()){
    throw std::length_error("Fragment-ion vectors must have same length");
  }
}

static bool const detect_end_na(double a, double b){
  return (a < 0) && !(b < 0);
};

static bool const detect_start_na(double a, double b){
  return !(a < 0) && (b < 0);
};


void interpolateZero(std::vector<double> & x){
  auto start = x.begin();
  auto end = x.end();
  if(x[0] < 0) start = std::adjacent_find(start, end, detect_start_na);
  while (true) {
    // Find transitions to and from NA values.  If we hit end of
    // vector whilst looking, our work is done.
    auto num_to_na = std::adjacent_find(start, end, detect_start_na);
    auto na_to_num = std::adjacent_find(start, end, detect_end_na);
    if (na_to_num == end) {
      break;
    }

    // At this point, num_to_na points to the last number before
    // an interpolation block, and na_to_num points to the last NA
    // of that block.

    ++na_to_num;            // Now, both iterators point to numbers.
    auto const base = *num_to_na;
    auto const target = *na_to_num;

    // To count rails rather than posts, measure difference before
    // incrementing the start position.
    auto const gaps = std::distance(num_to_na, na_to_num);

    ++num_to_na;
    // Now both iterators point immediately *after* transition.

    auto const make_value = [base, target, gaps, i = std::size_t{0}]()
      mutable { return base + (++i * (target - base) / gaps); };
    std::generate(num_to_na, na_to_num, make_value);

    // Advance onwards
    start = na_to_num;
  }
}
}
