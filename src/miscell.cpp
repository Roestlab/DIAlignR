#include <stdexcept>
#include <cstdlib>
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

static bool const lessZero(double a){
  return (a < 0);
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

std::vector<int> getNegIndices(const std::vector<double> & A){
  std::vector<int> B;
  std::vector<double>::const_iterator it = A.begin();
  while ((it= std::find_if(it, A.end(), [](int x){return x < 0; })) != A.end())
  {
    B.push_back(std::distance(A.begin(), it));
    it++;
  }
  return B;
}

std::vector<std::vector<double>> imputeChromatogram(const std::vector<std::vector<double>> & A,
                                                    const std::vector<double> & t,
                                                    const std::vector<int> & index){
  int nrow = index.size();
  // Expand time to indices
  std::vector<double> tnew(nrow, -1.0);
  for(int i= 0; i<nrow; i++){
    if(index[i] != 0){
      tnew[i] = t[index[i]-1];
    }
  }
  std::vector<int> s1 = getNegIndices(tnew);

  // Fill missing values like zoo::na.approx
  interpolateZero(tnew);
  std::vector<int> s2 = getNegIndices(tnew);
  std::vector<int> middle;
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(middle, middle.end()));
  std::vector<double> xout(middle.size());
  for(int i =0; i< middle.size(); i++){
    xout[i] = tnew[middle[i]];
  }

  // Interpolate intensity for each fragment.
  std::vector<std::vector<double>> Anew(A.size()+1);
  for(int i =0; i < (Anew.size()-1); i++){
    std::vector<double> intensity(nrow, -1.0);
    for(int j= 0; j<nrow; j++){
      if(index[j] != 0){
        intensity[j] = A[i][index[j]-1];
      }
    }
    std::vector<double> result = naturalSpline(t, A[i], xout);
    for(int j=0; j<result.size(); j++){
      intensity[middle[j]] = result[j];
    }
    Anew[i] = intensity;
  }

  // Append time vector with fragments intensities.
  Anew.back() = tnew;
  return Anew;
}

std::vector<int> getFlank(const std::vector<double> & t1, const std::vector<double> & t2){
  std::vector<int> s1 = getNegIndices(t1);
  std::vector<int> s2 = getNegIndices(t2);
  std::vector<int> flank(s1.size()+s2.size());
  std::vector<int>::iterator it;
  it = std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), flank.begin());
  flank.resize(it-flank.begin());
  return flank;
}

std::vector<int> getSkip(const std::vector<int> & index, const std::vector<int> & flank){
  int noKeep = std::count(index.begin(), index.end(), 0);
  std::vector<int> s3(noKeep);
  int j =0;
  for(int i = 0; i < index.size(); i++){
    if(index[i] == 0){
      s3[j] = i;
      ++j;
    }
  }
  std::vector<int> skip(flank.size()+s3.size());
  std::vector<int>::iterator it;
  it = std::set_union(flank.begin(), flank.end(), s3.begin(), s3.end(), skip.begin());
  skip.resize(it-skip.begin());
  return skip;
}

std::vector<int> getFlankN(const std::vector<double> & t, const std::vector<int> & flank){
  std::vector<int> s = getNegIndices(t);
  std::vector<int> flankN;
  std::set_intersection(s.begin(), s.end(), flank.begin(), flank.end(),
                        std::inserter(flankN, flankN.end()));
  return flankN;
}

std::vector<int> getKeep(const int length, const std::vector<int> & skip){
  std::vector<int> ivec(length);
  std::iota(ivec.begin(), ivec.end(), 0);
  std::vector<int> keep(length-skip.size());
  std::vector<int>::iterator it = std::set_difference(ivec.begin(), ivec.end(),
                                  skip.begin(), skip.end(), keep.begin());
  keep.resize(it-keep.begin());
  return keep;
}

void mergeTime(std::vector<double> & t1, const std::vector<double> & t2,
                              std::string mergeStrategy){
  int np = t1.size();
  if(mergeStrategy == "ref"){
    return;
  } else if (mergeStrategy == "avg" || mergeStrategy == "refStart" || mergeStrategy == "refEnd") {
    std::transform(t1.begin(), t1.end(), t2.begin(), t1.begin(), std::plus<double>());
    std::transform(t1.begin(), t1.end(), t1.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, 0.5));
    double delta = t1[np-1] - t1[0]/(np-1);
    double val;
    std::vector<double>::iterator it;
    if (mergeStrategy == "refStart"){
      for (it = t1.begin(), val = t1[0]; it != t1.end(); ++it, val += delta) {
        *it = val;  // changes the values in the vector
      }
    } else if (mergeStrategy == "refEnd"){
      for (it = t1.begin(), val = t1[np-1]-delta*(np-1); it != t1.end(); ++it, val += delta) {
        *it = val;  // changes the values in the vector
      }
    }
  }
}

void mergeIntensity(std::vector<std::vector<double>> & A,
                    std::vector<std::vector<double>> & B, double w){
  for(int i = 0; i< A.size(); i++){
    std::transform(A[i].begin(), A[i].end(), A[i].begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, w));
    std::transform(B[i].begin(), B[i].end(), B[i].begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, (1.0-w)));
    std::transform(A[i].begin(), A[i].end(), B[i].begin(), A[i].begin(), std::plus<double>());
  }
  return;
}

void addFlankToLeft(const std::vector<double> & t, std::vector<double> & tN, std::vector<double> & tA,
                    const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                    std::vector<int> & flank){
  auto it = std::adjacent_find(flank.begin(), flank.end(), [](int l, int r){return l+1<r;});
  int length = (it == flank.end()) ? it-flank.begin() : it-flank.begin()+1;
  std::vector<double> t_flank(length);
  double t_new = tN[0] - t[length];
  std::copy(t.begin(), t.begin()+length, t_flank.begin());
  std::transform(t_flank.begin(), t_flank.end(), t_flank.begin(),
                 std::bind2nd(std::plus<double>(), t_new));
  tN.insert(tN.begin(), t_flank.begin(), t_flank.end());
  std::copy(t_flank.begin(), t_flank.end(), tA.begin());
  for(int i = 0; i< intenN.size(); i++){
    intenN[i].insert(intenN[i].begin(), inten[i].begin(), inten[i].begin()+ length);
  }
  flank.erase(flank.begin(), flank.begin()+length);
  return;
}

void addFlankToRight(const std::vector<double> & t, std::vector<double> & tN, std::vector<double> & tA,
                     const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                    std::vector<int> & flank){
  int flankStart, eraseIdx;
  if(flank[0] != 0){
    eraseIdx = 0;
    flankStart = flank[0];
  } else{
    auto it = std::adjacent_find(flank.begin(), flank.end(), [](int l, int r){return l+1<r;});
    eraseIdx = std::distance(flank.begin(), it+1);
    flankStart = *(it+1);
  }
  int np = flank.back();
  double t_new = tN.back() - t[flankStart-1];
  std::vector<double> t_flank(np - flankStart +1);
  std::copy(t.begin()+flankStart, t.begin()+np+1, t_flank.begin());
  std::transform(t_flank.begin(), t_flank.end(), t_flank.begin(),
                 std::bind2nd(std::plus<double>(), t_new));
  tN.insert(tN.end(), t_flank.begin(), t_flank.end());
  std::copy(t_flank.begin(), t_flank.end(), tA.begin()+flankStart);
  for(int i = 0; i< intenN.size(); i++){
    intenN[i].insert(intenN[i].end(), inten[i].begin()+flankStart, inten[i].begin()+np+1);
  }
  flank.erase(flank.begin()+eraseIdx, flank.end());
  return;
}

void addFlankToLeft1(const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                    std::vector<int> & flank){
  auto it = std::adjacent_find(flank.begin(), flank.end(), [](int l, int r){return l+1<r;});
  int length = (it == flank.end()) ? it-flank.begin() : it-flank.begin()+1;
  for(int i = 0; i< intenN.size(); i++){
    intenN[i].insert(intenN[i].begin(), inten[i].begin(), inten[i].begin()+ length);
  }
  flank.erase(flank.begin(), flank.begin()+length);
  return;
}

void addFlankToRight1(const std::vector<std::vector<double>> & inten, std::vector<std::vector<double>> & intenN,
                     std::vector<int> & flank){
  int flankStart, eraseIdx;;
  if(flank[0] != 0){
    eraseIdx = 0;
    flankStart = flank[0];
  } else{
    auto it = std::adjacent_find(flank.begin(), flank.end(), [](int l, int r){return l+1<r;});
    eraseIdx = std::distance(flank.begin(), it+1);
    flankStart = *(it+1);
  }
  int np = flank.back();
  for(int i = 0; i< intenN.size(); i++){
    intenN[i].insert(intenN[i].end(), inten[i].begin()+flankStart, inten[i].begin()+np+1);
  }
  flank.erase(flank.begin()+eraseIdx, flank.end());
  return;
}

std::vector<int> getMatchingIdx(const std::vector<double> & tMain,
                                const std::vector<double> & t){
  // Collect matching indices of t in tMain
  std::vector<int> tIndex(tMain.size(), -1);
  for(int i =0, j=0; i< tMain.size(); i++){
    double tM = tMain[i];
    for(; j<t.size();){
      if(std::abs(tM - t[j])< 0.01){
        tIndex[i] = j;
        j++;
      }
      if(tM - t[j] < 0.0) break;
      ++j;
    }
  }
  return tIndex;
}

std::vector<std::vector<double>> imputeChromatogram1(const std::vector<std::vector<double>> & A,
                                                     const std::vector<int> & tIndex, const std::vector<double> & t,
                                                     const std::vector<double> & tnew){
  // Fill intensity for which there is a match as per tIndex.
  std::vector<std::vector<double>> intensity(A.size());
  for(int j = 0; j <A.size(); j++){
    std::vector<double> temp(tnew.size(), -1);
    for(int i = 0; i < temp.size(); i++){
      if(tIndex[i] != -1) temp[i] = A[j][tIndex[i]];
    }
    intensity[j] = temp;
  }

  // Get time for which no match as per tIndex.
  std::vector<double> doubleVec(tIndex.begin(), tIndex.end());
  std::vector<int> s1 = getNegIndices(doubleVec);
  std::vector<int> s2 = getNegIndices(tnew);
  std::vector<int> middle;
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(middle, middle.end()));
  std::vector<double> xout(middle.size());
  for(int i =0; i< middle.size(); i++){
    xout[i] = tnew[middle[i]];
  }

  // Interpolate intensity for each fragment.
  for(int i =0; i < intensity.size(); i++){
    std::vector<double> result = naturalSpline(t, A[i], xout);
    for(int j=0; j<result.size(); j++){
      intensity[i][middle[j]] = result[j];
    }
  }

  return intensity;
}

}
