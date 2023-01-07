#include <algorithm>
#include <numeric> // std::iota
#include <vector>

namespace std {
template <typename T> void argsort(const vector<T> &v, vector<size_t> &idx) {
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
}

template <typename T> void diff(const vector<T> &v, vector<T> &delta) {
  for (size_t i = 0; i < v.size() - 1; ++i) {
    delta.push_back(v[i + 1] - v[i]);
  }
}

template <typename T>
void concatenateSort(const vector<T> &A, const vector<T> &B, vector<T> &C) {
  C.reserve(A.size() + B.size());
  C.insert(C.end(), A.begin(), A.end());
  C.insert(C.end(), B.begin(), B.end());
  sort(C.begin(), C.end());
}

template <typename T>
void searchSorted(const vector<T> &v, const vector<size_t> &idx,
                  const vector<T> &allValues, vector<size_t> &cdfIdx) {
  // for every element in all indices, find id in A
  cdfIdx.reserve(allValues.size());
  for (auto elem = allValues.begin(); elem != allValues.end() - 1; ++elem) {
    auto up = upper_bound(idx.begin(), idx.end(), *elem,
                          [&v](T elem, size_t i) { return elem < v[i]; });
    if (up == idx.end()) {
      cdfIdx.push_back(idx.size());
    } else {
      cdfIdx.push_back(up - idx.begin());
    }
  }
}

template <typename T>
void computeCDF(const vector<T> &weights, vector<size_t> &idx,
                const vector<size_t> &cdfIdx, vector<T> &cdf) {

  vector<T> ordW;
  for (auto i : idx) {
    ordW.push_back(weights[i]);
  }

  cdf.reserve(1 + ordW.size());
  vector<T> sortedAccW = {0};
  partial_sum(ordW.begin(), ordW.end(), ordW.begin());
  sortedAccW.insert(sortedAccW.end(), ordW.begin(), ordW.end());
  T accum = sortedAccW.back();
  for (auto cIdx : cdfIdx) {
    cdf.push_back(sortedAccW[cIdx] / accum);
  }
}

template <typename T>
T computeDist(const vector<T> &cdfA, const vector<T> &cdfB,
              const vector<T> &deltas) {

  T dist = 0.0;
  for (size_t i = 0; i < deltas.size(); ++i) {
    dist += abs(cdfA[i] - cdfB[i]) * deltas[i];
  }
  return dist;
}

template <typename T>
T wasserstein(vector<T> &A, vector<T> AWeights, vector<T> &B,
              vector<T> BWeights) {

  vector<T> allValues;
  concatenateSort(A, B, allValues);

  vector<T> deltas;
  diff(allValues, deltas);

  vector<size_t> idxA(A.size());
  argsort(A, idxA);

  vector<size_t> idxB(B.size());
  argsort(B, idxB);

  vector<size_t> cdfIdxA;
  searchSorted(A, idxA, allValues, cdfIdxA);

  vector<size_t> cdfIdxB;
  searchSorted(B, idxB, allValues, cdfIdxB);

  vector<T> cdfA;
  computeCDF(AWeights, idxA, cdfIdxA, cdfA);

  vector<T> cdfB;
  computeCDF(BWeights, idxB, cdfIdxB, cdfB);

  return computeDist(cdfA, cdfB, deltas);
}
} // namespace std
