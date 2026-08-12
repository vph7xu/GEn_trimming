#ifndef HCAL_CLUSTER_SELECTOR_H
#define HCAL_CLUSTER_SELECTOR_H

#include <algorithm>
#include <cstddef>

namespace HCalClusterSelector {

struct BestResult {
  int  best_idx = 0;   // index in the re-ordered arrays
  bool pass     = true;
};

inline void swapd(double &a, double &b){
  double t = a; a = b; b = t;
}

// swap if both arrays exist
inline void swap_if(double* a, double* b, int i, int j){
  if(a && b) swapd(a[i], b[j]);
}
inline void swap_if(double* a, int i, int j){
  if(a) swapd(a[i], a[j]);
}

/**
 * Sort clusters by descending energy, swapping all provided arrays in lockstep.
 *
 * Required:
 *   E, atime, x, y must be valid pointers if nclus > 0.
 *
 * Optional:
 *   any of (tdctime, blkid, nblk, row, col, etc.) can be passed as nullptr.
 *
 * Notes:
 * - This is an in-place selection-sort (n is usually small, fine for HCAL clusters).
 * - Only indices [0, min(nclus, maxclus)) are considered.
 */
inline void SortByEnergy(int nclus,
                         double *E,
                         double *atime,
                         double *x,
                         double *y,
                         int maxclus = 1000,
                         // optional extras (swap in lockstep if provided)
                         double *tdctime = nullptr,
                         double *blkid   = nullptr,
                         double *nblk    = nullptr,
                         double *row     = nullptr,
                         double *col     = nullptr)
{
  const int n = std::max(0, std::min(nclus, maxclus));
  if(n <= 1) return;

  for(int i = 0; i < n - 1; i++){
    int maxIdx = i;
    for(int j = i + 1; j < n; j++){
      if(E[j] > E[maxIdx]) maxIdx = j;
    }
    if(maxIdx != i){
      // required arrays
      swapd(E[i],     E[maxIdx]);
      swapd(atime[i], atime[maxIdx]);
      swapd(x[i],     x[maxIdx]);
      swapd(y[i],     y[maxIdx]);

      // optional arrays
      if(tdctime) swapd(tdctime[i], tdctime[maxIdx]);
      if(blkid)   swapd(blkid[i],   blkid[maxIdx]);
      if(nblk)    swapd(nblk[i],    nblk[maxIdx]);
      if(row)     swapd(row[i],     row[maxIdx]);
      if(col)     swapd(col[i],     col[maxIdx]);
    }
  }
}

/**
 * Select best cluster among top Ntry (default 5) clusters (after SortByEnergy),
 * requiring: low < (clus_atime - sh_atime) < high where
 *   low  = mean - nsigma*sigma
 *   high = mean + nsigma*sigma
 */
inline BestResult SelectBestByCoincidence(int nclus,
                                         const double *atime,
                                         double sh_atime,
                                         double coin_mean,
                                         double coin_sigma,
                                         double nsigma,
                                         int ntry = 5)
{
  BestResult r;
  r.best_idx = 0;
  r.pass = true;

  const int n  = std::max(0, nclus);
  const int nt = std::min(ntry, n);

  const double low  = coin_mean - nsigma * coin_sigma;
  const double high = coin_mean + nsigma * coin_sigma;

  for(int i = 0; i < nt; i++){
    const double coint = atime[i] - sh_atime;
    if(low < coint && coint < high){
      r.best_idx = i;
      r.pass = true;
      return r;
    }
  }

  r.best_idx = 0;
  r.pass = false;
  return r;
}

} // namespace HCalClusterSelector

#endif

