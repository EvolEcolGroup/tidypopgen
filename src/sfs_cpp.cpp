/******************************************************************************/

#include "fbm_byte_acc.h"
#include <cmath>
#include <map>
#include <numeric>
#include <unordered_map>

/******************************************************************************/
// Site frequency spectra (SFS), with missing genotypes handled as in dadi
// (Gutenkunst et al. 2009) and easySFS: each population is projected down to
// a fixed number of allele copies with the hypergeometric distribution, and
// loci with fewer called allele copies than the projection in any population
// of a spectrum are dropped from that spectrum.
//
// - sfs_config_counts_cpp() (one pass over the FBM): for each requested
//   spectrum (a set of groups), the number of loci with each "configuration",
//   i.e. the number of called allele copies n and of alternate allele copies
//   k in each group of the spectrum (dadi's Spectrum.count_data_dict()).
// - sfs_project_cpp(): project the configurations of each spectrum into the
//   SFS (dadi's Spectrum._from_count_dict()).
// - sfs_preview_cpp(): the number of segregating sites in the projected 1D
//   SFS, for every projection (easySFS --preview).
/******************************************************************************/

// index of the pair (n, k), with 0 <= k <= n, in a triangular array
inline uint64_t sfs_tri(uint64_t n, uint64_t k) { return n * (n + 1) / 2 + k; }

// inverse of sfs_tri()
inline void sfs_untri(uint64_t t, int& n, int& k) {
  uint64_t nn = (uint64_t) ((std::sqrt(8.0 * (double) t + 1.0) - 1.0) / 2.0);
  while (nn > 0 && nn * (nn + 1) / 2 > t) nn--;
  while ((nn + 1) * (nn + 2) / 2 <= t) nn++;
  n = (int) nn;
  k = (int) (t - nn * (nn + 1) / 2);
}

/******************************************************************************/

// A counter of 64-bit keys in a flat hash table (open addressing with linear
// probing), which is much faster than std::unordered_map for the many
// increments of sfs_config_counts_cpp() (no allocation per key, and usually a
// single memory access per increment). UINT64_MAX marks empty slots, so it
// cannot be used as a key.
const uint64_t SFS_EMPTY_KEY = std::numeric_limits<uint64_t>::max();

class SfsFlatCounter {
public:
  SfsFlatCounter() : _size(0), _bits(10) {
    _keys.assign((size_t) 1 << _bits, SFS_EMPTY_KEY);
    _counts.assign((size_t) 1 << _bits, 0.0);
  }

  size_t size() const { return _size; }

  inline void increment(uint64_t key) {
    size_t mask = _keys.size() - 1;
    size_t i = slot(key);
    while (true) {
      if (_keys[i] == key) {
        _counts[i] += 1;
        return;
      }
      if (_keys[i] == SFS_EMPTY_KEY) {
        _keys[i] = key;
        _counts[i] = 1;
        _size++;
        // keep the load factor below 1/2
        if (2 * _size > _keys.size()) grow();
        return;
      }
      i = (i + 1) & mask;
    }
  }

  // the keys and their counts, sorted by key
  void sorted(std::vector<uint64_t>& keys, std::vector<double>& counts) const {
    std::vector<size_t> idx;
    idx.reserve(_size);
    for (size_t i = 0; i < _keys.size(); i++) {
      if (_keys[i] != SFS_EMPTY_KEY) idx.push_back(i);
    }
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      return _keys[a] < _keys[b];
    });
    keys.resize(idx.size());
    counts.resize(idx.size());
    for (size_t r = 0; r < idx.size(); r++) {
      keys[r] = _keys[idx[r]];
      counts[r] = _counts[idx[r]];
    }
  }

private:
  std::vector<uint64_t> _keys;
  std::vector<double> _counts;
  size_t _size;
  int _bits;

  // Fibonacci hashing of the key into the table
  inline size_t slot(uint64_t key) const {
    return (size_t) ((key * UINT64_C(0x9E3779B97F4A7C15)) >> (64 - _bits));
  }

  void grow() {
    std::vector<uint64_t> old_keys;
    std::vector<double> old_counts;
    old_keys.swap(_keys);
    old_counts.swap(_counts);
    _bits++;
    _keys.assign((size_t) 1 << _bits, SFS_EMPTY_KEY);
    _counts.assign((size_t) 1 << _bits, 0.0);
    size_t mask = _keys.size() - 1;
    for (size_t j = 0; j < old_keys.size(); j++) {
      if (old_keys[j] == SFS_EMPTY_KEY) continue;
      size_t i = slot(old_keys[j]);
      while (_keys[i] != SFS_EMPTY_KEY) i = (i + 1) & mask;
      _keys[i] = old_keys[j];
      _counts[i] = old_counts[j];
    }
  }
};

/******************************************************************************/

// Count the loci with each configuration of a spectrum. Each population i of
// the spectrum (group pops[i]) has at most max_n[i] allele copies; loci with
// fewer than min_n[i] called copies in any population are only counted as
// dropped. The configuration of a population is stored as the triangular
// index of (n, k), and the configurations of the populations are combined
// into a single key:
// - DENSE: a single population with few possible configurations, counted in
//   a dense vector indexed by the triangular index;
// - MAP64: the mixed radix number of the triangular indices of the
//   populations (the first population is the least significant), if it fits
//   in 64 bits (and is not UINT64_MAX), counted in an SfsFlatCounter;
// - MAPVEC: otherwise, a vector (n_1, k_1, ..., n_K, k_K) in an ordered map.
// Counts are integers (stored as doubles, exact up to 2^53), so they do not
// depend on the order in which loci are added.
class SfsConfigCounter {
public:
  SfsConfigCounter(const std::vector<int>& pops,
                   const std::vector<int>& max_n,
                   const std::vector<int>& min_n)
    : _pops(pops), _K(pops.size()), _min_n(min_n), _n_dropped(0),
      _n_unpolarised(0) {
    _radix.resize(_K);
    _n_tri.resize(_K);
    _mode = MAP64;
    uint64_t prod = 1;
    for (size_t i = 0; i < _K; i++) {
      // number of pairs (n, k) with 0 <= k <= n <= max_n
      _n_tri[i] = sfs_tri(max_n[i] + 1, 0);
      _radix[i] = prod;
      if (prod > (std::numeric_limits<uint64_t>::max() - 1) / _n_tri[i]) {
        _mode = MAPVEC;
      } else {
        prod *= _n_tri[i];
      }
    }
    if (_K == 1 && prod <= ((uint64_t) 1 << 21)) {
      _mode = DENSE;
      _dense.assign(prod, 0.0);
    }
  }

  // add the loci of a block; b_n and b_k hold the number of called and
  // derived allele copies of each locus (row) and group (column), with
  // ngroups values per locus; b_drop (if not NULL) flags the loci that cannot
  // be polarised, which are only counted as unpolarised
  void add_block(const int* b_n, const int* b_k, const char* b_drop,
                 size_t ngroups, size_t n_loci) {
    std::vector<int> key_vec(_mode == MAPVEC ? 2 * _K : 0);
    for (size_t j = 0; j < n_loci; j++) {
      if (b_drop != NULL && b_drop[j]) {
        _n_unpolarised += 1;
        continue;
      }
      const int* n = b_n + j * ngroups;
      const int* k = b_k + j * ngroups;
      bool keep = true;
      for (size_t i = 0; i < _K; i++) {
        if (n[_pops[i]] < _min_n[i]) {
          keep = false;
          break;
        }
      }
      if (!keep) {
        _n_dropped += 1;
        continue;
      }
      if (_mode == DENSE) {
        _dense[sfs_tri(n[_pops[0]], k[_pops[0]])] += 1;
      } else if (_mode == MAP64) {
        uint64_t key = 0;
        for (size_t i = 0; i < _K; i++) {
          key += sfs_tri(n[_pops[i]], k[_pops[i]]) * _radix[i];
        }
        _map64.increment(key);
      } else {
        for (size_t i = 0; i < _K; i++) {
          key_vec[2 * i] = n[_pops[i]];
          key_vec[2 * i + 1] = k[_pops[i]];
        }
        _mapvec[key_vec] += 1;
      }
    }
  }

  // the configurations (sorted by key) as a list with matrices n and k (one
  // row per configuration, one column per population), the number of loci
  // with each configuration (count), and the number of loci dropped for
  // missing data and because they could not be polarised
  List result() const {
    std::vector<int> n_out, k_out;
    std::vector<double> count_out;
    if (_mode == DENSE) {
      for (size_t t = 0; t < _dense.size(); t++) {
        if (_dense[t] > 0) {
          int n, k;
          sfs_untri(t, n, k);
          n_out.push_back(n);
          k_out.push_back(k);
          count_out.push_back(_dense[t]);
        }
      }
    } else if (_mode == MAP64) {
      std::vector<uint64_t> keys;
      std::vector<double> counts;
      _map64.sorted(keys, counts);
      for (size_t r = 0; r < keys.size(); r++) {
        for (size_t i = 0; i < _K; i++) {
          int n, k;
          sfs_untri((keys[r] / _radix[i]) % _n_tri[i], n, k);
          n_out.push_back(n);
          k_out.push_back(k);
        }
        count_out.push_back(counts[r]);
      }
    } else {
      for (const auto& kv : _mapvec) {
        for (size_t i = 0; i < _K; i++) {
          n_out.push_back(kv.first[2 * i]);
          k_out.push_back(kv.first[2 * i + 1]);
        }
        count_out.push_back(kv.second);
      }
    }
    size_t n_configs = count_out.size();
    IntegerMatrix n_mat(n_configs, _K), k_mat(n_configs, _K);
    NumericVector count(n_configs);
    for (size_t r = 0; r < n_configs; r++) {
      for (size_t i = 0; i < _K; i++) {
        n_mat(r, i) = n_out[r * _K + i];
        k_mat(r, i) = k_out[r * _K + i];
      }
      count[r] = count_out[r];
    }
    return List::create(_["n"] = n_mat, _["k"] = k_mat, _["count"] = count,
                        _["n_dropped"] = _n_dropped,
                        _["n_unpolarised"] = _n_unpolarised);
  }

private:
  enum Mode { DENSE, MAP64, MAPVEC };
  std::vector<int> _pops;
  size_t _K;
  std::vector<int> _min_n;
  std::vector<uint64_t> _radix, _n_tri;
  Mode _mode;
  std::vector<double> _dense;
  SfsFlatCounter _map64;
  std::map<std::vector<int>, double> _mapvec;
  double _n_dropped;
  double _n_unpolarised;
};

/******************************************************************************/

// Count the configurations of each spectrum in a single pass over the FBM.
//
// The number of called (n) and alternate (k) allele copies of each locus and
// group are computed with alt_valid_counts(), as in
// tajimas_d_window_sums_cpp(): every individual contributes as many allele
// copies as its ploidy, and pseudohaploids half their 0/2 dosage. Each
// element of spectra is an integer vector with the (0-based) groups of a
// spectrum; min_n gives, for each group, the smallest number of called
// allele copies for a locus to be counted (i.e. the projection). polarise is
// either empty (the reference allele is ancestral, and k counts alternate
// copies), or has a value for each locus: 0 if the reference allele is
// ancestral, 1 if the alternate allele is ancestral (k then counts reference
// copies), and -1 if the locus cannot be polarised (it is then dropped from
// all spectra). After each block of loci, the spectra are updated in
// parallel (each spectrum by a single thread). Returns a list with the
// result of SfsConfigCounter::result() for each spectrum.

// [[Rcpp::export]]
List sfs_config_counts_cpp(Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd,
                           const IntegerVector& groupIds,
                           int ngroups,
                           const NumericVector& ploidy,
                           bool is_pseudohaploid,
                           const List& spectra,
                           const IntegerVector& min_n,
                           const IntegerVector& polarise,
                           int ncores,
                           int block_size) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);
  if (block_size < 1) stop("block_size should be a positive integer");
  if (min_n.size() != ngroups) {
    stop("min_n should have one value per group");
  }
  bool has_polarise = polarise.size() > 0;
  if (has_polarise && (size_t) polarise.size() != m) {
    stop("polarise should be empty or have one value per locus");
  }
  for (R_xlen_t j = 0; j < polarise.size(); j++) {
    if (polarise[j] < -1 || polarise[j] > 1) {
      stop("polarise should only contain -1, 0 and 1");
    }
  }
  const int* p_polarise = polarise.begin();
  // allele counts need integer genotypes (not dosages)
  if (!acc.has_int_codes(255)) {
    stop("the site frequency spectrum needs integer genotype calls; it "
         "cannot be computed from dosages");
  }

  // pseudohaploid individuals contribute half their (0/2) dosage, and every
  // individual contributes as many allele copies as its ploidy
  std::vector<double> weight(n), copies(n);
  std::vector<int> group_max(ngroups, 0);
  for (size_t i = 0; i < n; i++) {
    weight[i] = (is_pseudohaploid && ploidy[i] == 1) ? 0.5 : 1;
    copies[i] = ploidy[i];
    group_max[group_ids[i]] += (int) ploidy[i];
  }

  // the counters of the spectra
  size_t n_spectra = spectra.size();
  std::vector<SfsConfigCounter> counters;
  counters.reserve(n_spectra);
  for (size_t s = 0; s < n_spectra; s++) {
    IntegerVector this_spectrum = spectra[s];
    if (this_spectrum.size() < 1) stop("each spectrum needs at least a group");
    std::vector<int> pops(this_spectrum.size()), max_n(pops.size()),
      this_min_n(pops.size());
    for (size_t i = 0; i < pops.size(); i++) {
      int g = this_spectrum[i];
      if (g < 0 || g >= ngroups) {
        stop("spectra should contain groups between 0 and ngroups - 1");
      }
      pops[i] = g;
      max_n[i] = group_max[g];
      this_min_n[i] = min_n[g];
    }
    counters.push_back(SfsConfigCounter(pops, max_n, this_min_n));
  }

  // called and derived allele copies of the loci of the current block, for
  // each group (element [(j - j0) * ngroups + g]; blocks start at multiples
  // of bs, so j - j0 = j % bs), whether a locus had a non-integer count, and
  // whether it cannot be polarised
  size_t bs = std::min((size_t) block_size, std::max(m, (size_t) 1));
  std::vector<int> b_n(bs * ngroups), b_k(bs * ngroups);
  std::vector<char> b_bad(bs, 0), b_drop(bs, 0);

  alt_valid_counts(acc, group_ids, ngroups, weight, copies, ncores, bs,
                   [&](size_t j, size_t g, double alt, double valid) {
                     size_t k = (j % bs) * ngroups + g;
                     b_n[k] = (int) valid;
                     b_k[k] = (int) alt;
                     // e.g. a heterozygous genotype in a pseudohaploid
                     if (alt != (double) b_k[k] || valid != (double) b_n[k]) {
                       b_bad[j % bs] = 1;
                     }
                     if (has_polarise) {
                       // the derived copies are the reference copies
                       if (p_polarise[j] == 1) b_k[k] = b_n[k] - b_k[k];
                       b_drop[j % bs] = (p_polarise[j] == -1);
                     }
                   },
                   [&](size_t j0, size_t j1) {
                     size_t n_block = j1 - j0;
                     for (size_t j = 0; j < n_block; j++) {
                       if (b_bad[j]) {
                         stop("non-integer allele counts (e.g. heterozygous "
                              "genotypes in pseudohaploid individuals)");
                       }
                     }
#pragma omp parallel for num_threads(ncores) schedule(dynamic, 1)
                     for (size_t s = 0; s < n_spectra; s++) {
                       counters[s].add_block(b_n.data(), b_k.data(),
                                             has_polarise ? b_drop.data() :
                                               NULL,
                                             (size_t) ngroups, n_block);
                     }
                   });

  List res(n_spectra);
  for (size_t s = 0; s < n_spectra; s++) {
    res[s] = counters[s].result();
  }
  return res;
}

/******************************************************************************/

// Projection of the configurations of each spectrum (dadi's
// Spectrum._from_count_dict()).
//
// A locus with n called copies in a population, k of which are alternate,
// contributes to cell j of the spectrum of that population the
// hypergeometric probability of drawing j alternate copies when sampling
// n_proj of its n copies without replacement (dhyper(j, k, n - k, n_proj)),
// which is non-zero for j in [max(0, n_proj - (n - k)), min(k, n_proj)]. Its
// contribution to a joint spectrum is the outer product of these vectors over
// the populations. Configurations with fewer than n_proj called copies in
// any population are dropped.
//
// Rather than adding the outer product of each configuration (which costs
// prod(n_proj + 1) operations per configuration), the configurations are
// sorted by the configuration of the last population, then the one before,
// and so on, and the outer products are factorised (a sparse tensor
// contraction): the vectors of the first population are summed over the
// configurations that share the configurations of all other populations,
// the result is multiplied by the vector of the second population and summed
// over the configurations that share the configurations of the third to last
// populations, and so on. Only the bounding box of the non-zero cells of each
// partial sum is visited. The result is the same (up to rounding), and the
// order of the operations only depends on the configurations.

// the hypergeometric probabilities of a configuration (n, k) of a population
struct SfsProjVec {
  int lo;        // first cell with a non-zero probability
  int hi;        // last cell with a non-zero probability
  size_t offset; // position of the probabilities in the flat storage
};

// the data needed to project a spectrum, as plain C++ objects
struct SfsSpectrumData {
  size_t K;                             // number of populations
  std::vector<int> p;                   // projection of each population
  std::vector<size_t> stride;           // stride of each dimension
  std::vector<uint64_t> keys;           // config of each kept row and pop
  std::vector<const SfsProjVec*> vecs;  // probabilities, idem
  std::vector<double> count;            // number of loci of each kept row
  double n_used;
  double n_dropped;
  std::vector<double> out;              // the spectrum (column-major)
};

// call f(i) for every linear index i of the box lo..hi (dimensions 0 to d)
// of an array with the given strides (dimension 0 is contiguous)
template <class F>
inline void sfs_box_for_each(int d, const int* lo, const int* hi,
                             const size_t* stride, F f) {
  std::vector<int> j(lo, lo + d + 1);
  while (true) {
    size_t base = 0;
    for (int i = 1; i <= d; i++) base += (size_t) j[i] * stride[i];
    for (int j0 = lo[0]; j0 <= hi[0]; j0++) f(base + j0);
    int i = 1;
    while (i <= d) {
      if (++j[i] <= hi[i]) break;
      j[i] = lo[i];
      i++;
    }
    if (i > d) break;
  }
}

// project the configurations of a spectrum (see above)
inline void sfs_contract(SfsSpectrumData& sd, const double* probs) {
  const size_t K = sd.K;
  const size_t n_rows = sd.count.size();
  // partial sums: T[d] has the dimensions 0 to d of the spectrum, and its
  // non-zero cells are within the box blo[d], bhi[d] (empty if
  // !has_box[d])
  std::vector<std::vector<double>> T(K);
  std::vector<std::vector<int>> blo(K), bhi(K);
  std::vector<char> has_box(K, 0);
  for (size_t d = 0; d < K; d++) {
    T[d].assign(sd.stride[d] * (sd.p[d] + 1), 0.0);
    blo[d].assign(d + 1, 0);
    bhi[d].assign(d + 1, 0);
  }
  // rows sorted by the configuration of the last population, then the one
  // before, ...
  std::vector<size_t> order(n_rows);
  std::iota(order.begin(), order.end(), 0);
  const uint64_t* keys = sd.keys.data();
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    for (size_t i = K; i-- > 0;) {
      if (keys[a * K + i] != keys[b * K + i]) {
        return keys[a * K + i] < keys[b * K + i];
      }
    }
    return a < b;
  });

  // multiply T[d] by the vector of population d + 1 of row r, add it to
  // T[d + 1], and reset T[d]
  auto flush = [&](size_t d, size_t r) {
    if (!has_box[d]) return;
    const SfsProjVec* v = sd.vecs[r * K + d + 1];
    const double* pr = probs + v->offset;
    double* t_in = T[d].data();
    double* t_out = T[d + 1].data();
    size_t st = sd.stride[d + 1];
    sfs_box_for_each((int) d, blo[d].data(), bhi[d].data(),
                     sd.stride.data(), [&](size_t idx) {
                       double val = t_in[idx];
                       if (val != 0) {
                         for (int j = v->lo; j <= v->hi; j++) {
                           t_out[idx + (size_t) j * st] += val * pr[j - v->lo];
                         }
                         t_in[idx] = 0;
                       }
                     });
    // the box of T[d + 1] now includes the box of T[d] times lo..hi
    for (size_t i = 0; i <= d + 1; i++) {
      int lo = (i <= d) ? blo[d][i] : v->lo;
      int hi = (i <= d) ? bhi[d][i] : v->hi;
      if (!has_box[d + 1]) {
        blo[d + 1][i] = lo;
        bhi[d + 1][i] = hi;
      } else {
        blo[d + 1][i] = std::min(blo[d + 1][i], lo);
        bhi[d + 1][i] = std::max(bhi[d + 1][i], hi);
      }
    }
    has_box[d + 1] = 1;
    has_box[d] = 0;
  };

  for (size_t o = 0; o < n_rows; o++) {
    size_t r = order[o];
    if (o > 0) {
      size_t prev = order[o - 1];
      // the last population whose configuration differs from the previous
      // row: the partial sums of the populations before it are complete
      size_t L = 0;
      for (size_t i = K; i-- > 0;) {
        if (keys[r * K + i] != keys[prev * K + i]) {
          L = i;
          break;
        }
      }
      for (size_t d = 0; d < L; d++) flush(d, prev);
    }
    // add the vector of the first population
    const SfsProjVec* v = sd.vecs[r * K];
    const double* pr = probs + v->offset;
    double w = sd.count[r];
    double* t0 = T[0].data();
    for (int j = v->lo; j <= v->hi; j++) t0[j] += w * pr[j - v->lo];
    if (!has_box[0]) {
      blo[0][0] = v->lo;
      bhi[0][0] = v->hi;
      has_box[0] = 1;
    } else {
      blo[0][0] = std::min(blo[0][0], v->lo);
      bhi[0][0] = std::max(bhi[0][0], v->hi);
    }
  }
  if (n_rows > 0) {
    for (size_t d = 0; d + 1 < K; d++) flush(d, order[n_rows - 1]);
  }
  sd.out.swap(T[K - 1]);
}

// Project the configurations of each spectrum. configs is a list with, for
// each spectrum, the output of sfs_config_counts_cpp() (n and k: one row per
// configuration and one column per population of the spectrum; count: the
// number of loci with each configuration); spectra gives the (0-based)
// groups of each spectrum, and n_proj the projection of each group. The
// hypergeometric probabilities are computed first (with R::dhyper()), and
// the spectra are then projected in parallel. Returns a list with, for each
// spectrum, its cells in column-major order (the first population changes
// fastest), the number of loci used and the number of loci dropped (with
// fewer than n_proj called copies in at least one population).

// [[Rcpp::export]]
List sfs_project_cpp(const List& configs,
                     const List& spectra,
                     const IntegerVector& n_proj,
                     int ncores) {
  size_t n_spectra = spectra.size();
  if ((size_t) configs.size() != n_spectra) {
    stop("configs and spectra should have the same length");
  }
  size_t ngroups = n_proj.size();
  for (size_t g = 0; g < ngroups; g++) {
    if (n_proj[g] < 0) stop("n_proj should not be negative");
  }

  // the probabilities of each configuration (n, k) of each group, in a flat
  // storage
  std::vector<std::unordered_map<uint64_t, SfsProjVec>> tables(ngroups);
  std::vector<double> probs;
  std::vector<SfsSpectrumData> data(n_spectra);

  // first pass: the probabilities needed, and the plain data of each
  // spectrum (without pointers, as the tables can still be rehashed)
  std::vector<std::vector<int>> spectrum_groups(n_spectra);
  for (size_t s = 0; s < n_spectra; s++) {
    IntegerVector groups = spectra[s];
    List conf = configs[s];
    IntegerMatrix n = conf["n"], k = conf["k"];
    NumericVector count = conf["count"];
    size_t K = groups.size();
    if (K < 1 || (size_t) n.ncol() != K || (size_t) k.ncol() != K ||
        n.nrow() != count.size() || k.nrow() != count.size()) {
      stop("invalid configurations");
    }
    SfsSpectrumData& sd = data[s];
    sd.K = K;
    sd.p.resize(K);
    sd.stride.resize(K);
    double n_cells = 1;
    size_t stride = 1;
    for (size_t i = 0; i < K; i++) {
      int g = groups[i];
      if (g < 0 || (size_t) g >= ngroups) stop("invalid group in spectra");
      spectrum_groups[s].push_back(g);
      sd.p[i] = n_proj[g];
      sd.stride[i] = stride;
      stride *= (size_t) n_proj[g] + 1;
      n_cells *= (double) n_proj[g] + 1;
    }
    if (n_cells > 1e10) stop("the spectrum has too many cells");
    sd.n_used = 0;
    sd.n_dropped = 0;
    for (size_t r = 0; r < (size_t) count.size(); r++) {
      bool keep = true;
      for (size_t i = 0; i < K; i++) {
        if (n(r, i) < sd.p[i]) {
          keep = false;
          break;
        }
      }
      if (!keep) {
        sd.n_dropped += count[r];
        continue;
      }
      for (size_t i = 0; i < K; i++) {
        int nn = n(r, i), kk = k(r, i), p = sd.p[i];
        if (kk < 0 || kk > nn) stop("invalid configuration");
        uint64_t key = sfs_tri(nn, kk);
        auto& table = tables[spectrum_groups[s][i]];
        if (table.find(key) == table.end()) {
          SfsProjVec v;
          v.lo = std::max(0, p - (nn - kk));
          v.hi = std::min(kk, p);
          v.offset = probs.size();
          for (int j = v.lo; j <= v.hi; j++) {
            probs.push_back(R::dhyper(j, kk, nn - kk, p, 0));
          }
          table.emplace(key, v);
        }
        sd.keys.push_back(key);
      }
      sd.count.push_back(count[r]);
      sd.n_used += count[r];
    }
  }
  // pointers to the probabilities of each row (the tables are now fixed)
  for (size_t s = 0; s < n_spectra; s++) {
    SfsSpectrumData& sd = data[s];
    size_t n_rows = sd.count.size();
    sd.vecs.resize(n_rows * sd.K);
    for (size_t r = 0; r < n_rows; r++) {
      for (size_t i = 0; i < sd.K; i++) {
        sd.vecs[r * sd.K + i] =
          &tables[spectrum_groups[s][i]].at(sd.keys[r * sd.K + i]);
      }
    }
  }

  // project the spectra in parallel
  const double* p_probs = probs.data();
#pragma omp parallel for num_threads(ncores) schedule(dynamic, 1)
  for (size_t s = 0; s < n_spectra; s++) {
    sfs_contract(data[s], p_probs);
  }

  List res(n_spectra);
  for (size_t s = 0; s < n_spectra; s++) {
    res[s] = List::create(
      _["sfs"] = NumericVector(data[s].out.begin(), data[s].out.end()),
      _["n_loci"] = data[s].n_used,
      _["n_dropped"] = data[s].n_dropped
    );
  }
  return res;
}

/******************************************************************************/

// The number of segregating sites in the 1D spectrum projected to each
// number of allele copies p from 0 to n_max, given the configurations of a
// population (n, k and count, as returned by sfs_config_counts_cpp()). A
// locus with n >= p called copies, k of which are alternate, is segregating
// in the projected sample with probability
// 1 - C(n - k, p) / C(n, p) - C(k, p) / C(n, p) (the probabilities of
// drawing only reference or only alternate copies), which is the sum of its
// contributions to the cells of the projected spectrum that are not
// corners. Also returns the number of loci with at least p called copies.
// The ratios of binomial coefficients are computed iteratively over p.

// [[Rcpp::export]]
List sfs_preview_cpp(const IntegerVector& n,
                     const IntegerVector& k,
                     const NumericVector& count,
                     int n_max) {
  size_t n_configs = count.size();
  if ((size_t) n.size() != n_configs || (size_t) k.size() != n_configs) {
    stop("n, k and count should have the same length");
  }
  if (n_max < 0) stop("n_max should not be negative");
  std::vector<long double> seg(n_max + 2, 0.0), loci(n_max + 2, 0.0);
  for (size_t r = 0; r < n_configs; r++) {
    int nn = n[r], kk = k[r];
    double c = count[r];
    if (nn < 1) continue;
    if (nn > n_max || kk < 0 || kk > nn) stop("invalid configuration");
    // the locus counts towards every p from 1 to nn
    loci[1] += c;
    loci[nn + 1] -= c;
    // a monomorphic locus is never segregating
    if (kk == 0 || kk == nn) continue;
    double h_ref = 1, h_alt = 1;
    for (int p = 1; p <= nn; p++) {
      // C(n - k, p) / C(n, p) and C(k, p) / C(n, p) (they reach 0, and stay
      // there, once p is larger than n - k or k)
      h_ref *= (double) (nn - kk - p + 1) / (double) (nn - p + 1);
      h_alt *= (double) (kk - p + 1) / (double) (nn - p + 1);
      if (h_ref < 0) h_ref = 0;
      if (h_alt < 0) h_alt = 0;
      seg[p] += c * (1 - h_ref - h_alt);
    }
  }
  NumericVector n_seg(n_max + 1), n_loci(n_max + 1);
  long double cum = 0;
  for (int p = 0; p <= n_max; p++) {
    cum += loci[p];
    n_loci[p] = (double) cum;
    n_seg[p] = (double) seg[p];
  }
  return List::create(_["n_seg_sites"] = n_seg, _["n_loci"] = n_loci);
}

/******************************************************************************/
