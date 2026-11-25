/******************************************************************************/

// read_packedancestry_to_fbm() reads the genotypes of a packedancestrymap
// file of EIGENSOFT into an FBM.code256. Its reading of GENO files
// (SNP-major) is adapted from readbina() in src/read-plink.cpp of the R
// package bigsnpr (https://github.com/privefl/bigsnpr), copyright (C) Florian
// Prive, released under the GNU General Public License version 3 (GPL-3),
// through read_bed_to_fbm() (src/read_bed.cpp): the file is read in large
// chunks, each decoded in parallel while the next one is read. TGENO files
// (individual-major) are transposed into the FBM in parallel, in blocks that
// fit in the cache.

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <memory>
#if defined(__linux__)
#include <sys/mman.h>
#endif

using namespace Rcpp;

/******************************************************************************/

// A packedancestrymap genotype file starts with a header, a string such as
// "GENO 5 16 be181b4d 9f75a306" (the numbers of individuals and of loci, and
// hashes of their names), padded with NUL bytes. Then:
// - GENO (SNP-major): one record per locus, with the genotypes of all the
//   individuals; records (and the header) have max(48, ceiling(n / 4)) bytes.
// - TGENO (individual-major, or "transposed"): one record per individual,
//   with its genotypes at all the loci; the header has 48 bytes, and records
//   ceiling(m / 4) bytes (with no minimum).
// In both, each byte holds 4 genotypes of 2 bits (the number of copies of the
// first allele of the .snp file, or 3 if missing), the first genotype in the
// highest bits. The code table of get_packedancestry_code() gives the 4 bytes
// of the FBM for each value of a byte (with the dosage of the second allele).
// (See geno.c in https://github.com/DReichLab/EIG.)

static const std::size_t TGENO_HEADER_BYTES = 48;
static const std::size_t GENO_MIN_RECORD_BYTES = 48;

/******************************************************************************/

// Decode the loci [j0, j1) of a chunk of a GENO file (whose first locus is
// chunk_j0, and whose records have bytes_per_locus bytes) into their columns
// of the FBM (n individuals), with a lookup and a 4-byte copy per byte.
static inline void decode_geno_loci(const unsigned char* chunk,
                                    std::size_t chunk_j0, std::size_t j0,
                                    std::size_t j1, std::size_t n,
                                    std::size_t bytes_per_locus,
                                    const unsigned char (*lut)[4],
                                    unsigned char* fbm) {
  const std::size_t n_full = n / 4;         // bytes with 4 individuals
  const std::size_t extra = n - 4 * n_full; // individuals in the last byte
  for (std::size_t j = j0; j < j1; j++) {
    const unsigned char* src = chunk + (j - chunk_j0) * bytes_per_locus;
    unsigned char* dst = fbm + j * n;
    for (std::size_t k = 0; k < n_full; k++) {
      std::memcpy(dst + 4 * k, lut[src[k]], 4);
    }
    if (extra) {
      // (the padding bits of the last byte are ignored)
      std::memcpy(dst + 4 * n_full, lut[src[n_full]], extra);
    }
  }
}

// Read the m loci of a GENO file (positioned after its header) into the FBM.
// The file is read in chunks of loci (of ~64MB); while the threads decode a
// chunk, one of them reads the next chunk into a second buffer (and then
// joins the decoding). The loci of a chunk are decoded in blocks (of ~1MB in
// the FBM), handed out to the threads in order, so that the FBM is written
// roughly from start to end.
static void read_geno_records(std::ifstream& geno, const std::string& filename,
                              std::size_t n, std::size_t m,
                              std::size_t bytes_per_locus,
                              const unsigned char (*lut)[4],
                              unsigned char* fbm, int n_threads) {
  const std::size_t chunk_loci = std::max<std::size_t>(
    1, std::min<std::size_t>(m, (64u << 20) / bytes_per_locus));
  const std::size_t block_loci =
    std::max<std::size_t>(1, (1u << 20) / std::max<std::size_t>(n, 1));
  const std::size_t n_chunks = (m + chunk_loci - 1) / chunk_loci;

  std::vector<unsigned char> buf[2];
  buf[0].resize(chunk_loci * bytes_per_locus);
  if (n_chunks > 1) {
    buf[1].resize(chunk_loci * bytes_per_locus);
  }

  // read chunk c into its buffer (returns false on a short read)
  auto read_chunk = [&](std::size_t c) -> bool {
    std::size_t j0 = c * chunk_loci;
    std::size_t n_loci = std::min(chunk_loci, m - j0);
    std::vector<unsigned char>& b = buf[c & 1];
    geno.read(reinterpret_cast<char*>(b.data()),
              static_cast<std::streamsize>(n_loci * bytes_per_locus));
    return static_cast<bool>(geno);
  };

  if (m > 0 && !read_chunk(0)) {
    stop("error reading the .geno file " + filename);
  }
  for (std::size_t c = 0; c < n_chunks; c++) {
    const std::size_t chunk_j0 = c * chunk_loci;
    const std::size_t chunk_j1 = std::min(m, chunk_j0 + chunk_loci);
    const std::size_t n_blocks =
      (chunk_j1 - chunk_j0 + block_loci - 1) / block_loci;
    const unsigned char* chunk = buf[c & 1].data();
    bool read_ok = true;
#pragma omp parallel num_threads(n_threads)
{
  // one thread reads the next chunk (into the other buffer), and then joins
  // the others in decoding this one
#pragma omp single nowait
{
  if (c + 1 < n_chunks) {
    read_ok = read_chunk(c + 1);
  }
}
#pragma omp for schedule(dynamic, 1)
for (long blk = 0; blk < static_cast<long>(n_blocks); blk++) {
  std::size_t j0 = chunk_j0 + blk * block_loci;
  std::size_t j1 = std::min(chunk_j1, j0 + block_loci);
  decode_geno_loci(chunk, chunk_j0, j0, j1, n, bytes_per_locus, lut, fbm);
}
}
    if (!read_ok) {
      stop("error reading the .geno file " + filename);
    }
    Rcpp::checkUserInterrupt();
  }
}

/******************************************************************************/

// On little-endian platforms, the genotypes of a TGENO file are transposed
// in 64-bit words ("SIMD within a register"): the bytes of 8 individuals at 8
// bytes of their records (32 loci) are loaded into 8 words, their 2-bit
// genotypes are recoded, the 8 x 8 bytes are transposed, and the 4 genotypes
// of each byte are written to the FBM 8 individuals at a time. This takes
// about 40% less CPU time than decoding one byte at a time with the code
// table.
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define TGENO_WORDS 1
#endif
#endif
#ifndef TGENO_WORDS
#define TGENO_WORDS 0
#endif

#if defined(__GNUC__) || defined(__clang__)
#define TGENO_PREFETCH(p) __builtin_prefetch(p)
#else
#define TGENO_PREFETCH(p)
#endif

// Recode each 2-bit genotype of a word as codes[genotype] (each code is at
// most 3, so that it fits in the 2 bits)
static inline uint64_t recode_genotypes(uint64_t x, const uint64_t* codes) {
  const uint64_t LOW = 0x5555555555555555ULL; // the low bit of each genotype
  const uint64_t lo = x & LOW;
  const uint64_t hi = (x >> 1) & LOW;
  // (for each genotype, exactly one of the four indicators is 1)
  return (~(lo | hi) & LOW) * codes[0] + (lo & ~hi) * codes[1] +
    (hi & ~lo) * codes[2] + (lo & hi) * codes[3];
}

// Transpose the 8 x 8 bytes of 8 words: byte c of word r becomes byte r of
// word c (swapping blocks of 4 x 4, then 2 x 2, then 1 x 1 bytes)
static inline void transpose_8x8_bytes(uint64_t* a) {
  for (int r = 0; r < 4; r++) {
    uint64_t t = ((a[r] >> 32) ^ a[r + 4]) & 0x00000000FFFFFFFFULL;
    a[r] ^= t << 32;
    a[r + 4] ^= t;
  }
  static const int first_rows[4] = {0, 1, 4, 5};
  for (int p = 0; p < 4; p++) {
    const int r = first_rows[p];
    uint64_t t = ((a[r] >> 16) ^ a[r + 2]) & 0x0000FFFF0000FFFFULL;
    a[r] ^= t << 16;
    a[r + 2] ^= t;
  }
  for (int r = 0; r < 8; r += 2) {
    uint64_t t = ((a[r] >> 8) ^ a[r + 1]) & 0x00FF00FF00FF00FFULL;
    a[r] ^= t << 8;
    a[r + 1] ^= t;
  }
}

// Transpose 8 bytes (32 loci) of the records of 8 individuals: src points to
// the first byte of the first individual (whose records are slab_bytes
// apart), and dst to the FBM, at the first individual of the first locus
static inline void transpose_8x8_block(const unsigned char* src,
                                       std::size_t slab_bytes,
                                       const uint64_t* codes,
                                       unsigned char* dst, std::size_t n) {
  const uint64_t LOW2 = 0x0303030303030303ULL; // one genotype per byte
  uint64_t a[8];
  for (int r = 0; r < 8; r++) {
    std::memcpy(&a[r], src + r * slab_bytes, 8);
    a[r] = recode_genotypes(a[r], codes);
  }
  transpose_8x8_bytes(a);
  // word c holds byte c of the 8 individuals, whose 4 genotypes (the first
  // in the highest bits) are the loci 4c to 4c + 3
  for (int c = 0; c < 8; c++) {
    for (int k = 0; k < 4; k++) {
      uint64_t g = (a[c] >> (6 - 2 * k)) & LOW2;
      std::memcpy(dst + (4 * c + k) * n, &g, 8);
    }
  }
}

// Transpose the bytes [b0, b1) of the records of a slab of a TGENO file into
// the FBM (n individuals, m loci). The slab holds, for each of the n
// individuals in turn, slab_bytes bytes of its record, starting at byte
// first_byte of the record; byte b of a record holds the loci 4b to 4b + 3.
// codes gives the code of the FBM for each 2-bit genotype, if the code table
// is made of them, and is NULL otherwise.
// The bytes are transposed in blocks of 64 bytes (256 loci) by 512
// individuals; within a block, 8 bytes (32 loci of the FBM) at a time for all
// the individuals, prefetching the records of the individuals to come (whose
// bytes are far apart, so that the prefetchers of the CPU miss them). The
// individuals left over (fewer than 8) and the bytes left over (fewer than 8
// in a block), and all the bytes if codes is NULL or the platform is
// big-endian, are decoded one byte at a time with the code table.
static inline void transpose_tgeno_bytes(const unsigned char* slab,
                                         std::size_t slab_bytes,
                                         std::size_t first_byte,
                                         std::size_t b0, std::size_t b1,
                                         std::size_t n, std::size_t m,
                                         const unsigned char (*lut)[4],
                                         const uint64_t* codes,
                                         unsigned char* fbm) {
  const std::size_t BLOCK_BYTES = 64;
  const std::size_t BLOCK_INDIVS = 512;
  const std::size_t PREFETCH_INDIVS = 32;
  // the bytes whose 4 loci all exist (the last byte of a record can be
  // partly padding)
  const std::size_t full_bytes = m / 4;
  const std::size_t b1_full = std::min(b1, (full_bytes > first_byte) ?
                                       full_bytes - first_byte : 0);
  // individuals transposed 8 at a time
  const std::size_t n8 = (TGENO_WORDS && codes != NULL) ? n - n % 8 : 0;
  for (std::size_t bb0 = b0; bb0 < b1_full; bb0 += BLOCK_BYTES) {
    const std::size_t bb1 = std::min(b1_full, bb0 + BLOCK_BYTES);
    // bytes transposed 8 at a time
    const std::size_t bb8 = (n8 > 0) ? bb0 + ((bb1 - bb0) / 8) * 8 : bb0;
    for (std::size_t i0 = 0; i0 < n8; i0 += BLOCK_INDIVS) {
      const std::size_t i1 = std::min(n8, i0 + BLOCK_INDIVS);
      for (std::size_t b = bb0; b < bb8; b += 8) {
        unsigned char* dst = fbm + 4 * (first_byte + b) * n;
        for (std::size_t i = i0; i < i1; i += 8) {
          const unsigned char* src = slab + i * slab_bytes + b;
          if (i + PREFETCH_INDIVS + 8 <= n) {
            for (std::size_t r = 0; r < 8; r++) {
              TGENO_PREFETCH(src + (PREFETCH_INDIVS + r) * slab_bytes);
            }
          }
          transpose_8x8_block(src, slab_bytes, codes, dst + i, n);
        }
      }
    }
    // the rest, one byte at a time
    for (std::size_t b = bb0; b < bb1; b++) {
      const unsigned char* src = slab + b;
      unsigned char* d0 = fbm + 4 * (first_byte + b) * n;
      unsigned char* d1 = d0 + n;
      unsigned char* d2 = d1 + n;
      unsigned char* d3 = d2 + n;
      for (std::size_t i = (b < bb8) ? n8 : 0; i < n; i++) {
        const unsigned char* g = lut[src[i * slab_bytes]];
        d0[i] = g[0];
        d1[i] = g[1];
        d2[i] = g[2];
        d3[i] = g[3];
      }
    }
  }
  // the last byte of the records, if it holds fewer than 4 loci
  if (b1_full < b1) {
    const std::size_t b = b1_full;
    const std::size_t j0 = 4 * (first_byte + b);
    for (std::size_t i = 0; i < n; i++) {
      const unsigned char* g = lut[slab[i * slab_bytes + b]];
      for (std::size_t j = j0; j < m; j++) {
        fbm[j * n + i] = g[j - j0];
      }
    }
  }
}

// A buffer of bytes (not initialised), aligned to 2MB; on Linux, it is backed
// by transparent huge pages where the system allows it (e.g. with
// /sys/kernel/mm/transparent_hugepage/enabled set to "madvise", the default
// of Ubuntu), so that the transposition of a TGENO file, which reads its
// records far apart, makes fewer TLB misses (it then takes about 20% less
// CPU time).
class HugePageBuffer {
public:
  explicit HugePageBuffer(std::size_t size)
    : raw_(new unsigned char[size + ALIGN]) {
    std::uintptr_t p = reinterpret_cast<std::uintptr_t>(raw_.get());
    data_ = reinterpret_cast<unsigned char*>((p + ALIGN - 1) & ~(ALIGN - 1));
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (size > 0) {
      madvise(data_, size, MADV_HUGEPAGE); // (a failure is harmless)
    }
#endif
  }
  unsigned char* data() { return data_; }
private:
  static const std::size_t ALIGN = std::size_t(1) << 21;
  std::unique_ptr<unsigned char[]> raw_;
  unsigned char* data_;
};

// Read the n records of a TGENO file (positioned after its header) into the
// FBM, i.e. transpose them. The FBM is written in slabs of loci, in order
// (so that an FBM larger than the memory is written sequentially), each
// slab from a buffer with the bytes of those loci for all the individuals.
// If the records fit in max_buffer_bytes, there is a single slab, read with
// one sequential read of the whole file (the fastest on any disk). Otherwise
// each slab needs a read for each individual, and the next slab is read (by
// one thread, into a second buffer of the same size) while the others decode
// this one. The bytes of a slab are transposed in parallel, in tasks of at
// least ~1MB of the FBM, handed out in order.
static void read_tgeno_records(std::ifstream& geno, const std::string& filename,
                               std::size_t n, std::size_t m,
                               std::size_t bytes_per_indiv,
                               const unsigned char (*lut)[4],
                               const uint64_t* codes,
                               unsigned char* fbm, int n_threads,
                               double max_buffer_bytes) {
  if (n == 0 || m == 0) {
    return;
  }
  // bytes of each record in a slab (a multiple of 64 unless it is the
  // whole record)
  std::size_t slab_bytes = bytes_per_indiv;
  if (static_cast<double>(n) * static_cast<double>(bytes_per_indiv) >
        max_buffer_bytes) {
    double per_indiv = max_buffer_bytes / 2 / static_cast<double>(n);
    slab_bytes = std::max<std::size_t>(
      64, static_cast<std::size_t>(per_indiv / 64) * 64);
    slab_bytes = std::min(slab_bytes, bytes_per_indiv);
  }
  const std::size_t n_slabs =
    (bytes_per_indiv + slab_bytes - 1) / slab_bytes;
  // bytes of the records per task (at least ~1MB of the FBM)
  const std::size_t task_bytes = std::max<std::size_t>(
    64, ((((1u << 20) / (4 * n)) + 63) / 64) * 64);

  std::unique_ptr<HugePageBuffer> buf[2];
  buf[0].reset(new HugePageBuffer(n * slab_bytes));
  if (n_slabs > 1) {
    buf[1].reset(new HugePageBuffer(n * slab_bytes));
  }
  const std::streamoff records_start = geno.tellg();

  // read slab s into its buffer (returns false on a short read)
  auto read_slab = [&](std::size_t s) -> bool {
    HugePageBuffer& b = *buf[s & 1];
    if (n_slabs == 1) {
      // the records, one after the other
      geno.read(reinterpret_cast<char*>(b.data()),
                static_cast<std::streamsize>(n * bytes_per_indiv));
      return static_cast<bool>(geno);
    }
    const std::size_t first_byte = s * slab_bytes;
    const std::size_t n_bytes =
      std::min(slab_bytes, bytes_per_indiv - first_byte);
    for (std::size_t i = 0; i < n; i++) {
      geno.seekg(records_start +
        static_cast<std::streamoff>(i * bytes_per_indiv + first_byte));
      geno.read(reinterpret_cast<char*>(b.data() + i * slab_bytes),
                static_cast<std::streamsize>(n_bytes));
      if (!geno) {
        return false;
      }
    }
    return true;
  };

  if (!read_slab(0)) {
    stop("error reading the .geno file " + filename);
  }
  for (std::size_t s = 0; s < n_slabs; s++) {
    const std::size_t first_byte = s * slab_bytes;
    const std::size_t n_bytes =
      std::min(slab_bytes, bytes_per_indiv - first_byte);
    const std::size_t n_tasks = (n_bytes + task_bytes - 1) / task_bytes;
    const unsigned char* slab = buf[s & 1]->data();
    bool read_ok = true;
#pragma omp parallel num_threads(n_threads)
{
#pragma omp single nowait
{
  if (s + 1 < n_slabs) {
    read_ok = read_slab(s + 1);
  }
}
#pragma omp for schedule(dynamic, 1)
for (long t = 0; t < static_cast<long>(n_tasks); t++) {
  std::size_t b0 = t * task_bytes;
  std::size_t b1 = std::min(n_bytes, b0 + task_bytes);
  transpose_tgeno_bytes(slab, slab_bytes, first_byte, b0, b1, n, m, lut, codes,
                        fbm);
}
}
    if (!read_ok) {
      stop("error reading the .geno file " + filename);
    }
    Rcpp::checkUserInterrupt();
  }
}

/******************************************************************************/

// Read a packedancestrymap genotype file (GENO or TGENO, as transposed says)
// into an FBM.code256 (with as many rows as individuals and as many columns
// as loci), using n_threads threads. For a TGENO file, the records are read
// into memory in one go if they fit in max_buffer_bytes (and otherwise in
// slabs, using up to max_buffer_bytes).
// The caller (fbm_read_packedancestry()) checks the header and the size of the
// file before creating the FBM, so that a bad file leaves no backing file
// behind; the checks of the size are repeated here in case the file changed
// in between (a larger file is accepted: the extra bytes at its end are
// ignored).
// [[Rcpp::export]]
void read_packedancestry_to_fbm(std::string filename, Environment BM,
                                const RawMatrix& tab, bool transposed,
                                int n_threads, double max_buffer_bytes) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* fbm = static_cast<unsigned char*>(xpBM->matrix());
  const std::size_t n = xpBM->nrow();  // individuals
  const std::size_t m = xpBM->ncol();  // loci

  if ((tab.nrow() != 4) || (tab.ncol() != 256)) {
    stop("the code table should have 4 rows and 256 columns");
  }
  if (n_threads < 1) {
    n_threads = 1;
  }

  std::size_t header_bytes, record_bytes, n_records;
  if (transposed) {
    header_bytes = TGENO_HEADER_BYTES;
    record_bytes = (m + 3) / 4;
    n_records = n;
  } else {
    record_bytes = std::max(GENO_MIN_RECORD_BYTES, (n + 3) / 4);
    header_bytes = record_bytes;
    n_records = m;
  }

  std::ifstream geno(filename.c_str(), std::ios::in | std::ios::binary);
  if (!geno) {
    stop("could not open the .geno file " + filename);
  }
  geno.seekg(0, std::ios::end);
  const unsigned long long file_size =
    static_cast<unsigned long long>(geno.tellg());
  const unsigned long long needed = header_bytes +
    static_cast<unsigned long long>(n_records) * record_bytes;
  if (file_size < needed) {
    stop("the .geno file " + filename + " is smaller than expected for the " +
         "number of individuals in the .ind file and of loci in the .snp " +
         "file; it might be truncated or corrupted");
  }
  geno.seekg(0, std::ios::beg);
  char magic[5];
  geno.read(magic, 5);
  const bool magic_ok = transposed ? (std::memcmp(magic, "TGENO", 5) == 0) :
    (std::memcmp(magic, "GENO", 4) == 0);
  if (!geno || !magic_ok) {
    stop("the file " + filename + " is not a packedancestrymap " +
         (transposed ? "TGENO" : "GENO") + " file");
  }
  geno.seekg(static_cast<std::streamoff>(header_bytes), std::ios::beg);

  // the 4 bytes of the FBM for each value of a byte (column b of tab holds
  // the codes of the 4 genotypes of byte b)
  unsigned char lut[256][4];
  for (int b = 0; b < 256; b++) {
    std::memcpy(lut[b], &tab(0, b), 4);
  }

  if (transposed) {
    // the code of the FBM for each 2-bit genotype, if the code table is made
    // of them (as that of get_packedancestry_code()), with codes of at most 3
    // (as for CODE_012), so that TGENO files can be transposed in words
    uint64_t codes[4];
    bool word_codes = true;
    for (int v = 0; v < 4; v++) {
      codes[v] = lut[v << 6][0];
      word_codes = word_codes && (codes[v] <= 3);
    }
    for (int b = 0; b < 256; b++) {
      for (int k = 0; k < 4; k++) {
        word_codes = word_codes &&
          (lut[b][k] == codes[(b >> (6 - 2 * k)) & 3]);
      }
    }
    read_tgeno_records(geno, filename, n, m, record_bytes, lut,
                       word_codes ? codes : NULL, fbm, n_threads,
                       max_buffer_bytes);
  } else {
    read_geno_records(geno, filename, n, m, record_bytes, lut, fbm,
                      n_threads);
  }
}

/******************************************************************************/
