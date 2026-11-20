/******************************************************************************/

// read_bed_to_fbm() is adapted from readbina() in src/read-plink.cpp of the R
// package bigsnpr (https://github.com/privefl/bigsnpr), copyright (C) Florian
// Prive, released under the GNU General Public License version 3 (GPL-3).
// readbina() reads a PLINK .bed file into an FBM.code256 one locus at a time,
// on a single core; here the file is read in large chunks, each chunk is
// decoded in parallel (while the next chunk is read), the size of the file is
// checked before reading it, and all three bytes of its magic number are
// checked.

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <algorithm>

using namespace Rcpp;

/******************************************************************************/

// Decode the loci [j0, j1) of a chunk of the .bed file (whose first locus is
// chunk_j0, and whose bytes are in chunk) into their columns of the FBM.
// Each byte of the .bed holds the genotypes of 4 individuals (2 bits each,
// the first individual in the lowest bits); lut holds the 4 bytes of the FBM
// for each of the 256 values of a byte (from the code table of bigsnpr), so
// that each byte is decoded with a lookup and a 4-byte copy.
static inline void decode_bed_loci(const unsigned char* chunk,
                                   std::size_t chunk_j0, std::size_t j0,
                                   std::size_t j1, std::size_t n,
                                   std::size_t bytes_per_locus,
                                   const uint32_t* lut, unsigned char* fbm) {
  const std::size_t n_full = n / 4;       // bytes with 4 individuals
  const std::size_t extra = n - 4 * n_full; // individuals in the last byte
  for (std::size_t j = j0; j < j1; j++) {
    const unsigned char* src = chunk + (j - chunk_j0) * bytes_per_locus;
    unsigned char* dst = fbm + j * n;
    for (std::size_t k = 0; k < n_full; k++) {
      std::memcpy(dst + 4 * k, &lut[src[k]], 4);
    }
    if (extra) {
      // (the padding bits of the last byte are ignored)
      std::memcpy(dst + 4 * n_full, &lut[src[n_full]], extra);
    }
  }
}

/******************************************************************************/

// Read a PLINK .bed file (SNP-major) into an FBM.code256 (with as many rows
// as individuals and as many columns as loci), using n_threads threads.
// Returns true if the end of the file was reached (i.e. the file holds no
// more bytes than the loci of the FBM need), false otherwise (the caller
// gives a warning, as for bigsnpr). A file with fewer bytes than needed is an
// error.
//
// The file is read in chunks of loci (of ~64MB); while the threads decode a
// chunk into the FBM, one of them reads the next chunk into a second buffer
// (and then joins the decoding), so that reading from the disk overlaps with
// decoding. The loci of a chunk are decoded in blocks (of ~1MB in the FBM),
// handed out to the threads in order, so that the FBM is written roughly
// from start to end.
// [[Rcpp::export]]
bool read_bed_to_fbm(std::string filename, Environment BM,
                     const RawMatrix& tab, int n_threads) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* fbm = static_cast<unsigned char*>(xpBM->matrix());
  const std::size_t n = xpBM->nrow();  // individuals
  const std::size_t m = xpBM->ncol();  // loci
  const std::size_t bytes_per_locus = (n + 3) / 4;

  if ((tab.nrow() != 4) || (tab.ncol() != 256)) {
    stop("the code table should have 4 rows and 256 columns");
  }
  if (n_threads < 1) {
    n_threads = 1;
  }

  std::ifstream bed(filename.c_str(), std::ios::in | std::ios::binary);
  if (!bed) {
    stop("could not open the .bed file " + filename);
  }
  // the size of the file must be at least that of the magic number and of
  // the loci (a smaller file would be truncated or corrupted)
  bed.seekg(0, std::ios::end);
  const unsigned long long file_size =
    static_cast<unsigned long long>(bed.tellg());
  const unsigned long long needed =
    3ULL + static_cast<unsigned long long>(m) * bytes_per_locus;
  if (file_size < needed) {
    stop("the .bed file " + filename + " is smaller than expected for the " +
         "number of individuals in the .fam file and of loci in the .bim " +
         "file; it might be truncated or corrupted");
  }
  bed.seekg(0, std::ios::beg);
  // magic number: 0x6C 0x1B, then 0x01 for a SNP-major file
  unsigned char magic[3];
  bed.read(reinterpret_cast<char*>(magic), 3);
  if (!bed || (magic[0] != 0x6C) || (magic[1] != 0x1B)) {
    stop("the file " + filename + " is not a PLINK .bed file (wrong magic " +
         "number)");
  }
  if (magic[2] != 0x01) {
    stop("the .bed file " + filename + " is in individual-major mode, which " +
         "is not supported; convert it to SNP-major mode with PLINK (e.g. " +
         "plink --bfile <prefix> --make-bed --out <new prefix>)");
  }

  // the 4 bytes of the FBM for each value of a byte of the .bed (column b of
  // tab holds the codes of the 4 individuals of byte b)
  uint32_t lut[256];
  for (int b = 0; b < 256; b++) {
    std::memcpy(&lut[b], &tab(0, b), 4);
  }

  // loci per chunk (~64MB of the .bed) and per block (~1MB of the FBM)
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
    bed.read(reinterpret_cast<char*>(b.data()),
             static_cast<std::streamsize>(n_loci * bytes_per_locus));
    return static_cast<bool>(bed);
  };

  if (m > 0 && !read_chunk(0)) {
    stop("error reading the .bed file " + filename);
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
      // one thread reads the next chunk (into the other buffer), and then
      // joins the others in decoding this one
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
        decode_bed_loci(chunk, chunk_j0, j0, j1, n, bytes_per_locus, lut,
                        fbm);
      }
    }
    if (!read_ok) {
      stop("error reading the .bed file " + filename);
    }
    Rcpp::checkUserInterrupt();
  }

  return file_size == needed;
}

/******************************************************************************/
