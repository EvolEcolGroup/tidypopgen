// Prototype: ways of decompressing a BGZF (bgzip) file, each returning the
// number of uncompressed bytes and of newlines (to check that they agree).
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <zlib.h>
#include <libdeflate.h>
#include <cstdio>
#include <cstring>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static size_t count_newlines(const char* p, size_t n) {
  size_t count = 0;
  const char* end = p + n;
  while ((p = static_cast<const char*>(std::memchr(p, '\n', end - p)))) {
    count++;
    p++;
  }
  return count;
}

// 1. zlib streaming (as gzgets does now), reading 4 MB at a time
// [[Rcpp::export]]
NumericVector zlib_stream(std::string filename) {
  gzFile f = gzopen(filename.c_str(), "rb");
  gzbuffer(f, 1 << 17);
  std::vector<char> buf(1 << 22);
  double bytes = 0, lines = 0;
  int n;
  while ((n = gzread(f, buf.data(), buf.size())) > 0) {
    bytes += n;
    lines += count_newlines(buf.data(), n);
  }
  gzclose(f);
  return NumericVector::create(bytes, lines);
}

// read the next BGZF block from f into block; returns its size (0 at EOF)
static size_t read_bgzf_block(FILE* f, std::vector<unsigned char> &block) {
  unsigned char header[18];
  size_t n = std::fread(header, 1, 18, f);
  if (n == 0) return 0;
  if (n < 18 || header[0] != 31 || header[1] != 139 || !(header[3] & 4) ||
      header[12] != 'B' || header[13] != 'C') {
    stop("not a BGZF file");
  }
  size_t bsize = (header[16] | (header[17] << 8)) + 1;
  block.resize(bsize);
  std::memcpy(block.data(), header, 18);
  if (std::fread(block.data() + 18, 1, bsize - 18, f) != bsize - 18) {
    stop("truncated BGZF block");
  }
  return bsize;
}

// 2./3. read batches of BGZF blocks and inflate them in parallel, with zlib
// (raw inflate) or libdeflate
// [[Rcpp::export]]
NumericVector bgzf_parallel(std::string filename, int n_threads,
                            bool use_libdeflate, int blocks_per_batch = 256) {
  FILE* f = std::fopen(filename.c_str(), "rb");
  if (!f) stop("cannot open");
  std::vector<std::vector<unsigned char> > blocks(blocks_per_batch);
  std::vector<size_t> out_offset(blocks_per_batch + 1);
  std::vector<char> out;
  double bytes = 0, lines = 0;
  bool error = false;
  while (true) {
    int n_blocks = 0;
    out_offset[0] = 0;
    for (; n_blocks < blocks_per_batch; n_blocks++) {
      size_t bsize = read_bgzf_block(f, blocks[n_blocks]);
      if (bsize == 0) break;
      // uncompressed size: last 4 bytes of the block (ISIZE)
      const unsigned char* t = blocks[n_blocks].data() + bsize - 4;
      size_t isize = t[0] | (t[1] << 8) | (t[2] << 16) | ((size_t)t[3] << 24);
      out_offset[n_blocks + 1] = out_offset[n_blocks] + isize;
    }
    if (n_blocks == 0) break;
    out.resize(out_offset[n_blocks]);
#pragma omp parallel num_threads(n_threads)
    {
      struct libdeflate_decompressor* d =
        use_libdeflate ? libdeflate_alloc_decompressor() : NULL;
      z_stream zs;
      if (!use_libdeflate) {
        std::memset(&zs, 0, sizeof(zs));
        inflateInit2(&zs, -15);
      }
#pragma omp for schedule(dynamic, 8)
      for (int i = 0; i < n_blocks; i++) {
        std::vector<unsigned char> &b = blocks[i];
        size_t xlen = b[10] | (b[11] << 8);
        size_t cdata_start = 12 + xlen;
        size_t cdata_len = b.size() - cdata_start - 8;
        size_t isize = out_offset[i + 1] - out_offset[i];
        if (isize == 0) continue;
        if (use_libdeflate) {
          size_t actual;
          if (libdeflate_deflate_decompress(d, b.data() + cdata_start, cdata_len,
              &out[out_offset[i]], isize, &actual) != LIBDEFLATE_SUCCESS ||
              actual != isize) {
            error = true;
          }
        } else {
          inflateReset(&zs);
          zs.next_in = b.data() + cdata_start;
          zs.avail_in = cdata_len;
          zs.next_out = reinterpret_cast<Bytef*>(&out[out_offset[i]]);
          zs.avail_out = isize;
          if (inflate(&zs, Z_FINISH) != Z_STREAM_END) error = true;
        }
      }
      if (use_libdeflate) libdeflate_free_decompressor(d);
      else inflateEnd(&zs);
    }
    bytes += out.size();
    lines += count_newlines(out.data(), out.size());
  }
  std::fclose(f);
  if (error) stop("decompression error");
  return NumericVector::create(bytes, lines);
}

// CRC32 of the whole (uncompressed) file, a 64 KB block at a time, with zlib
// or libdeflate, in parallel over the blocks of each batch
// [[Rcpp::export]]
double crc_cost(std::string filename, int n_threads, bool use_libdeflate) {
  gzFile f = gzopen(filename.c_str(), "rb");
  gzbuffer(f, 1 << 17);
  std::vector<unsigned char> buf(1 << 24);
  double total = 0;
  int n;
  double secs = 0;
  while ((n = gzread(f, buf.data(), buf.size())) > 0) {
    int n_blocks = (n + 65279) / 65280;
    double t0 = omp_get_wtime();
    unsigned long acc = 0;
#pragma omp parallel for num_threads(n_threads) reduction(^:acc)
    for (int i = 0; i < n_blocks; i++) {
      size_t start = (size_t) i * 65280;
      size_t len = std::min((size_t) 65280, (size_t) n - start);
      acc ^= use_libdeflate ? libdeflate_crc32(0, &buf[start], len) :
        crc32(0L, &buf[start], len);
    }
    secs += omp_get_wtime() - t0;
    total += acc;
  }
  gzclose(f);
  return secs;
}
