// What would it cost to bgzip the locus sidecar? A standalone benchmark of
// BGZF-style compression (independent deflate blocks of 65,280 bytes of
// input, as bgzip and htslib write them) of a real sidecar file, with zlib
// and libdeflate, against reading the plain file.
//
//   g++ -O2 -fopenmp -o sidecar_compress sidecar_compress.cpp -lz -ldeflate
//   ./sidecar_compress <sidecar.loci.tsv> <out_dir> [threads=4] [reps=3]
//                      [lookups=2000]
//
// For each codec (zlib and libdeflate, levels 1 and 6) it writes
// <out_dir>/sidecar.<codec>.gz (a valid BGZF file, with the EOF block: check
// it with `gzip -dc file | cmp - sidecar.loci.tsv`) and prints a TSV row
// with:
// - the compressed size;
// - compression time (CRC32 included, writing the file not included) with 1
//   thread and with `threads` threads (OpenMP over the blocks, which are
//   independent);
// - full decompression time (inflate and CRC32 check of every block, into
//   one buffer) with 1 thread and with `threads` threads, and the time to
//   decompress and scan every line (find the newlines and parse the integer
//   key of each line, which is what building a sparse index needs);
// - the cost of one random lookup: pread() of one random block of the file
//   (in the page cache) and its inflation with CRC32 check, averaged over
//   `lookups` blocks.
// The "plain" row gives the same for the uncompressed file: the scan of the
// lines, and a random pread() of 1 MB (what SidecarReader reads after each
// seek, see src/loci_sidecar.cpp) and of 64 KB.
// Times are the median of `reps` repeats.

#include <zlib.h>
#include <libdeflate.h>
#include <omp.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static const size_t BLOCK_IN = 0xff00;     // input per block (as htslib)
static const size_t BLOCK_MAX = 0x10000;   // largest BGZF block
static const size_t HEADER = 18;
static const size_t FOOTER = 8;
static const unsigned char BGZF_EOF[28] = {
  0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0,
  0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static double now() {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count();
}

static double median(std::vector<double> x) {
  std::sort(x.begin(), x.end());
  size_t n = x.size();
  return (n % 2) ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

static void fail(const std::string &msg) {
  std::fprintf(stderr, "error: %s\n", msg.c_str());
  std::exit(1);
}

static std::vector<unsigned char> read_file(const std::string &path) {
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) fail("cannot open " + path);
  std::vector<unsigned char> buf;
  std::fseek(f, 0, SEEK_END);
  long size = std::ftell(f);
  std::fseek(f, 0, SEEK_SET);
  buf.resize(size);
  if (std::fread(buf.data(), 1, size, f) != (size_t) size) fail("read " + path);
  std::fclose(f);
  return buf;
}

static void put16(unsigned char* p, uint32_t x) {
  p[0] = x & 0xff; p[1] = (x >> 8) & 0xff;
}
static void put32(unsigned char* p, uint32_t x) {
  p[0] = x & 0xff; p[1] = (x >> 8) & 0xff; p[2] = (x >> 16) & 0xff;
  p[3] = (x >> 24) & 0xff;
}
static uint32_t get32(const unsigned char* p) {
  return p[0] | (p[1] << 8) | (p[2] << 16) | ((uint32_t) p[3] << 24);
}

// scan the lines of text: count them and parse the integer key at the start
// of each (as sidecar_index() does); returns a checksum of the keys
static uint64_t scan_lines(const unsigned char* p, size_t n,
                           uint64_t &n_lines) {
  const unsigned char* end = p + n;
  uint64_t sum = 0;
  while (p < end) {
    int64_t key = 0;
    while (p < end && *p >= '0' && *p <= '9') key = key * 10 + (*p++ - '0');
    sum += key;
    const unsigned char* nl = (const unsigned char*) std::memchr(
      p, '\n', end - p);
    if (!nl) break;
    p = nl + 1;
    n_lines++;
  }
  return sum;
}

struct Codec {
  std::string name;
  bool libdeflate;
  int level;
};

// compress one block (in, len) into out (BLOCK_MAX bytes): returns the size
// of the BGZF block
struct Compressor {
  bool use_ld;
  z_stream zs;
  libdeflate_compressor* ld;
  Compressor(const Codec &c) : use_ld(c.libdeflate), ld(NULL) {
    if (use_ld) {
      ld = libdeflate_alloc_compressor(c.level);
      if (!ld) fail("libdeflate_alloc_compressor");
    } else {
      std::memset(&zs, 0, sizeof(zs));
      if (deflateInit2(&zs, c.level, Z_DEFLATED, -15, 8,
                       Z_DEFAULT_STRATEGY) != Z_OK) fail("deflateInit2");
    }
  }
  ~Compressor() {
    if (use_ld) libdeflate_free_compressor(ld); else deflateEnd(&zs);
  }
  size_t compress(const unsigned char* in, size_t len, unsigned char* out) {
    size_t cap = BLOCK_MAX - HEADER - FOOTER;
    size_t clen;
    uint32_t crc;
    if (use_ld) {
      clen = libdeflate_deflate_compress(ld, in, len, out + HEADER, cap);
      if (clen == 0) fail("block does not fit");
      crc = libdeflate_crc32(0, in, len);
    } else {
      deflateReset(&zs);
      zs.next_in = (Bytef*) in;
      zs.avail_in = len;
      zs.next_out = out + HEADER;
      zs.avail_out = cap;
      if (deflate(&zs, Z_FINISH) != Z_STREAM_END) fail("block does not fit");
      clen = cap - zs.avail_out;
      crc = crc32(crc32(0L, Z_NULL, 0), in, len);
    }
    static const unsigned char h[16] = {31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6,
                                        0, 'B', 'C', 2, 0};
    std::memcpy(out, h, 16);
    size_t bsize = HEADER + clen + FOOTER;
    put16(out + 16, bsize - 1);
    put32(out + HEADER + clen, crc);
    put32(out + HEADER + clen + 4, len);
    return bsize;
  }
};

// inflate one BGZF block (block, bsize) into out (isize bytes), checking its
// CRC32
struct Decompressor {
  bool use_ld;
  z_stream zs;
  libdeflate_decompressor* ld;
  Decompressor(bool libdeflate) : use_ld(libdeflate), ld(NULL) {
    if (use_ld) {
      ld = libdeflate_alloc_decompressor();
    } else {
      std::memset(&zs, 0, sizeof(zs));
      if (inflateInit2(&zs, -15) != Z_OK) fail("inflateInit2");
    }
  }
  ~Decompressor() {
    if (use_ld) libdeflate_free_decompressor(ld); else inflateEnd(&zs);
  }
  void inflate_block(const unsigned char* block, size_t bsize,
                     unsigned char* out, size_t isize) {
    const unsigned char* data = block + HEADER;
    size_t dlen = bsize - HEADER - FOOTER;
    uint32_t crc;
    if (use_ld) {
      size_t actual;
      if (libdeflate_deflate_decompress(ld, data, dlen, out, isize, &actual) !=
          LIBDEFLATE_SUCCESS || actual != isize) fail("inflate");
      crc = libdeflate_crc32(0, out, isize);
    } else {
      inflateReset(&zs);
      zs.next_in = (Bytef*) data;
      zs.avail_in = dlen;
      zs.next_out = out;
      zs.avail_out = isize;
      if (inflate(&zs, Z_FINISH) != Z_STREAM_END || zs.avail_out != 0)
        fail("inflate");
      crc = crc32(crc32(0L, Z_NULL, 0), out, isize);
    }
    if (crc != get32(block + bsize - 8)) fail("CRC mismatch");
  }
};

struct BlockTable {
  std::vector<size_t> offset, bsize, ioffset, isize;
};

// the blocks of a BGZF file in memory (without the EOF block)
static BlockTable parse_blocks(const std::vector<unsigned char> &f) {
  BlockTable t;
  size_t pos = 0, ipos = 0;
  while (pos + HEADER <= f.size()) {
    const unsigned char* b = f.data() + pos;
    if (b[0] != 31 || b[1] != 139 || b[12] != 'B' || b[13] != 'C')
      fail("not BGZF");
    size_t bsize = (b[16] | (b[17] << 8)) + 1;
    size_t isize = get32(b + bsize - 4);
    if (isize > 0) {
      t.offset.push_back(pos);
      t.bsize.push_back(bsize);
      t.ioffset.push_back(ipos);
      t.isize.push_back(isize);
    }
    pos += bsize;
    ipos += isize;
  }
  return t;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <file.tsv> <out_dir> [threads] [reps] "
                 "[lookups]\n", argv[0]);
    return 1;
  }
  std::string in_path = argv[1];
  std::string out_dir = argv[2];
  int threads = argc > 3 ? std::atoi(argv[3]) : 4;
  int reps = argc > 4 ? std::atoi(argv[4]) : 3;
  int lookups = argc > 5 ? std::atoi(argv[5]) : 2000;

  std::vector<unsigned char> input = read_file(in_path);
  size_t n = input.size();
  size_t n_blocks = (n + BLOCK_IN - 1) / BLOCK_IN;
  std::mt19937_64 rng(42);

  std::printf("codec\tthreads\tin_bytes\tout_bytes\tratio\tcomp_1t_s\t"
              "comp_nt_s\tdecomp_1t_s\tdecomp_nt_s\tdecomp_scan_1t_s\t"
              "lookup_us\tlookup64k_us\tn_blocks\tcheck\n");

  // the plain file: scan its lines, and random reads of 1 MB and 64 KB
  {
    std::vector<double> t_scan;
    uint64_t sum = 0, n_lines = 0;
    for (int r = 0; r < reps; r++) {
      n_lines = 0;
      double t0 = now();
      sum = scan_lines(input.data(), n, n_lines);
      t_scan.push_back(now() - t0);
    }
    int fd = open(in_path.c_str(), O_RDONLY);
    std::vector<unsigned char> buf(1 << 20);
    std::uniform_int_distribution<size_t> pick(0, n > buf.size() ?
                                               n - buf.size() : 0);
    double t_1m = 0, t_64k = 0;
    for (int r = 0; r < 2; r++) {
      // (the first pass warms up)
      double t0 = now();
      for (int i = 0; i < lookups; i++) {
        if (pread(fd, buf.data(), 1 << 20, pick(rng)) < 0) fail("pread");
      }
      t_1m = (now() - t0) / lookups;
      t0 = now();
      for (int i = 0; i < lookups; i++) {
        if (pread(fd, buf.data(), 1 << 16, pick(rng)) < 0) fail("pread");
      }
      t_64k = (now() - t0) / lookups;
    }
    close(fd);
    std::printf("plain\t%d\t%zu\t%zu\t1\tNA\tNA\tNA\tNA\t%.4f\t%.1f\t%.1f\t"
                "NA\t%llu lines, key sum %llu\n", threads, n, n,
                median(t_scan), t_1m * 1e6, t_64k * 1e6,
                (unsigned long long) n_lines, (unsigned long long) sum);
    std::fflush(stdout);
  }

  std::vector<Codec> codecs = {
    {"zlib-1", false, 1}, {"zlib-6", false, 6},
    {"libdeflate-1", true, 1}, {"libdeflate-6", true, 6}};
  // one output slot of BLOCK_MAX bytes per block
  std::vector<unsigned char> slots(n_blocks * BLOCK_MAX);
  std::vector<size_t> sizes(n_blocks);
  std::vector<unsigned char> output(n);

  for (const Codec &codec : codecs) {
    // compression, 1 thread and `threads` threads
    auto compress_all = [&](int nt) {
      double t0 = now();
      #pragma omp parallel num_threads(nt)
      {
        Compressor c(codec);
        #pragma omp for schedule(dynamic, 8)
        for (size_t b = 0; b < n_blocks; b++) {
          size_t len = std::min(BLOCK_IN, n - b * BLOCK_IN);
          sizes[b] = c.compress(input.data() + b * BLOCK_IN, len,
                                slots.data() + b * BLOCK_MAX);
        }
      }
      return now() - t0;
    };
    std::vector<double> t_c1, t_cn;
    for (int r = 0; r < reps; r++) t_c1.push_back(compress_all(1));
    for (int r = 0; r < reps; r++) t_cn.push_back(compress_all(threads));

    // write the BGZF file (with the EOF block)
    std::string out_path = out_dir + "/sidecar." + codec.name + ".gz";
    FILE* f = std::fopen(out_path.c_str(), "wb");
    if (!f) fail("cannot create " + out_path);
    size_t out_bytes = 0;
    for (size_t b = 0; b < n_blocks; b++) {
      std::fwrite(slots.data() + b * BLOCK_MAX, 1, sizes[b], f);
      out_bytes += sizes[b];
    }
    std::fwrite(BGZF_EOF, 1, 28, f);
    out_bytes += 28;
    if (std::fclose(f) != 0) fail("write " + out_path);

    // full decompression, from the file in memory
    std::vector<unsigned char> comp = read_file(out_path);
    BlockTable t = parse_blocks(comp);
    if (t.offset.size() != n_blocks) fail("block count");
    // (the checksum of the scan is stored, so that the compiler can not drop
    // the scan from the timed path)
    volatile uint64_t scan_sink = 0;
    auto decompress_all = [&](int nt, bool scan) {
      uint64_t n_lines = 0;
      double t0 = now();
      #pragma omp parallel num_threads(nt)
      {
        Decompressor d(codec.libdeflate);
        #pragma omp for schedule(dynamic, 8)
        for (size_t b = 0; b < t.offset.size(); b++) {
          d.inflate_block(comp.data() + t.offset[b], t.bsize[b],
                          output.data() + t.ioffset[b], t.isize[b]);
        }
      }
      if (scan) scan_sink = scan_sink + scan_lines(output.data(), n, n_lines);
      return now() - t0;
    };
    std::vector<double> t_d1, t_dn, t_ds;
    std::fill(output.begin(), output.end(), 0);
    for (int r = 0; r < reps; r++) t_d1.push_back(decompress_all(1, false));
    bool same = std::memcmp(output.data(), input.data(), n) == 0;
    for (int r = 0; r < reps; r++) t_dn.push_back(decompress_all(threads,
                                                                 false));
    for (int r = 0; r < reps; r++) t_ds.push_back(decompress_all(1, true));

    // random lookups: pread() one block of the file and inflate it
    int fd = open(out_path.c_str(), O_RDONLY);
    std::vector<unsigned char> blk(BLOCK_MAX), out(BLOCK_IN);
    std::uniform_int_distribution<size_t> pick(0, t.offset.size() - 1);
    Decompressor d(codec.libdeflate);
    double t_lookup = 0;
    for (int r = 0; r < 2; r++) {
      double t0 = now();
      for (int i = 0; i < lookups; i++) {
        size_t b = pick(rng);
        if (pread(fd, blk.data(), t.bsize[b], t.offset[b]) !=
            (ssize_t) t.bsize[b]) fail("pread");
        d.inflate_block(blk.data(), t.bsize[b], out.data(), t.isize[b]);
      }
      t_lookup = (now() - t0) / lookups;
    }
    close(fd);

    std::printf("%s\t%d\t%zu\t%zu\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\t"
                "NA\t%zu\t%s\n", codec.name.c_str(), threads, n, out_bytes,
                (double) n / out_bytes, median(t_c1), median(t_cn),
                median(t_d1), median(t_dn), median(t_ds), t_lookup * 1e6,
                n_blocks, same ? "identical" : "DIFFERENT");
    std::fflush(stdout);
  }
  return 0;
}
