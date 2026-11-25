// Microbenchmark of the TGENO -> FBM transposition kernels, outside R.
// g++ -O2 -fopenmp -o transpose_micro transpose_micro.cpp
// ./transpose_micro <file.tgeno> <threads> <kernel: 1 byte at a time, 4 or 5 words>
//   <out: 0 pre-faulted memory, 1 fresh mmap file (as the FBM)> <IB> <BB> <mode: 0 full,
//   1 no stores, 2 no loads> <prefetch distance (individuals)> <stage: 0/1> <huge input: 0/1>
// The dimensions are fixed to those of make_packed aadr_like 4866 1233013.
// microbenchmark of TGENO -> FBM transpose kernels
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>
#include <sys/mman.h>
static int STAGE = 0, HUGE_IN = 0;
static void* huge_alloc(size_t sz, bool huge) { void* p = 0; if (posix_memalign(&p, 1 << 21, sz)) return 0; if (huge) madvise(p, sz, MADV_HUGEPAGE); return p; }

static unsigned char lut[256][4];
static unsigned char map4[4] = {2, 1, 0, 3};
static size_t IB = 64, BB = 64; static int MODE = 0; static size_t PD = 0; static volatile uint64_t SINK = 0;

// K1: current kernel
static void k1(const unsigned char* slab, size_t sb, size_t fb, size_t b0, size_t b1, size_t n, size_t m, unsigned char* fbm) {
  const size_t full = m / 4;
  const size_t b1f = std::min(b1, full > fb ? full - fb : 0);
  for (size_t bb0 = b0; bb0 < b1f; bb0 += 64) {
    size_t bb1 = std::min(b1f, bb0 + 64);
    for (size_t i0 = 0; i0 < n; i0 += 64) {
      size_t i1 = std::min(n, i0 + 64);
      for (size_t b = bb0; b < bb1; b++) {
        const unsigned char* src = slab + b;
        unsigned char* d0 = fbm + 4 * (fb + b) * n; unsigned char* d1 = d0 + n; unsigned char* d2 = d1 + n; unsigned char* d3 = d2 + n;
        for (size_t i = i0; i < i1; i++) { const unsigned char* c = lut[src[i * sb]]; d0[i] = c[0]; d1[i] = c[1]; d2[i] = c[2]; d3[i] = c[3]; }
      }
    }
  }
}

// SWAR: map 2-bit fields, 8x8 byte transpose, expand
static inline uint64_t map_fields(uint64_t x, const uint64_t* mv) {
  const uint64_t F = 0x5555555555555555ULL;
  uint64_t lo = x & F, hi = (x >> 1) & F;
  return (~(lo | hi) & F) * mv[0] + (lo & ~hi) * mv[1] + (hi & ~lo) * mv[2] + (lo & hi) * mv[3];
}
static inline void tr8(uint64_t* a) {
  for (int r = 0; r < 4; r++) { uint64_t t = ((a[r] >> 32) ^ a[r + 4]) & 0x00000000FFFFFFFFULL; a[r] ^= t << 32; a[r + 4] ^= t; }
  for (int r : {0, 1, 4, 5}) { uint64_t t = ((a[r] >> 16) ^ a[r + 2]) & 0x0000FFFF0000FFFFULL; a[r] ^= t << 16; a[r + 2] ^= t; }
  for (int r : {0, 2, 4, 6}) { uint64_t t = ((a[r] >> 8) ^ a[r + 1]) & 0x00FF00FF00FF00FFULL; a[r] ^= t << 8; a[r + 1] ^= t; }
}
static void k4(const unsigned char* slab, size_t sb, size_t fb, size_t b0, size_t b1, size_t n, size_t m, unsigned char* fbm) {
  uint64_t mv[4]; for (int v = 0; v < 4; v++) mv[v] = map4[v];
  const size_t full = m / 4;
  const size_t b1f = std::min(b1, full > fb ? full - fb : 0);
  const size_t n8 = n & ~size_t(7);
  const uint64_t M3 = 0x0303030303030303ULL;
  for (size_t bb0 = b0; bb0 < b1f; bb0 += BB) {
    size_t bb1 = std::min(b1f, bb0 + BB);
    size_t bb8 = bb0 + ((bb1 - bb0) & ~size_t(7));
    for (size_t i0 = 0; i0 < n8; i0 += IB) {
      size_t i1 = std::min(n8, i0 + IB);
      for (size_t i = i0; i < i1; i += 8) {
        for (size_t b = bb0; b < bb8; b += 8) {
          uint64_t a[8];
          for (int r = 0; r < 8; r++) { std::memcpy(&a[r], slab + (i + r) * sb + b, 8); a[r] = map_fields(a[r], mv); }
          tr8(a);
          for (int c = 0; c < 8; c++) {
            unsigned char* d = fbm + 4 * (fb + b + c) * n + i;
            uint64_t v0 = (a[c] >> 6) & M3, v1 = (a[c] >> 4) & M3, v2 = (a[c] >> 2) & M3, v3 = a[c] & M3;
            std::memcpy(d, &v0, 8); std::memcpy(d + n, &v1, 8); std::memcpy(d + 2 * n, &v2, 8); std::memcpy(d + 3 * n, &v3, 8);
          }
        }
      }
    }
    // leftovers: individuals n8..n for all bytes; bytes bb8..bb1 for all individuals
    for (size_t b = bb0; b < bb1; b++) {
      size_t ia = (b < bb8) ? n8 : 0;
      unsigned char* d0 = fbm + 4 * (fb + b) * n;
      for (size_t i = ia; i < n; i++) { const unsigned char* c = lut[slab[i * sb + b]]; d0[i] = c[0]; d0[n + i] = c[1]; d0[2 * n + i] = c[2]; d0[3 * n + i] = c[3]; }
    }
  }
}


// K5: tiles of IB individuals x BB bytes; within a tile, for each group of 8
// bytes (32 output rows), sweep the individuals 8 at a time
static void k5(const unsigned char* slab, size_t sb, size_t fb, size_t b0, size_t b1, size_t n, size_t m, unsigned char* fbm) {
  uint64_t acc = 0;
  uint64_t mv[4]; for (int v = 0; v < 4; v++) mv[v] = map4[v];
  const size_t full = m / 4;
  const size_t b1f = std::min(b1, full > fb ? full - fb : 0);
  const size_t n8 = n & ~size_t(7);
  const uint64_t M3 = 0x0303030303030303ULL;
  for (size_t bb0 = b0; bb0 < b1f; bb0 += BB) {
    size_t bb1 = std::min(b1f, bb0 + BB);
    size_t bb8 = bb0 + ((bb1 - bb0) & ~size_t(7));
    for (size_t i0 = 0; i0 < n8; i0 += IB) {
      size_t i1 = std::min(n8, i0 + IB);
      for (size_t b = bb0; b < bb8; b += 8) {
        unsigned char* dcol = fbm + 4 * (fb + b) * n;
        for (size_t i = i0; i < i1; i += 8) {
          uint64_t a[8];
          const unsigned char* s = slab + i * sb + b;
          if (PD && i + PD + 8 <= n) { const unsigned char* p = s + PD * sb; for (int r = 0; r < 8; r++) __builtin_prefetch(p + r * sb); }
          if (MODE == 2) { for (int r = 0; r < 8; r++) a[r] = i * 0x0101010101010101ULL + r; }
          else { for (int r = 0; r < 8; r++) { std::memcpy(&a[r], s + r * sb, 8); a[r] = map_fields(a[r], mv); } }
          tr8(a);
          unsigned char* d = dcol + i;
          if (MODE == 1) { uint64_t x = 0; for (int c = 0; c < 8; c++) x ^= a[c]; acc ^= x; continue; }
          for (int c = 0; c < 8; c++) {
            uint64_t v0 = (a[c] >> 6) & M3, v1 = (a[c] >> 4) & M3, v2 = (a[c] >> 2) & M3, v3 = a[c] & M3;
            std::memcpy(d, &v0, 8); std::memcpy(d + n, &v1, 8); std::memcpy(d + 2 * n, &v2, 8); std::memcpy(d + 3 * n, &v3, 8);
            d += 4 * n;
          }
        }
      }
    }
    for (size_t b = bb0; b < bb1; b++) {
      size_t ia = (b < bb8) ? n8 : 0;
      unsigned char* d0 = fbm + 4 * (fb + b) * n;
      for (size_t i = ia; i < n; i++) { const unsigned char* c = lut[slab[i * sb + b]]; d0[i] = c[0]; d0[n + i] = c[1]; d0[2 * n + i] = c[2]; d0[3 * n + i] = c[3]; }
    }
  }
  SINK ^= acc;
}

int main(int argc, char** argv) {
  const char* file = argv[1]; int threads = atoi(argv[2]); int kern = atoi(argv[3]); int out_mode = atoi(argv[4]); // 0 prefaulted malloc, 1 fresh mmap file
  const char* outfile = "out.bk"; if (argc > 6) { IB = atoi(argv[5]); BB = atoi(argv[6]); } if (argc > 7) MODE = atoi(argv[7]); if (argc > 8) PD = atoi(argv[8]); if (argc > 9) STAGE = atoi(argv[9]); if (argc > 10) HUGE_IN = atoi(argv[10]);
  size_t n = 4866, m = 1233013, rec = (m + 3) / 4;
  for (int b = 0; b < 256; b++) for (int k = 0; k < 4; k++) lut[b][k] = map4[(b >> (6 - 2 * k)) & 3];
  unsigned char* inp = (unsigned char*)huge_alloc(n * rec, HUGE_IN);
  FILE* f = fopen(file, "rb"); fseek(f, 48, SEEK_SET); if (fread(inp, 1, n * rec, f) != n * rec) return 1; fclose(f);
  struct { unsigned char* data() { return p; } unsigned char* p; } in = {inp};
  unsigned char* out; size_t osz = n * m;
  int fd = -1;
  if (out_mode == 0) { out = (unsigned char*)malloc(osz); memset(out, 0, osz); }
  else { fd = open(outfile, O_RDWR | O_CREAT | O_TRUNC, 0644); if (ftruncate(fd, osz)) return 2; out = (unsigned char*)mmap(0, osz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0); }
  size_t task = std::max<size_t>(BB, (((1u << 20) / (4 * n)) + BB - 1) / BB * BB);
  size_t ntask = (rec + task - 1) / task;
  auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(threads)
  {
    unsigned char* stage = STAGE ? (unsigned char*)huge_alloc(4 * task * n + (1 << 21), true) : 0;
#pragma omp for schedule(dynamic, 1)
  for (long t = 0; t < (long)ntask; t++) {
    size_t b0 = t * task, b1 = std::min(rec, b0 + task);
    if (STAGE) {
      // transpose into the staging buffer as if it were the FBM starting at locus 4*b0 (fb = b0, base shifted)
      size_t bfull = std::min(b1, m / 4);
      if (bfull > b0) {
        unsigned char* base = stage - 4 * b0 * n;  // so that fbm + 4*(fb+b)*n lands in stage
        k5(in.data() + b0, rec, b0, 0, bfull - b0, n, m, base);
        std::memcpy(out + 4 * b0 * n, stage, 4 * (bfull - b0) * n);
      }
    } else {
    if (kern == 1) k1(in.data(), rec, 0, b0, b1, n, m, out); else if (kern == 4) k4(in.data(), rec, 0, b0, b1, n, m, out); else k5(in.data(), rec, 0, b0, b1, n, m, out);
    }
  }
    free(stage);
  }
  // last partial byte
  if (m % 4) { size_t b = m / 4; for (size_t i = 0; i < n; i++) for (size_t j = 4 * b; j < m; j++) out[j * n + i] = lut[inp[i * rec + b]][j - 4 * b]; }
  auto t1 = std::chrono::steady_clock::now();
  unsigned long long cs = 0; for (size_t x = 0; x < osz; x += 4093) cs = cs * 31 + out[x];
  printf("kernel %d threads %d out_mode %d: %.3f s  checksum %llu\n", kern, threads, out_mode, std::chrono::duration<double>(t1 - t0).count(), cs);
  if (fd >= 0) { munmap(out, osz); close(fd); unlink(outfile); }
  return 0;
}
