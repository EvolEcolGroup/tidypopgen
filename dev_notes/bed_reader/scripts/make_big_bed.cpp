// Make a synthetic PLINK 1 binary fileset (.bed, .bim, .fam) for timing the
// readers of .bed files.
// usage: make_big_bed <out prefix> <n individuals> <n snps> [seed]
// Genotypes: for each snp, an allele frequency p ~ U(0.02, 0.98), genotypes in
// Hardy-Weinberg proportions, and 1% missing; each snp has its own random
// stream (seeded from the seed and its index), so the files do not depend on
// the number of threads. Loci are on chromosomes 1-22 (in proportion to
// their index), with positions 50-249 bp apart, names rs<index> and two
// different random alleles of A, C, G, T. Individuals are in 20 populations.
// build: g++ -O2 -fopenmp make_big_bed.cpp -o make_big_bed
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static inline double unif(uint64_t &x) {
  return (splitmix64(x) >> 11) * (1.0 / 9007199254740992.0);
}

int main(int argc, char **argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s <out prefix> <n indiv> <n snps> [seed]\n",
                 argv[0]);
    return 1;
  }
  std::string prefix = argv[1];
  std::size_t n = std::strtoull(argv[2], NULL, 10);
  std::size_t m = std::strtoull(argv[3], NULL, 10);
  uint64_t seed = (argc > 4) ? std::strtoull(argv[4], NULL, 10) : 20261120;
  std::size_t bytes_per_snp = (n + 3) / 4;

  // .fam
  FILE *fam = std::fopen((prefix + ".fam").c_str(), "w");
  for (std::size_t i = 0; i < n; i++) {
    std::fprintf(fam, "pop%zu ind%zu 0 0 %d -9\n", i % 20 + 1, i + 1,
                 (int) (i % 2) + 1);
  }
  std::fclose(fam);

  // .bim
  FILE *bim = std::fopen((prefix + ".bim").c_str(), "w");
  const char nuc[4] = {'A', 'C', 'G', 'T'};
  uint64_t rs = seed ^ 0xABCDEF;
  long pos = 0;
  int last_chr = 0;
  for (std::size_t j = 0; j < m; j++) {
    int chr = 1 + (int) ((j * 22) / m);
    if (chr != last_chr) {
      pos = 0;
      last_chr = chr;
    }
    pos += 50 + (long) (splitmix64(rs) % 200);
    int a1 = (int) (splitmix64(rs) % 4);
    int a2 = (a1 + 1 + (int) (splitmix64(rs) % 3)) % 4;
    std::fprintf(bim, "%d\trs%zu\t0\t%ld\t%c\t%c\n", chr, j + 1, pos,
                 nuc[a1], nuc[a2]);
  }
  std::fclose(bim);

  // .bed, written in chunks of snps generated in parallel
  FILE *bed = std::fopen((prefix + ".bed").c_str(), "wb");
  const unsigned char magic[3] = {0x6C, 0x1B, 0x01};
  std::fwrite(magic, 1, 3, bed);
  std::size_t chunk = std::max<std::size_t>(1, (256u << 20) / bytes_per_snp);
  std::vector<unsigned char> buf(chunk * bytes_per_snp);
  for (std::size_t j0 = 0; j0 < m; j0 += chunk) {
    std::size_t nj = std::min(chunk, m - j0);
#pragma omp parallel for schedule(static)
    for (long jj = 0; jj < (long) nj; jj++) {
      uint64_t s = seed ^ ((j0 + jj + 1) * 0xD1B54A32D192ED03ULL);
      double p = 0.02 + 0.96 * unif(s);
      double p00 = p * p, p10 = p00 + 2 * p * (1 - p);
      unsigned char *out = &buf[jj * bytes_per_snp];
      std::fill(out, out + bytes_per_snp, 0);
      for (std::size_t i = 0; i < n; i++) {
        double u = unif(s);
        unsigned code;
        if (u < 0.01) {
          code = 1;  // missing
        } else {
          double v = (u - 0.01) / 0.99;
          code = (v < p00) ? 0 : ((v < p10) ? 2 : 3);
        }
        out[i / 4] |= (unsigned char) (code << (2 * (i % 4)));
      }
    }
    std::fwrite(buf.data(), 1, nj * bytes_per_snp, bed);
  }
  std::fclose(bed);
  return 0;
}
