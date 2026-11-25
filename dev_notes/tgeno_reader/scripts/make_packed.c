// make_packed <prefix> <n> <m>: synthetic packedancestrymap files with the
// same genotypes in GENO (<prefix>.geno) and TGENO (<prefix>.tgeno) layouts,
// plus <prefix>.snp and <prefix>.ind. Genotype (i, j) is a hash of (i, j):
// 10% missing (3), else 0, 1 or 2.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static inline uint64_t mix(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
static inline unsigned g(uint64_t i, uint64_t j, uint64_t m) {
  uint64_t h = mix(i * m + j);
  if ((h & 0xff) < 26) return 3;
  return (unsigned)((h >> 8) % 3);
}
int main(int argc, char **argv) {
  const char *p = argv[1];
  uint64_t n = strtoull(argv[2], 0, 10), m = strtoull(argv[3], 0, 10);
  char fn[4096];
  // GENO
  size_t rec = (n + 3) / 4; if (rec < 48) rec = 48;
  unsigned char *buf = calloc(rec > 48 ? rec : 48, 1);
  snprintf(fn, sizeof fn, "%s.geno", p);
  FILE *f = fopen(fn, "wb");
  memset(buf, 0, rec);
  sprintf((char *)buf, "GENO %7llu %7llu %x %x", (unsigned long long)n, (unsigned long long)m, 0, 0);
  fwrite(buf, 1, rec, f);
  for (uint64_t j = 0; j < m; j++) {
    memset(buf, 0, rec);
    for (uint64_t i = 0; i < n; i++) buf[i / 4] |= g(i, j, m) << (6 - 2 * (i % 4));
    fwrite(buf, 1, rec, f);
  }
  fclose(f); free(buf);
  // TGENO
  size_t trec = (m + 3) / 4;
  buf = calloc(trec > 48 ? trec : 48, 1);
  snprintf(fn, sizeof fn, "%s.tgeno", p);
  f = fopen(fn, "wb");
  sprintf((char *)buf, "TGENO %7llu %7llu %x %x", (unsigned long long)n, (unsigned long long)m, 0, 0);
  fwrite(buf, 1, 48, f);
  for (uint64_t i = 0; i < n; i++) {
    memset(buf, 0, trec);
    for (uint64_t j = 0; j < m; j++) buf[j / 4] |= g(i, j, m) << (6 - 2 * (j % 4));
    fwrite(buf, 1, trec, f);
  }
  fclose(f);
  // .snp (22 chromosomes) and .ind
  snprintf(fn, sizeof fn, "%s.snp", p);
  f = fopen(fn, "w");
  uint64_t per_chr = m / 22 + 1;
  for (uint64_t j = 0; j < m; j++)
    fprintf(f, "rs%llu\t%llu\t0.0\t%llu\tA\tG\n", (unsigned long long)j + 1,
            (unsigned long long)(j / per_chr + 1), (unsigned long long)(j % per_chr) * 100 + 1);
  fclose(f);
  snprintf(fn, sizeof fn, "%s.ind", p);
  f = fopen(fn, "w");
  for (uint64_t i = 0; i < n; i++) fprintf(f, "I%llu U Pop%llu\n", (unsigned long long)i + 1, (unsigned long long)(i % 50));
  fclose(f);
  return 0;
}
