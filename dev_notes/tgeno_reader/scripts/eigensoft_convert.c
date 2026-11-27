// convert a GENO file to TGENO (or back) with EIGENSOFT's own geno.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geno.h"
int genotype_matrix_read_file(genotype_matrix *this, const char * const genotype_filename, int flags, const unsigned int *individual_file_hash, const unsigned int *snp_file_hash);
int genotype_matrix_write_packed_file(const genotype_matrix *this, const char *filename, bool transpose);
int main(int argc, char **argv) {
  if (argc != 4) { fprintf(stderr, "usage: convert in out T|G\n"); return 1; }
  genotype_matrix g;
  int r = genotype_matrix_read_file(&g, argv[1], 0, NULL, NULL);
  if (r != 0) { fprintf(stderr, "read failed %d\n", r); return 1; }
  genotype_matrix_write_packed_file(&g, argv[2], argv[3][0] == 'T');
  return 0;
}
