/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/
// Function to copy a matrix into a BM, in parallel

// [[Rcpp::export]]
void write_to_FBM (Environment BM,
                   IntegerMatrix& allele_counts,
                   const int col_start,
                   const int n_loci,
                   int ncores) {
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  BMAcc_RW<unsigned char> macc_fbm(xpBM);

  size_t n_indiv = macc_fbm.nrow();

  #pragma omp parallel for num_threads(ncores)
  for (int locus_i = 0; locus_i < n_loci; locus_i++){
    for (size_t indiv_i = 0; indiv_i < n_indiv; indiv_i++)
      macc_fbm(indiv_i,col_start+locus_i) = allele_counts(indiv_i, locus_i);
  }

}
