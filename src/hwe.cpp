#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

#include <bigstatsr/BMCodeAcc.h>


// The functions below are used for HWE. They come from PLINK 1.90 under GPL3, based on
// #3848a39 on https://github.com/chrchang/plink-ng/1.9/plink_stats.c
// Full credit goes to the original authors of these functions (see below for details)

// This file is part of PLINK 1.90, copyright (C) 2005-2023 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#define SMALL_EPSILON 0.00000000000005684341886080801486968994140625
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 000 - 000.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // - Proper handling of >64k genotypes.  Previously, there was a potential
  //   integer overflow.
  // - Detection and efficient handling of floating point overflow and
  //   underflow.  E.g. instead of summing a tail all the way down, the loop
  //   stops once the latest increment underflows the partial sum's 53-bit
  //   precision; this results in a large speedup when max heterozygote count
  //   >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  //   sums.
  // - Support for the mid-p variant of this test.  See Graffelman J, Moreno V
  //   (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium.
  //
  // Note that the SNPHWE_t() function below is a lot more efficient for
  // testing against a p-value inclusion threshold.  SNPHWE2() should only be
  // used if you need the actual p-value.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  int32_t tie_ct = 1;
  double curr_hets_t2 = obs_hets;
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;
  double preaddp;
  if (!genotypes2) {
    if (midp) {
      return 0.5;
    } else {
      return 1;
    }
  }

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    while (curr_hets_t2 > 1.5) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
        if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
          tie_ct++;
        }
        tailp += lastp2;
        break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
        return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_hets_t2 > 1.5) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
        break;
      }
    }
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_homr_t1 > 0.5) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
        break;
      }
      curr_hets_t1 += 2;
      curr_homr_t1 -= 1;
      curr_homc_t1 -= 1;
    }
  } else {
    // tail 1 = lower
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
        if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
          tie_ct++;
        }
        tailp += lastp2;
        break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
        return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
        break;
      }
    }
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_hets_t1 > 1.5) {
      curr_homr_t1 += 1;
      curr_homc_t1 += 1;
      lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
        break;
      }
      curr_hets_t1 -= 2;
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  } else {
    return (tailp - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (tailp + centerp);
  }
}

// [[Rcpp::export]]
double SNPHWE2_R(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp) {
  // obs_hets: number of heterozygotes
  // obs_hom1: number of homozygotes for allele 1
  // obs_hom2: number of homozygotes for allele 2
  // midp: if true, use mid-p method
  return SNPHWE2(obs_hets, obs_hom1, obs_hom2, midp);
}

// Function to apply the hwe on each column of a matrix generated by big_counts
// [[Rcpp::export]]
NumericVector hwe_on_matrix(IntegerMatrix geno_counts, uint32_t midp) {
  // geno_counts: a matrix of genotype counts, with row for obs_hom1, obs_hets,  obs_hom2
  // midp: if true, use mid-p method
  int n = geno_counts.ncol();
  NumericVector p_values(n);
  for (int i = 0; i < n; i++) {
    p_values[i] = SNPHWE2_R(geno_counts(1, i), geno_counts(0, i), geno_counts(2, i), midp);
  }
  return p_values;
}

/******************************************************************************/
// Estimate HWE p on a grouped gen_tibble
// this only works for diploid data


// [[Rcpp::export]]
NumericMatrix gt_grouped_hwe(Environment BM,
                                     const IntegerVector& rowInd,
                                     const IntegerVector& colInd,
                                     const IntegerVector& groupIds,
                                     size_t ngroups,
                                     uint32_t midp) {
  // midp: if true, use mid-p method
  
  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);
  
  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci
  
  NumericMatrix hwe_p(m, ngroups);
  IntegerMatrix genotypes(3, ngroups);
  
  for (size_t j = 0; j < m; j++) { // for each locus
    genotypes.fill(0); // set all to zero
    for (size_t i = 0; i < n; i++) { // for each individual
      double x = macc(i, j);
      if (x > -1){
        genotypes(x, groupIds[i]) += 1;
      }
    }
    // now for each group, estimate HWE p
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      hwe_p(j, group_i) = SNPHWE2(genotypes(1,group_i), genotypes(0,group_i), genotypes(2,group_i), midp);
    }
  }
  
  return hwe_p;
}

/******************************************************************************/

