# Package index

## `gen_tibble`

Functions for creating, saving, and loading `gen_tibble` objects.

- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  :

  Constructor for a `gen_tibble`

- [`gt_save()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_save.md)
  : Save a gen_tibble

- [`gt_load()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_load.md)
  : Load a gen_tibble

- [`gt_get_file_names()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_get_file_names.md)
  :

  Get the names of files storing the genotypes of a `gen_tibble`

- [`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_order_loci.md)
  : Order the loci table of a gen_tibble

- [`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_update_backingfile.md)
  : Update the backing matrix

- [`gt_add_sf()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_add_sf.md)
  :

  Add an simple feature geometry to a `gen_tibble`

- [`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pseudohaploid.md)
  :

  Set the ploidy of a `gen_tibble` to include pseudohaploids

## Show

Functions to view loci and genotype information stored in a
`gen_tibble`.

- [`show_genotypes()`](https://evolecolgroup.github.io/tidypopgen/reference/show_genotypes.md)
  :

  Show the genotypes of a `gen_tibble`

- [`show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md)
  [`` `show_loci<-`() ``](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md)
  :

  Show the loci information of a `gen_tibble`

- [`show_ploidy()`](https://evolecolgroup.github.io/tidypopgen/reference/show_ploidy.md)
  :

  Show the ploidy information of a `gen_tibble`

- [`count_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/count_loci.md)
  :

  Count the number of loci in a `gen_tibble`

- [`augment_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/augment_loci.md)
  : Augment the loci table with information from a analysis object

- [`is_loci_table_ordered()`](https://evolecolgroup.github.io/tidypopgen/reference/is_loci_table_ordered.md)
  : Test if the loci table is ordered

- [`find_duplicated_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/find_duplicated_loci.md)
  : Find duplicates in the loci table

## Selecting loci

Functions to select loci from a `gen_tibble`, either based on loci
information, or genotype information.

- [`select_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/select_loci.md)
  :

  The `select` verb for `loci`

- [`select_loci_if()`](https://evolecolgroup.github.io/tidypopgen/reference/select_loci_if.md)
  :

  The `select_if` verb for `loci`

## Imputing

Functions for imputing missing data in a `gen_tibble`.

- [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_impute_simple.md)
  : Simple imputation based on allele frequencies

- [`gt_impute_xgboost()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_impute_xgboost.md)
  : Imputation based XGBoost

- [`gt_has_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_has_imputed.md)
  :

  Checks if a `gen_tibble` has been imputed

- [`gt_uses_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_uses_imputed.md)
  :

  Checks if a `gen_tibble` uses imputed data

- [`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_set_imputed.md)
  :

  Sets a `gen_tibble` to use imputed data

## Merging

Functions for merging `gen_tibble` objects.

- [`rbind_dry_run()`](https://evolecolgroup.github.io/tidypopgen/reference/rbind_dry_run.md)
  : Generate a report of what would happen to each SNP in a merge
- [`rbind(`*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/rbind.gen_tbl.md)
  : Combine two gen_tibbles
- [`summary(`*`<rbind_report>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/summary_rbind_dry_run.md)
  : Print a summary of a merge report
- [`cbind(`*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/cbind.gen_tbl.md)
  : Combine a gen_tibble to a data.frame or tibble by column

## Exporting

Functions for exporting a `gen_tibble` to another format.

- [`gt_as_genind()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_genind.md)
  :

  Convert a `gen_tibble` to a `genind` object from `adegenet`

- [`gt_as_genlight()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_genlight.md)
  :

  Convert a `gen_tibble` to a `genlight` object from `adegenet`

- [`gt_as_geno_lea()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_geno_lea.md)
  :

  Convert a `gentibble` to a .geno file for sNMF from the LEA package

- [`gt_as_hierfstat()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_hierfstat.md)
  :

  Convert a `gen_tibble` to a data.frame compatible with `hierfstat`

- [`gt_as_plink()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_plink.md)
  :

  Export a `gen_tibble` object to PLINK bed format

- [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_as_vcf.md)
  :

  Convert a `gen_tibble` to a VCF

## Importing

Functions for importing another format into a `gen_tibble`.

- [`gt_from_genlight()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_from_genlight.md)
  :

  Convert a `genlight` object from adegenet to a `gen_tibble`

## Quality control

Functions reporting summary statistics for a `gen_tibble`.

- [`qc_report_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/qc_report_loci.md)
  : Create a Quality Control report for loci

- [`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/reference/qc_report_indiv.md)
  : Create a Quality Control report for individuals

- [`autoplot(`*`<qc_report_loci>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot.qc_report_loci.md)
  :

  Autoplots for `qc_report_loci` objects

- [`autoplot(`*`<qc_report_indiv>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot.qc_report_indiv.md)
  :

  Autoplots for `qc_report_indiv` objects

- [`filter_high_relatedness()`](https://evolecolgroup.github.io/tidypopgen/reference/filter_high_relatedness.md)
  : Filter individuals based on a relationship threshold

## Loci

Functions operating across loci.

- [`loci_alt_freq()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_alt_freq.md)
  [`loci_maf()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_alt_freq.md)
  : Estimate allele frequencies at each locus

- [`loci_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_chromosomes.md)
  :

  Get the chromosomes of loci in a `gen_tibble`

- [`loci_hwe()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_hwe.md)
  : Test Hardy-Weinberg equilibrium at each locus

- [`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_ld_clump.md)
  : Clump loci based on a Linkage Disequilibrium threshold

- [`loci_missingness()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_missingness.md)
  : Estimate missingness at each locus

- [`loci_names()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_names.md)
  :

  Get the names of loci in a `gen_tibble`

- [`loci_pi()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_pi.md)
  : Estimate nucleotide diversity (pi) at each locus

- [`loci_transitions()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_transitions.md)
  : Find transitions

- [`loci_transversions()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_transversions.md)
  : Find transversions

## Individual

Functions operating across individuals.

- [`indiv_het_obs()`](https://evolecolgroup.github.io/tidypopgen/reference/indiv_het_obs.md)
  : Estimate individual observed heterozygosity
- [`indiv_inbreeding()`](https://evolecolgroup.github.io/tidypopgen/reference/indiv_inbreeding.md)
  : Individual inbreeding coefficient
- [`indiv_missingness()`](https://evolecolgroup.github.io/tidypopgen/reference/indiv_missingness.md)
  : Estimate individual missingness
- [`indiv_ploidy()`](https://evolecolgroup.github.io/tidypopgen/reference/indiv_ploidy.md)
  : Return individual ploidy

## Population

Functions operating on a `gen_tibble` grouped by population.

- [`pop_fis()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_fis.md)
  : Compute population specific FIS
- [`pop_fst()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_fst.md)
  : Compute population specific Fst
- [`pop_global_stats()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_global_stats.md)
  : Compute basic population global statistics
- [`pop_het_exp()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_het_exp.md)
  [`pop_gene_div()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_het_exp.md)
  : Compute the population expected heterozygosity
- [`pop_het_obs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_het_obs.md)
  : Compute the population observed heterozygosity
- [`pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_tajimas_d.md)
  : Estimate Tajima's D for the whole genome

## Pairwise

Functions comparing individuals or populations.

- [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_allele_sharing.md)
  :

  Compute the Pairwise Allele Sharing Matrix for a `gen_tibble` object

- [`pairwise_grm()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_grm.md)
  :

  Compute the Genomic Relationship Matrix for a `gen_tibble` object

- [`pairwise_ibs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_ibs.md)
  :

  Compute the Identity by State Matrix for a `gen_tibble` object

- [`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_king.md)
  :

  Compute the KING-robust Matrix for a `gen_tibble` object

- [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_fst.md)
  : Compute pairwise population Fst

- [`nwise_pop_pbs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_pbs.md)
  : Compute the Population Branch Statistics for each combination of
  populations

## Windows

Functions to compute statistics in sliding windows across loci.

- [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)
  [`gt_roh_window()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)
  : Detect runs of homozygosity using a sliding-window approach
- [`windows_nwise_pop_pbs()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_nwise_pop_pbs.md)
  : Compute the Population Branch Statistics over a sliding window
- [`windows_pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_pairwise_pop_fst.md)
  : Compute pairwise Fst for a sliding window
- [`windows_pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_pop_tajimas_d.md)
  : Compute Tajima's D for a sliding window
- [`windows_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_stats_generic.md)
  : Estimate window statistics from per locus estimates

## PCA

Functions for performing Principal Component Analysis.

- [`gt_pca`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md)
  :

  Principal Component Analysis for `gen_tibble` objects

- [`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_autoSVD.md)
  :

  PCA controlling for LD for `gen_tibble` objects

- [`gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_partialSVD.md)
  :

  PCA for `gen_tibble` objects by partial SVD

- [`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_randomSVD.md)
  :

  PCA for `gen_tibble` objects by randomized partial SVD

- [`predict(`*`<gt_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/predict_gt_pca.md)
  : Predict scores of a PCA

- [`tidy(`*`<gt_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/tidy_gt_pca.md)
  :

  Tidy a `gt_pca` object

- [`augment(`*`<gt_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/augment_gt_pca.md)
  : Augment data with information from a gt_pca object

- [`augment_loci(`*`<gt_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/augment_loci_gt_pca.md)
  : Augment the loci table with information from a gt_pca object

- [`autoplot(`*`<gt_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_gt_pca.md)
  :

  Autoplots for `gt_pca` objects

## DAPC

Functions for performing Discriminant Analysis of Principal Components.

- [`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_dapc.md)
  : Discriminant Analysis of Principal Components for gen_tibble

- [`gt_cluster_pca()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_cluster_pca.md)
  : Run K-clustering on principal components

- [`gt_cluster_pca_best_k()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_cluster_pca_best_k.md)
  : Find the best number of clusters based on principal components

- [`tidy(`*`<gt_dapc>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/tidy.gt_dapc.md)
  :

  Tidy a `gt_dapc` object

- [`augment(`*`<gt_dapc>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/augment.gt_dapc.md)
  : Augment data with information from a gt_dapc object

- [`autoplot(`*`<gt_dapc>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_gt_dapc.md)
  :

  Autoplots for `gt_dapc` objects

- [`autoplot(`*`<gt_cluster_pca>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_gt_cluster_pca.md)
  :

  Autoplots for `gt_cluster_pca` objects

## K-clustering and ADMXITURE

Functions for creating, tidying, and visualising `gt_admix` and
`q_matrix` objects.

- [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.md)
  : Run ADMIXTURE from R

- [`gt_snmf()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_snmf.md)
  : Run SNMF from R in tidypopgen

- [`autoplot(`*`<gt_admix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_gt_admix.md)
  :

  Autoplots for `gt_admix` objects

- [`summary(`*`<gt_admix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/summary.gt_admix.md)
  : Summary method for gt_admix objects

- [`c(`*`<gt_admix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/c.gt_admix.md)
  : Combine method for gt_admix objects

- [`gt_admix_reorder_q()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admix_reorder_q.md)
  : Reorder the q matrices based on the grouping variable

- [`q_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/q_matrix.md)
  :

  Convert a standard matrix to a `q_matrix` object

- [`read_q_files()`](https://evolecolgroup.github.io/tidypopgen/reference/read_q_files.md)
  :

  Read and structure .Q files or existing matrices as `q_matrix` or
  `gt_admix` objects.

- [`get_q_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/get_q_matrix.md)
  :

  Return a single Q matrix from a `gt_admix` object

- [`get_p_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/get_p_matrix.md)
  :

  Return a single P matrix from a `gt_admix` object

- [`tidy(`*`<q_matrix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/tidy.q_matrix.md)
  : Tidy a Q matrix

- [`augment(`*`<q_matrix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/augment_q_matrix.md)
  : Augment data with information from a q_matrix object

- [`autoplot(`*`<q_matrix>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_q_matrix.md)
  :

  Autoplots for `q_matrix` objects

- [`distruct_colours`](https://evolecolgroup.github.io/tidypopgen/reference/distruct_colours.md)
  : Distruct colours

- [`scale_fill_distruct()`](https://evolecolgroup.github.io/tidypopgen/reference/scale_fill_distruct.md)
  : Scale constructor using the distruct colours

- [`theme_distruct()`](https://evolecolgroup.github.io/tidypopgen/reference/theme_distruct.md)
  : A theme to match the output of distruct

## Other analyses (e.g. ADMIXTOOLS2, pcadapt)

Functions for other analyses within `R`.

- [`gt_extract_f2()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_extract_f2.md)
  : Compute and store blocked f2 statistics for ADMIXTOOLS 2

- [`gt_pcadapt()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pcadapt.md)
  :

  pcadapt analysis on a `gen_tibble` object

- [`autoplot(`*`<gt_pcadapt>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_gt_pcadapt.md)
  :

  Autoplots for `gt_pcadapt` objects

## `bigSNP`

Mostly for developers: Functions for operating directly on `bigSNP`
objects.

- [`snp_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/reference/snp_allele_sharing.md)
  : Compute the Pairwise Allele Sharing Matrix for a bigSNP object
- [`snp_ibs()`](https://evolecolgroup.github.io/tidypopgen/reference/snp_ibs.md)
  : Compute the Identity by State Matrix for a bigSNP object
- [`snp_king()`](https://evolecolgroup.github.io/tidypopgen/reference/snp_king.md)
  : Compute the KING-robust Matrix for a bigSNP object

## Example data

Function for loading an example `gen_tibble`.

- [`load_example_gt()`](https://evolecolgroup.github.io/tidypopgen/reference/load_example_gt.md)
  : Load example gen_tibble

## tidyverse methods

Methods for tidyverse verbs for `gen_tibble` objects

- [`arrange(`*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/arrange.gen_tbl.md)
  :

  An arrange method for `gen_tibble` objects

- [`arrange(`*`<grouped_gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/arrange.grouped_gen_tbl.md)
  :

  An arrange method for grouped `gen_tibble` objects

- [`filter(`*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/filter.gen_tbl.md)
  : Tidyverse methods for gt objects

- [`filter(`*`<grouped_gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/filter.grouped_gen_tbl.md)
  :

  A filter method for grouped `gen_tibble` objects

- [`mutate(`*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/mutate.gen_tbl.md)
  :

  A mutate method for `gen_tibble` objects

- [`mutate(`*`<grouped_gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/mutate.grouped_gen_tbl.md)
  :

  A mutate method for grouped `gen_tibble` objects

- [`` `$<-`( ``*`<gen_tbl>`*`)`](https://evolecolgroup.github.io/tidypopgen/reference/cash-set-.gen_tbl.md)
  :

  A \$ method for `gen_tibble` objects
