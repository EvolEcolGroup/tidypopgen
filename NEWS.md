# tidypopgen `dev`
* fix printing a `gen_tibble` with no loci, and `format()` of a genotype column
  with no individuals returning a non-empty vector
* fix chromosome names being renumbered in lexicographic order (e.g. chromosome
  10 coded as 2) when some chromosome names are not numeric, which affected
  `gt_as_plink(chromosomes_as_int = TRUE)`, `gt_order_loci()` and
  `gt_pca_autoSVD()`
* fix `gt_pseudohaploid()` resetting the ploidy of individuals in the backing
  matrix that are not in the `gen_tibble`, and make it error on data that are
  not diploid
* fix `tidy()` and `augment()` of a `q_matrix` ignoring its individual ids when
  it has no group information, which misassigned proportions when `data` was in
  a different order
* correct the documentation of the columns returned by `tidy()` for the loadings
  of `gt_dapc` objects
* fix `gt_dapc()` using k instead of k - 1 principal components by default when
  the number of available components equals the number of groups k
* fix `gt_dapc()` reporting the proportion of variance retained by the PCs
  (`$var`) from singular values instead of eigenvalues
* `gt_dapc()` now gives an informative error when `pop` is not a factor or a
  value of k explored by `gt_cluster_pca()`
* fix `gt_cluster_pca_best_k()` returning the position of the best k rather than
  its value, and `gt_dapc()` selecting the wrong clustering when `k_clusters`
  did not start at 1
* fix `gt_cluster_pca()` failing when only one value of k above 1 (or only k =
  1) is explored, which includes the default `k_clusters` for datasets with 25
  or fewer individuals, and when `n_pca = 1`
* `augment_loci()` now works on `gt_pca_autoSVD()` results, giving NA loadings
  to the loci removed by LD pruning
* fix `augment()` for `gt_pca` and `gt_dapc` objects (and `augment_loci()` for
  `gt_pca`) failing when a single component is requested, when individual ids
  are integers, or when ids do not match (now giving an informative error)
* fix `gt_pca_autoSVD()` mislabelling loadings, misreporting variance explained
  and breaking `predict()` when the loci table is not in the order of the
  backing file
* fix `windows_pairwise_pop_fst()` and `windows_nwise_pop_pbs()` always using
  Hudson's Fst, ignoring `method` and `fst_method` respectively
* `pop_fst()` and `pop_fis(method = "WG17")` now give an informative error,
  instead of "subscript out of bounds", when there is only one group
* fix `loci_alt_freq()` on polyploids ignoring `as_counts = TRUE` and failing
  when only one locus is selected
* fix `loci_alt_freq()` and `loci_maf()` returning all NA for a single polyploid
  individual
* fix `loci_alt_freq(as_counts = TRUE)` on a single individual returning NA
  alternate and non-zero valid allele counts for missing genotypes, and unnamed
  columns
* fix `pop_tajimas_d()` and `windows_pop_tajimas_d()` returning NA when a locus
  has no valid genotypes, and ignoring segregating sites with pi equal to 1 (a
  single heterozygote genotyped)
* fix `gt_as_vcf()` and `gt_as_plink(type = "ped")` writing positions stored as
  doubles in scientific notation (e.g. `1e+05`)
* fix `gt_as_plink(type = "bed")` not writing the population (as family ID),
  parents, sex and phenotype to the .fam file, and `gt_as_plink(type =
  "ped")`/`type = "raw"` ignoring the grouping variable as family ID
* fix `gt_as_plink(type = "ped")` and `gt_as_plink(type = "raw")` dropping
  parents (`paternal_ID`/`maternal_ID`, or `father`/`mother`) and numeric sex
  codes read from PLINK files
* fix `gt_as_plink(type = "ped")` writing a genetic distance of 0 for all loci
  in the .map file
* fix `gt_as_plink()` failing to write ped and raw files when the phenotype
  column is numeric
* fix `gt_as_plink()` failing to write ped and raw files for `gen_tibble`s
  without a population column (e.g. read from a vcf)
* fix `gt_as_hierfstat()` turning numeric-like population labels (e.g. "1", "2")
  into NA, and keep loci names as column names
* fix the `scale_fill_distruct()` example, and the documentation of
  `qc_report_indiv()` and its `autoplot()` method
* fix the description of the upset plot in the quality control vignette
* the `heatmap_dendro` object returned by `heatmap_add_dendro()` is now a proper
  `gtable`, so it can be combined with other plots using patchwork (previously
  this failed on printing or saving)
* fix `heatmap_add_dendro()` drawing dendrograms that did not match the order of
  the heatmap when there are tied distances
* fix `heatmap_add_dendro()` failing with ggplot2 versions earlier than 4.0.0
* fix `heatmap_pairwise()` (and `autoplot()` of tidy pairwise results) always
  failing with "Duplicate pairs found" when the id columns are factors, e.g.
  from `pairwise_pop_fst()` on data grouped by a factor
* fix grouped `loci_alt_freq()`, `loci_maf()`, `loci_hwe()`,
  `loci_missingness()`, `loci_pi()` and `pop_tajimas_d()` failing when the last
  group is empty (factor grouping with `.drop = FALSE`)
* fix `gt_as_plink()`, `gt_as_geno_lea()` and `gt_as_hierfstat()` (and hence
  `gt_admixture()` and `gt_snmf()`) erroring on pseudohaploid data, which they
  document as supported
* fix grouped `qc_report_loci()` returning HWE p-values greater than 1 (the
  Bonferroni-corrected p-value is now capped at 1)
* grouped `qc_report_indiv()` now gives informative errors when data are grouped
  by more than one variable, or when `kings_threshold` is used with
  pseudohaploid data (previously silently ignored)
* `autoplot()` of a `qc_report_indiv` with `type = "relatedness"` now gives an
  informative error if the report was created without `kings_threshold`
* fix `qc_report_indiv()` on grouped `gen_tibble`s attaching `to_keep` flags to
  the wrong individuals when groups are not contiguous, and failing when
  grouping by a factor with unused levels or a numeric variable
* fix grouped `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`,
  `loci_hwe()` and `loci_pi()` failing when the grouping variable is a factor,
  or numeric with values other than 1 to the number of groups
* fix `pop_global_stats()` (and the global values of `pop_het_exp()` and
  `pop_fis()`) returning NaN at loci where a population has no valid genotypes,
  which also biased the genome-wide estimates
* fix `rbind()` attaching individual ploidies to the wrong individuals (and
  mislabelling the merged ploidy) when the input `gen_tibble`s had been filtered
  or reordered
* fix `loci_transversions()` and `loci_transitions()` misclassifying A<->C and
  G<->T substitutions as transitions
* fix `gen_tibble()` on VCF files with partially missing genotypes (e.g. "1/."):
  the "cpp" parser silently coded them as non-missing dosages, and the "vcfR"
  parser failed; they are now coded as missing
* fix `gen_tibble()` with `parser = "vcfR"` failing when a sample has a missing
  genotype (".") on the second variant of the VCF
* fix `gt_update_backingfile()` (and `gt_order_loci()`) failing when the number
  of loci is an exact multiple of `chunk_size`
* fix `snp_allele_sharing()` (and thus `pairwise_allele_sharing()`,
  `pairwise_grm()`, `pop_fst()` and `indiv_inbreeding()`), which overestimated
  allele sharing when loci were processed in blocks of different sizes
* fix `indiv_het_obs()` to divide by the number of selected loci rather than
  the number of loci in the backing file (results were deflated after
  `select_loci()`)
* fix the order of the loci in the P matrices returned by `gt_admixture()`:
  ADMIXTURE sorts loci by chromosome and position, so rows did not match the
  loci of the input when these were not already sorted (e.g. when chromosome
  names are not numeric)
* implement an `autoplot()` method for distance matrices created by `pairwise_*`
  functions. This requires the matrices to have an appropriate class, so it will
  not work with objects generated by older versions of `tidypopgen`
* add option to use numerical identifiers for loci names when working with 
  large datasets, to reduce memory usage and speed up processing. This is
  useful when working with full genomes (tens of millions of positions).
* generalise several functions to work with polyploid data: `indiv_het_obs()`,
  `pop_het_obs()`, `loci_pi()` and `gt_as_genind()`.
* `gt_as_hierfstat()`, `gt_as_geno_lea()` and `gt_as_plink()` now error clearly
  on polyploid data, as `hierfstat`, LEA's `.geno` format and PLINK's bed/ped/
  raw formats only support haploid or diploid data.

# tidypopgen 0.4.4
* fix seed argument of `gt_admixture()`
* implement our own upset plot function to avoid dependency on `UpSetR`
  package, which is no longer maintained
* fix reading of VCFs when genotype separator symbol is present in the
  sample information
* fix `loci_` functions output when calculating across groups for a single loci
* fix `gt_update_backingfile()` to preserve imputation metadata when backing
  files are updated
* add validation in `gt_impute_simple()` to detect and report SNP dimension
  mismatches caused by excluded SNPs
* add `outdir` argument in `gt_admixture()` to store partial outputs for large jobs
* fix `gt_admixture()` output file names when n_runs > 1

# tidypopgen 0.4.3
* fix occasional test failure when temporary file names contained
  certain patterns

# tidypopgen 0.4.2
* allow least squares projection of multiple components at once in
  `gt_project_pca()`

# tidypopgen 0.4.1
* duplicate-distance checks removed in `is_loci_table_ordered()`
* make `qc_report_*` functions work with pseudohaploids

# tidypopgen 0.4.0
* store chromosomes as factors, and remove redundant column storing
  chromosome as integer

# tidypopgen 0.3.3
* only use the FBM backend (rather than a complete bigsnpr object)
* fix a few minor bugs

# tidypopgen 0.3.2
* tidy up for CRAN submission

# tidypopgen 0.3.0
* add `gt_pseudohaploid()` and adapt appropriate functions to work with
  pseudohaploids

# tidypopgen 0.2.0
* Lots of optimisations to speed up processing of large datasets.
* Added an `sf` backend to `tidypopgen` for spatial data.
* Windowed statistics.

# tidypopgen 0.1.0
* Initial public version of `tidypopgen` released on r-universe.
