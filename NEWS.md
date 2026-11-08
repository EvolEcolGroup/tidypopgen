# tidypopgen (development version)

* `gt_as_vcf()` is much faster, and uses less memory, on data with many
  chromosomes (e.g. RAD loci assembled de novo, with one chromosome per locus)
* `gt_as_vcf()`: fix loci with missing chromosomes (now written as ".")
* `gt_as_genlight()` reads the genotypes in chunks of loci (new argument
  `chunk_size`), instead of loading them all in memory
* `select_loci()` resolves function calls that give positions or names (e.g.
  `which(...)` or `sample(...)`) and single values directly, without the slow
  general `tidyselect` machinery
* `select_loci(.engine = "tidyselect")` always uses `tidyselect`, in case a
  selection is not resolved by the faster path as `tidyselect` would
* `gen_tibble()` reads VCF files in a single pass with the "cpp" parser
  (decompressing a *.vcf.gz* file only once), and gives a clear error for VCF
  files with no biallelic SNPs
* `gen_tibble()` uses `n_cores` (previously ignored) to read VCF files with
  the "cpp" parser in parallel: the genotypes are parsed in parallel, and
  *.vcf.gz* files compressed with `bgzip` are also decompressed in parallel
* `gen_tibble()` reads VCF files with many chromosomes (e.g. RAD loci
  assembled de novo) faster, sorting the chromosome names in C order when
  this gives the order of the locale
* `gen_tibble()` and `gt_save()` save the *.gt* file twice as fast, with the
  fastest gzip compression
* `gen_tibble()` on VCF files: the "cpp" parser reads a missing chromosome
  (".") as `NA`, as the "vcfR" parser does

# tidypopgen 1.0.0

## New features and changes

* optimisations to work on whole genomes, both speed and memory
* new functions to compute site frequency spectra (SFS) for demographic
  inference: `pop_sfs()`,`pairwise_pop_sfs()` and `nwise_pop_sfs()`
* native `windows_indiv_roh()` matching the algorithm used by PLINK 1.9 
  (`--homozyg`); the old version, based on `detectRUNS`, is deprecated as
  `gt_roh_window()`
* `gt_roh_froh()`, `gt_roh_loci()` and `gt_roh_islands()` summarise runs of
  homozygosity (F_ROH, loci in runs, and ROH islands)
* `jackknife_*` functions for genome-wide statistics with block jackknife
  standard errors
* `pop_tajimas_d()` and related functions handle missing genotypes as `pixy`,
  and work with polyploid and pseudohaploid data
* `autoplot()` method for the distance matrices of `pairwise_*` functions
* option to use integer loci names (`names_as_int`) for large datasets, when
  reading VCF and PLINK bed files
* polyploid support in `indiv_het_obs()`, `pop_het_obs()`, `loci_pi()` and
  `gt_as_genind()`; clear errors where polyploids are not supported
* `pairwise_ibs(type = "adjusted_counts")` now returns PLINK's allele count
  distance (`--distance flat-missing`)
* `pairwise_*` functions use imputed genotypes when set with
  `gt_set_imputed()`, and use several cores by default
* `gen_tibble()` on ped files stores individual information in the same
  columns as for bed files (see `?gen_tibble`)
* `gen_tibble()` copies the backing file of bigSNP objects (never modifying
  them), and works on read-only files
* `rbind(use_position = TRUE)` keeps the loci names of the reference, and
  needs much less memory; integer names are matched by position
* the `as_is` argument of `rbind()` is deprecated (it was never implemented)
* `cbind()` renames duplicated columns as `dplyr::bind_cols()`; `id` and
  `genotypes` can not be added
* `distruct_colours` holds the 56 distinct distruct colours, recycled with a
  warning by `scale_fill_distruct()`
* `gt_as_genlight()` keeps alleles, chromosomes and positions, so that
  `gt_from_genlight()` can convert back
* classes are preserved by all dplyr verbs; grouped functions keep the order
  of factor levels, and work with several grouping variables and empty groups
* `heatmap_add_dendro()` returns a `gtable`, which can be combined with
  patchwork
* `augment_loci()` works on `gt_pca_autoSVD()` results
* `gt_admixture()` uses ADMIXTURE on the PATH when conda is not installed
* faster "cpp" VCF parser, and `select_loci_if()` no longer loads all the
  genotypes
* informative errors for invalid inputs to `gen_tibble()`, `show_loci<-`,
  `gt_order_loci()`, `gt_pcadapt()`, `gt_dapc()`, `pop_fst()` and others

## Bug fixes

* `gen_tibble()`: fix backing file names (truncated at dots, or ".vcf.bk"),
  and backing files left behind on errors
* `gen_tibble()`: fix files with trailing blank lines, CRLF line endings or
  upper case extensions
* `gen_tibble()` on bed files keeps .fam IDs as text, and quantitative
  phenotypes as numbers
* `gen_tibble()` on VCF files: fix the ploidy of individuals missing on the
  first variant, partially missing genotypes, and genotypes above the ploidy
* `gen_tibble()`: "0" and "." are always missing alleles; fix the "vcfR"
  parser on indels and on missing genotypes on the second variant
* `gt_load()`: fix the ploidies of gen_tibbles saved in the old bigSNP format
* `gt_as_plink()`: fix the .fam and .map information, missing alleles, and
  ped/raw files with numeric phenotypes, no population or pseudohaploids
* `gt_as_vcf()` and `gt_as_plink()`: fix positions in scientific notation, and
  the date in the VCF header
* `gt_admixture()` and `gt_snmf()`: fix seeds with several runs, paths
  containing "geno", and the order of loci in P matrices
* `gt_snmf()`: fix continued projects, and the order of cross-entropies
* `gt_impute_xgboost()`: fix unsorted loci; `gt_impute_simple()`: restore the
  BLAS option
* `filter_high_relatedness()`: fix missing and strongly negative kinship
  coefficients
* `snp_allele_sharing()`: fix blocks of different sizes; `indiv_het_obs()`:
  fix the number of loci after `select_loci()`
* `pop_global_stats()`: fix loci where a population has no data;
  `pop_tajimas_d()`: fix loci with no valid genotypes
* `loci_alt_freq()` and `loci_maf()`: fix polyploids and single individuals
* grouped `loci_*` functions and `qc_report_*`: fix factor, numeric,
  non-contiguous and empty groups, and p-values above 1
* `loci_transitions()` and `loci_transversions()`: fix misclassified
  substitutions
* `windows_*` functions: fix the last loci of chromosomes;
  `windows_pairwise_pop_fst()` and `windows_nwise_pop_pbs()` ignoring `method`
* `gt_roh_window()`: fix homozygotes for the alternate allele with
  `heterozygosity = TRUE`
* `gt_dapc()` and `gt_cluster_pca()`: fix the default number of components,
  the variance retained, and the choice of k
* `gt_pca_autoSVD()`: fix unsorted loci; `augment()` of `gt_pca` and
  `gt_dapc`: fix single components and integer ids
* `rbind()`: fix ploidies after filtering, and missing `allele_ref`;
  `rbind_dry_run()`: fix paths to .bim files
* dplyr verbs on grouped gen_tibbles: fix `select()`, `rename()`,
  `relocate()`, `[`, `ungroup()`, `group_by()` and `mutate()`
* `autoplot()` of `q_matrix` and `gt_admix`: fix the sorting options, numeric
  and unused groups, and `type = "cv"`; ordering in other `autoplot()`s
* `gt_admix_reorder_q()`: fix K = 1
* `heatmap_add_dendro()`: fix ties and older ggplot2; `heatmap_pairwise()`:
  fix factor ids
* fix chromosome names being renumbered in lexicographic order
* `gt_pseudohaploid()`: fix other individuals of the backing file, and error
  on data that are not diploid
* `tidy()` and `augment()` of `q_matrix`: fix ignoring individual ids
* `gt_as_genind()`: fix loci names with "."; `gt_from_genlight()`: fix
  objects without `pop`; `gt_as_hierfstat()`: fix numeric labels
* `gt_update_backingfile()`: fix numbers of loci that are multiples of
  `chunk_size`; fix printing gen_tibbles with no loci
* fix documentation and error messages of several functions

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
