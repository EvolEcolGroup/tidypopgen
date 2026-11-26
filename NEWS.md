# tidypopgen 1.1.0

## New features and changes

* new `gt_set_chromosome_types()`, `gt_chromosome_types()` and
  `loci_chromosome_types()` tell autosomes from sex chromosomes (X, Y, Z, W),
  pseudo-autosomal regions and mtDNA, recognised from their names or declared
  (e.g. PLINK's numeric codes with `chr_set = "human"`)
* the VCF readers (`gen_tibble()`, `gt_merge_vcfs()`, `gt_concat_vcfs()`)
  store a haploid genotype of a diploid individual (e.g. the X of a male) as a
  homozygote (0 or 2), as PLINK does, rather than as 0 or 1; haploid genotypes
  on sex chromosomes and mtDNA no longer set the ploidy of individuals
* `gen_tibble(phased = TRUE)` reads haploid genotypes as homozygotes
* `gen_tibble()` keeps the sex of individuals from packedancestrymap `.ind`
  files
* new `indiv_sex_check()` checks the recorded sex of individuals against
  their X (and Y) genotypes, as `plink --check-sex` (PLINK 1.9)
* functions that assume diploid autosomes (e.g. `gt_pca()`, `pairwise_king()`,
  `indiv_het_obs()`, `pop_fst()`, `gt_extract_f2()`) use autosomal loci only
  and say how many loci on X, Y or MT they leave out, as PLINK; declare a
  chromosome as an autosome with `gt_set_chromosome_types()` to keep it
* new `gt_pca()` computes PCAs, by default with missing genotypes handled as
  `plink --pca` (PLINK 1.9), without imputation; `gt_pca_partialSVD()` and
  `gt_pca_randomSVD()` are deprecated
* `gt_pca()` accepts pseudohaploid data, coded as homozygous diploids (as in
  PLINK)
* `gt_pcadapt()`, `gt_dapc()` and `loci_ld_clump()` work with missing
  genotypes without imputation, using the individuals genotyped at each locus
  (as `pcadapt` and PLINK)
* new `loci_ld_prune()` prunes loci as `plink --indep-pairwise` (PLINK 1.9)
* `predict()` of a `gt_pca` uses `project_method = "least_squares"` on all
  components by default, and is 2-3 times faster; new
  `project_method = "plink"` projects as PLINK 1.9
* `pairwise_grm(method = "plink")` gives PLINK 1.9's relationship matrix
* `pairwise_ibs(missing_correction = "allele_freq")` corrects for missing
  genotypes as PLINK 1.9's `--distance`
* new `pairwise_dist()` computes genetic distances between individuals as
  ngsDist, e.g. for neighbour-joining trees
* `pairwise_king()` on a grouped `gen_tibble` computes kinship within each
  group (it used to ignore the groups)
* `pairwise_king()`, `pairwise_ibs()`, `pairwise_allele_sharing()`,
  `pairwise_dist()`, `indiv_inbreeding()` and `windows_indiv_roh()` are S3
  generics, so that other packages can add methods
* new `bootstrap_blocks()`, for developers: the blocks of the block bootstrap,
  so that other packages can resample the same blocks
* new `group_labels()` and `group_columns()`, for developers: the labels and
  tidy columns of groups, so that other packages name groups as tidypopgen
* new `gt_merge_vcfs()` and `gt_concat_vcfs()` read several VCF files into a
  single `gen_tibble`, as `bcftools merge` and `bcftools concat`
* integer chromosomes (`chromosomes_as_int`) for datasets with many
  chromosomes; the original names of integer loci and chromosomes are kept in
  sidecar files (`gt_lookup_loci()`, `gt_restore_loci()` and
  `gt_lookup_chromosomes()`)
* phased genotypes (experimental): `gen_tibble(phased = TRUE)` reads phased
  VCF files; `is_phased()` and `show_haplotypes()` give the phase,
  `gt_as_vcf()` writes it, and `rbind()` keeps it
* new `gt_simulate_bn()` simulates toy datasets (Balding-Nichols model)
* `qc_report_loci(hwe = FALSE)` gives the MAF and missingness of loci in a
  single pass
* `gt_impute_simple()` works after removing loci (but not individuals)
* `gen_tibble()` reads *packedancestry* genotype files in the TGENO layout
  (one record per individual, as in recent releases of the AADR) as well as
  GENO; the layout is detected from the header of the file, whatever its
  extension (*.geno* or *.tgeno*)
* faster input and output: VCF, *.bed* and *packedancestry* (*.geno*) files
  are read in parallel (`n_cores`), *.ped* files are read faster, and *.gt*
  files are saved twice as fast
* faster `pairwise_king()` (up to 1.8 times on CPUs with AVX2), `pairwise_ibs()`
  and `pairwise_allele_sharing()`, and much faster `indiv_inbreeding()` and
  `qc_report_indiv()` on grouped data
* faster `select_loci()`, and faster `gt_as_vcf()`, `gt_update_backingfile()`
  and `gt_order_loci()` with many chromosomes; `gt_as_genlight()` uses less
  memory
* `indiv_inbreeding()` on a grouped `gen_tibble` returns one value per
  individual, in the order of the rows (computed within each group); the
  list of one vector per group (`type = "list"`) is deprecated. It also
  works on the `genotypes` column, e.g. within `mutate()`
* `indiv_het_obs()`, `indiv_missingness()` and `indiv_ploidy()` have methods
  for grouped `gen_tibble`s (all individuals in one pass, in the order of the
  rows) and for tibbles that lost their class
* `loci_alt_freq(as_counts = TRUE)` on a grouped `gen_tibble` follows `type`:
  by default a tibble with the columns `n_alt` and `n_valid`;
  `type = "matrix"` gives the two matrices of before, now with loci and group
  names
* new article on large RAD-seq datasets; the big data article uses
  `gt_merge_vcfs()`

## Bug fixes

* fix occasional "File '...bk' already exists" errors from functions using
  temporary backing files (PCA, LD clumping, KING, IBS and allele sharing)
* `loci_ld_clump()`, `gt_pcadapt()` and the PCA functions: fix false "You
  can't have missing values" errors (or hangs) after removing loci or
  individuals with missing genotypes
* `gt_pca(method = "random")` and `gt_pca_autoSVD()`: fix hangs or cryptic
  errors with a custom `fun_scaling` and missing genotypes, and use at most
  `n_cores` cores
* `gt_pcadapt()`: clear error for loci with the same genotype in every
  individual (#179)
* `gen_tibble()`: error on *.bed* files whose size does not match the *.bim*
  and *.fam* files, and on individual-major *.bed* files
* `gen_tibble()`: *packedancestry* files that are truncated, or whose header
  does not match the *.ind* and *.snp* files, give an error before the
  backing file is created, and EIGENSTRAT (text) *.geno* files a clear error
* `gen_tibble()` on *.ped* files: a third allele always flags the locus as
  multi-allelic
* `gen_tibble()` with the "cpp" parser: fix malformed VCF files read as wrong
  genotypes (now an error), missing chromosomes (now `NA`), and backing files
  left behind on errors
* `gt_as_vcf()`: fix loci with missing chromosomes (now written as ".")
* `loci_ld_clump()`: error on `exclude` values that are not indices of loci
* `rbind()`: fix grouped `gen_tibble`s, which now keep the grouping of the
  first one
* `gt_add_sf()`: no "Adding missing grouping variables" message on grouped
  `gen_tibble`s
* `loci_pi()` on a grouped `gen_tibble` with a single individual: fix the
  output, which was a plain vector whatever `type`

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
