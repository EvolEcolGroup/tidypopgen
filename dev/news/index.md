# Changelog

## tidypopgen (development version)

- `gen_tibble(names_as_int = TRUE)` on VCF files writes the original
  loci names to a sidecar file (*.loci.tsv*); new
  [`gt_lookup_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  and
  [`gt_restore_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  retrieve them.
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md),
  [`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md),
  [`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md),
  [`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_order_loci.md)
  and [`rbind()`](https://rdrr.io/r/base/cbind.html) keep the sidecar
- new argument `chromosomes_as_int` of
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
  and
  [`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md)
  stores chromosomes as integer codes, saving memory with millions of
  chromosomes (e.g. RAD loci assembled de novo); new
  [`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  translates codes back into names
- new article on large RAD-seq datasets
- new
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
  reads several VCF files with different samples into a single
  `gen_tibble`, as `bcftools merge` would, in one go (without reading
  each file into its own `gen_tibble` and binding them in pairs with
  [`rbind()`](https://rdrr.io/r/base/cbind.html)): loci are matched by
  chromosome and position, keeping the sites found in all files
  (`variants = "intersection"`) or all sites (`"union"`, with missing or
  homozygous reference genotypes for the samples of the files without a
  site); loci with different alleles are dropped, kept apart, or
  swapped, following the options of `bcftools merge`. Unlike
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  it names loci with integers by default (`names_as_int = TRUE`), to
  save memory with large datasets
- new
  [`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md)
  reads several VCF files with the same samples and different loci
  (e.g. one file per chromosome) into a single `gen_tibble`, as
  `bcftools concat` would, in one go; the samples can be in a different
  order in each file. As
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md),
  it names loci with integers by default (`names_as_int = TRUE`)
- the big data article merges the VCFs of the SGDP individuals with
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md),
  rather than with `bcftools merge`
- [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md)
  is much faster, and uses less memory, on data with many chromosomes
  (e.g. RAD loci assembled de novo, with one chromosome per locus)
- [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md):
  fix loci with missing chromosomes (now written as “.”)
- [`gt_as_genlight()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_genlight.md)
  reads the genotypes in chunks of loci (new argument `chunk_size`),
  instead of loading them all in memory
- [`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md)
  resolves function calls that give positions or names (e.g.
  `which(...)` or `sample(...)`) and single values directly, without the
  slow general `tidyselect` machinery
- [`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md)
  resolves a variable that holds positions or names (e.g.
  `select_loci(x, idx)`) directly, as it does `all_of(idx)`: with
  millions of loci, `tidyselect` took minutes and gigabytes of memory.
  As in `tidyselect`, a locus named like the variable is selected
  instead, and the deprecation warning of `tidyselect` (in favour of
  `all_of(idx)`) is kept
- `select_loci(.engine = "tidyselect")` always uses `tidyselect`, in
  case a selection is not resolved by the faster path as `tidyselect`
  would
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  reads VCF files in a single pass with the “cpp” parser (decompressing
  a *.vcf.gz* file only once), and gives a clear error for VCF files
  with no biallelic SNPs
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  uses `n_cores` (previously ignored) to read VCF files with the “cpp”
  parser in parallel: the genotypes are parsed in parallel, and
  *.vcf.gz* files compressed with `bgzip` are also decompressed in
  parallel
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  with the “cpp” parser gives a clear error, and leaves no backing
  files, when the first variant of a VCF does not have a genotype for
  each sample in its header
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  with the “cpp” parser gives a clear error, rather than reading wrong
  genotypes, on malformed VCF files: genotypes of biallelic SNPs with
  alleles other than 0, 1 or “.”, or empty, a FORMAT that does not start
  with GT, positions that are not integers (e.g. “12x”), and corrupted
  BGZF blocks; and on an invalid `n_cores` (which could crash R)
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  reads VCF files with many chromosomes (e.g. RAD loci assembled de
  novo) faster: the “cpp” parser looks up chromosome names without
  allocating memory for each, and the chromosome names are sorted in C
  order when this gives the order of the locale
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  and
  [`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md)
  save the *.gt* file twice as fast, with the fastest gzip compression
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  on VCF files: the “cpp” parser reads a missing chromosome (“.”) as
  `NA`, as the “vcfR” parser does

## tidypopgen 1.0.0

CRAN release: 2026-11-07

### New features and changes

- optimisations to work on whole genomes, both speed and memory
- new functions to compute site frequency spectra (SFS) for demographic
  inference:
  [`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_sfs.md),[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_sfs.md)
  and
  [`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/nwise_pop_sfs.md)
- native
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)
  matching the algorithm used by PLINK 1.9 (`--homozyg`); the old
  version, based on `detectRUNS`, is deprecated as
  [`gt_roh_window()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_window.md)
- [`gt_roh_froh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_froh.md),
  [`gt_roh_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_loci.md)
  and
  [`gt_roh_islands()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_islands.md)
  summarise runs of homozygosity (F_ROH, loci in runs, and ROH islands)
- `jackknife_*` functions for genome-wide statistics with block
  jackknife standard errors
- [`pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_tajimas_d.md)
  and related functions handle missing genotypes as `pixy`, and work
  with polyploid and pseudohaploid data
- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method for the distance matrices of `pairwise_*` functions
- option to use integer loci names (`names_as_int`) for large datasets,
  when reading VCF and PLINK bed files
- polyploid support in
  [`indiv_het_obs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_het_obs.md),
  [`pop_het_obs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_het_obs.md),
  [`loci_pi()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_pi.md)
  and
  [`gt_as_genind()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_genind.md);
  clear errors where polyploids are not supported
- `pairwise_ibs(type = "adjusted_counts")` now returns PLINK’s allele
  count distance (`--distance flat-missing`)
- `pairwise_*` functions use imputed genotypes when set with
  [`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_set_imputed.md),
  and use several cores by default
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  on ped files stores individual information in the same columns as for
  bed files (see
  [`?gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md))
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  copies the backing file of bigSNP objects (never modifying them), and
  works on read-only files
- `rbind(use_position = TRUE)` keeps the loci names of the reference,
  and needs much less memory; integer names are matched by position
- the `as_is` argument of [`rbind()`](https://rdrr.io/r/base/cbind.html)
  is deprecated (it was never implemented)
- [`cbind()`](https://rdrr.io/r/base/cbind.html) renames duplicated
  columns as
  [`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html);
  `id` and `genotypes` can not be added
- `distruct_colours` holds the 56 distinct distruct colours, recycled
  with a warning by
  [`scale_fill_distruct()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/scale_fill_distruct.md)
- [`gt_as_genlight()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_genlight.md)
  keeps alleles, chromosomes and positions, so that
  [`gt_from_genlight()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_from_genlight.md)
  can convert back
- classes are preserved by all dplyr verbs; grouped functions keep the
  order of factor levels, and work with several grouping variables and
  empty groups
- [`heatmap_add_dendro()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/heatmap_add_dendro.md)
  returns a `gtable`, which can be combined with patchwork
- [`augment_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/augment_loci.md)
  works on
  [`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md)
  results
- [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
  uses ADMIXTURE on the PATH when conda is not installed
- faster “cpp” VCF parser, and
  [`select_loci_if()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci_if.md)
  no longer loads all the genotypes
- informative errors for invalid inputs to
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  `show_loci<-`,
  [`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_order_loci.md),
  [`gt_pcadapt()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pcadapt.md),
  [`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_dapc.md),
  [`pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_fst.md)
  and others

### Bug fixes

- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md):
  fix backing file names (truncated at dots, or “.vcf.bk”), and backing
  files left behind on errors
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md):
  fix files with trailing blank lines, CRLF line endings or upper case
  extensions
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  on bed files keeps .fam IDs as text, and quantitative phenotypes as
  numbers
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  on VCF files: fix the ploidy of individuals missing on the first
  variant, partially missing genotypes, and genotypes above the ploidy
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md):
  “0” and “.” are always missing alleles; fix the “vcfR” parser on
  indels and on missing genotypes on the second variant
- [`gt_load()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_load.md):
  fix the ploidies of gen_tibbles saved in the old bigSNP format
- [`gt_as_plink()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_plink.md):
  fix the .fam and .map information, missing alleles, and ped/raw files
  with numeric phenotypes, no population or pseudohaploids
- [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md)
  and
  [`gt_as_plink()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_plink.md):
  fix positions in scientific notation, and the date in the VCF header
- [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
  and
  [`gt_snmf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_snmf.md):
  fix seeds with several runs, paths containing “geno”, and the order of
  loci in P matrices
- [`gt_snmf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_snmf.md):
  fix continued projects, and the order of cross-entropies
- [`gt_impute_xgboost()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_xgboost.md):
  fix unsorted loci;
  [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_simple.md):
  restore the BLAS option
- [`filter_high_relatedness()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/filter_high_relatedness.md):
  fix missing and strongly negative kinship coefficients
- [`snp_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/snp_allele_sharing.md):
  fix blocks of different sizes;
  [`indiv_het_obs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_het_obs.md):
  fix the number of loci after
  [`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md)
- [`pop_global_stats()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_global_stats.md):
  fix loci where a population has no data;
  [`pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_tajimas_d.md):
  fix loci with no valid genotypes
- [`loci_alt_freq()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_alt_freq.md)
  and
  [`loci_maf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_alt_freq.md):
  fix polyploids and single individuals
- grouped `loci_*` functions and `qc_report_*`: fix factor, numeric,
  non-contiguous and empty groups, and p-values above 1
- [`loci_transitions()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_transitions.md)
  and
  [`loci_transversions()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_transversions.md):
  fix misclassified substitutions
- `windows_*` functions: fix the last loci of chromosomes;
  [`windows_pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_pairwise_pop_fst.md)
  and
  [`windows_nwise_pop_pbs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_nwise_pop_pbs.md)
  ignoring `method`
- [`gt_roh_window()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_window.md):
  fix homozygotes for the alternate allele with `heterozygosity = TRUE`
- [`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_dapc.md)
  and
  [`gt_cluster_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_cluster_pca.md):
  fix the default number of components, the variance retained, and the
  choice of k
- [`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md):
  fix unsorted loci;
  [`augment()`](https://generics.r-lib.org/reference/augment.html) of
  `gt_pca` and `gt_dapc`: fix single components and integer ids
- [`rbind()`](https://rdrr.io/r/base/cbind.html): fix ploidies after
  filtering, and missing `allele_ref`;
  [`rbind_dry_run()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind_dry_run.md):
  fix paths to .bim files
- dplyr verbs on grouped gen_tibbles: fix
  [`select()`](https://dplyr.tidyverse.org/reference/select.html),
  [`rename()`](https://dplyr.tidyverse.org/reference/rename.html),
  [`relocate()`](https://dplyr.tidyverse.org/reference/relocate.html),
  `[`,
  [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html),
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
  and [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  of `q_matrix` and `gt_admix`: fix the sorting options, numeric and
  unused groups, and `type = "cv"`; ordering in other
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)s
- [`gt_admix_reorder_q()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admix_reorder_q.md):
  fix K = 1
- [`heatmap_add_dendro()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/heatmap_add_dendro.md):
  fix ties and older ggplot2;
  [`heatmap_pairwise()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/heatmap_pairwise.md):
  fix factor ids
- fix chromosome names being renumbered in lexicographic order
- [`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pseudohaploid.md):
  fix other individuals of the backing file, and error on data that are
  not diploid
- [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html) of
  `q_matrix`: fix ignoring individual ids
- [`gt_as_genind()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_genind.md):
  fix loci names with “.”;
  [`gt_from_genlight()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_from_genlight.md):
  fix objects without `pop`;
  [`gt_as_hierfstat()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_hierfstat.md):
  fix numeric labels
- [`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md):
  fix numbers of loci that are multiples of `chunk_size`; fix printing
  gen_tibbles with no loci
- fix documentation and error messages of several functions

## tidypopgen 0.4.4

CRAN release: 2026-05-19

- fix seed argument of
  [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
- implement our own upset plot function to avoid dependency on `UpSetR`
  package, which is no longer maintained
- fix reading of VCFs when genotype separator symbol is present in the
  sample information
- fix `loci_` functions output when calculating across groups for a
  single loci
- fix
  [`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md)
  to preserve imputation metadata when backing files are updated
- add validation in
  [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_simple.md)
  to detect and report SNP dimension mismatches caused by excluded SNPs
- add `outdir` argument in
  [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
  to store partial outputs for large jobs
- fix
  [`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
  output file names when n_runs \> 1

## tidypopgen 0.4.3

CRAN release: 2026-01-23

- fix occasional test failure when temporary file names contained
  certain patterns

## tidypopgen 0.4.2

CRAN release: 2026-01-22

- allow least squares projection of multiple components at once in
  `gt_project_pca()`

## tidypopgen 0.4.1

CRAN release: 2025-12-20

- duplicate-distance checks removed in
  [`is_loci_table_ordered()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/is_loci_table_ordered.md)
- make `qc_report_*` functions work with pseudohaploids

## tidypopgen 0.4.0

CRAN release: 2025-10-24

- store chromosomes as factors, and remove redundant column storing
  chromosome as integer

## tidypopgen 0.3.3

- only use the FBM backend (rather than a complete bigsnpr object)
- fix a few minor bugs

## tidypopgen 0.3.2

CRAN release: 2025-08-27

- tidy up for CRAN submission

## tidypopgen 0.3.0

- add
  [`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pseudohaploid.md)
  and adapt appropriate functions to work with pseudohaploids

## tidypopgen 0.2.0

- Lots of optimisations to speed up processing of large datasets.
- Added an `sf` backend to `tidypopgen` for spatial data.
- Windowed statistics.

## tidypopgen 0.1.0

- Initial public version of `tidypopgen` released on r-universe.
