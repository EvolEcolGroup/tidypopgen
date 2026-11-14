# Changelog

## tidypopgen 1.1.0

### New features and changes

- new
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
  and
  [`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md)
  read several VCF files into a single `gen_tibble`, as `bcftools merge`
  and `bcftools concat`
- integer chromosomes (`chromosomes_as_int`) for datasets with many
  chromosomes; the original names of integer loci and chromosomes are
  kept in sidecar files
  ([`gt_lookup_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md),
  [`gt_restore_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  and
  [`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)),
  which
  [`gt_get_file_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_get_file_names.md)
  and
  [`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md)
  list
- phased genotypes (experimental): `gen_tibble(phased = TRUE)` reads
  phased VCF files;
  [`is_phased()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/is_phased.md)
  and
  [`show_haplotypes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_haplotypes.md)
  give the phase,
  [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md)
  writes it, and [`rbind()`](https://rdrr.io/r/base/cbind.html) keeps it
- `pairwise_ibs(missing_correction = "allele_freq")` corrects for
  missing genotypes as PLINK 1.9’s `--distance` without `flat-missing`
- `qc_report_loci(hwe = FALSE)` gives the MAF and missingness of loci in
  a single pass
- faster VCF reading (in parallel with `n_cores`), and faster reading of
  *.ped* files straight into the backing file
- faster
  [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md),
  [`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md)
  and
  [`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_order_loci.md)
  with many chromosomes;
  [`gt_as_genlight()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_genlight.md)
  works in chunks of loci
- faster
  [`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md)
  with variables and function calls; `.engine = "tidyselect"` always
  uses `tidyselect`
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  and
  [`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md)
  save the *.gt* file twice as fast
- [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_simple.md)
  works after removing loci (but not individuals)
- new article on large RAD-seq datasets; the big data article uses
  [`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
- new
  [`bootstrap_blocks()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/bootstrap_blocks.md),
  for developers (not in the reference index): the blocks of loci of the
  block bootstrap and the blocks drawn by each replicate, so that other
  packages can resample the same blocks as the bootstrap functions of
  tidypopgen
- new
  [`pairwise_dist()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_dist.md):
  genetic distances between individuals, as ngsDist, with the genotype
  or nucleotide score, corrected for saturation (`-log(1 - d)` by
  default, Jukes-Cantor, or none), e.g. to build neighbour-joining
  trees; an S3 generic, so that other packages (e.g. for genotype
  likelihoods) can add methods
- [`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md),
  [`pairwise_ibs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_ibs.md),
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
  (with methods for `gen_tibble` objects and their `genotypes` column)
  and
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)
  (with a method for `gen_tibble` objects) are now S3 generics, and the
  generic
  [`indiv_inbreeding()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_inbreeding.md)
  only has the arguments `.x` and `...` (its method for `gen_tibble`
  objects is a `gen_tbl` method), so that other packages can provide
  methods for their own objects, with their own arguments.
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)
  now also works on `gen_tibble` objects that lost their class
  (e.g. with
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)),
  as the other generics do

### Bug fixes

- functions using temporary backing files (PCA, LD clumping, KING, IBS
  and allele sharing) delete them, fixing occasional “File ‘…bk’ already
  exists” errors
- [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)
  and the `gt_roh_*()` functions give chromosome names with
  `chromosomes_as_int`
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  on *.ped* files: a third allele always flags the locus as
  multi-allelic
- [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  with the “cpp” parser: fix malformed VCF files read as wrong genotypes
  (now an error), missing chromosomes (now `NA`), and backing files left
  behind on errors
- [`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md):
  fix loci with missing chromosomes (now written as “.”)

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
