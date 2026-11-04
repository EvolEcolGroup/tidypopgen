# tidypopgen `dev`
* new functions to compute site frequency spectra (SFS) for demographic
  inference: `pop_sfs()` (1D spectrum of each population),
  `pairwise_pop_sfs()` (joint spectrum of each pair of populations) and
  `nwise_pop_sfs()` (joint spectrum of all populations). Missing genotypes
  are handled as in easySFS (and dadi), by projecting each population down to
  a fixed number of allele copies. `pop_sfs_preview()` gives the number of
  segregating sites for each projection (as `easySFS --preview`) and the
  projection that maximises it, which is used by default. Spectra are
  returned as tidy tibbles (with the order of the cells of joint spectra set
  by `fastest`), or as `sfs` arrays (`type = "array"`). `write_sfs()` writes
  a spectrum in dadi or fastsimcoal2 format, and `gt_as_easysfs()` writes all
  the spectra written by easySFS (1D, 2D and multidimensional, in both
  formats, with the same file names), optionally setting the number of
  monomorphic sites from the total length of the sequences. Unfolded spectra
  can be polarised with the ancestral allele of each locus (`ancestral`, a
  column of the loci table or a vector); loci with an unknown ancestral
  allele are then dropped, as in dadi.
* grouped functions now work with gen_tibbles grouped by more than one
  variable (e.g. `group_by(population, sex)`), including `loci_alt_freq()`,
  `loci_maf()`, `loci_missingness()`, `loci_hwe()`, `loci_pi()`,
  `qc_report_indiv()`, `nwise_pop_pbs()`, the windowed and jackknife versions
  of `pop_tajimas_d()` and `nwise_pop_pbs()`, `gt_admixture()`, `gt_snmf()`
  and `pairwise_pop_fst(type = "pairwise")`, which used to stop with an
  error. Groups are labelled by the values of the grouping variables pasted
  with "." (e.g. "pop1.F"), and "tidy" outputs have one column per grouping
  variable. This also fixes silently wrong labels with several grouping
  variables: the names of the results of `pop_het_obs()`, `pop_het_exp()`,
  `pop_fis()` and `pop_fst()`, the family ID written by `gt_as_plink()`, and
  the groups of `windows_indiv_roh()` and of `tidy()`/`augment()` of a
  `q_matrix`, which only used the first grouping variable; and the Bonferroni
  correction of `qc_report_loci()`, which counted the values of the first
  grouping variable rather than the groups.
  `heatmap_pairwise()` and `autoplot()` accept the tidy output of
  `pairwise_pop_fst()` with several grouping variables
* empty groups, which only exist if requested (with `group_by(.drop = FALSE)`,
  or `filter()`/`slice()` with `.preserve = TRUE`), are kept in the results
  as NA/NaN, without changing the results of the other groups. `pop_fst()` and
  `pop_fis(method = "WG17")` used to fail with an empty group, and
  `qc_report_indiv(kings_threshold = )` now skips empty groups
* `gt_extract_f2()` always drops empty groups (from `group_by(.drop = FALSE)`,
  or `filter()`/`slice()` with `.preserve = TRUE`), rather than passing them to
  *ADMIXTOOLS 2* as populations without data
* grouped `loci_hwe()` now returns NA for empty groups, rather than the p-value
  of a test on no genotypes (0.5 with `mid_p = TRUE`); grouped
  `qc_report_loci()` ignores empty groups accordingly, as they could lower the
  smallest HWE p-value to 0.5
* functions applied to a `gen_tibble` without individuals (as `group_map()`
  does for an empty group) now return NA for each locus (`loci_alt_freq()`,
  `loci_maf()` and `loci_hwe()`, and NaN for `loci_missingness()`, which
  returned an empty vector), or empty results (`pairwise_king()`,
  `pairwise_ibs()`, `pairwise_allele_sharing()`, `pairwise_grm()` and
  `indiv_inbreeding()`), rather than failing, as the grouped methods do for
  empty groups. Grouped `indiv_inbreeding()`, which uses `group_map()`, now
  works with empty groups
* fix `pairwise_allele_sharing()`, `pairwise_grm()` and `indiv_inbreeding()`
  failing ("'dimnames' applied to non-array") with a single individual, and
  thus grouped `indiv_inbreeding()` with a group of one individual;
  `indiv_inbreeding()` is NaN for a single individual, as it is relative to the
  other individuals
* fix `windows_pop_tajimas_d(type = "matrix")` splitting a window over several
  rows when the number of loci with data in the window differs between groups
  (e.g. with an empty group, or a group without data at a locus); `n_loci` is
  now the number of loci in the window, and the number with data in each group
  is given by the "tidy" and "list" outputs
* the `group` column of the "tidy" outputs of grouped `loci_alt_freq()`,
  `loci_maf()`, `loci_missingness()`, `loci_hwe()`, `loci_pi()`,
  `windows_pop_tajimas_d()` and `jackknife_pop_tajimas_d()` now keeps the type
  of the grouping variable (e.g. a factor keeps its levels, and numbers stay
  numbers), as `dplyr::summarise()` does and as the columns of several
  grouping variables already did, rather than being character; groups thus
  stay in their order (e.g. the levels of the factor) when the output is
  sorted or plotted
* with several grouping variables, the group of each individual (the `group`
  column of `qc_report_indiv()`, and the `group` element of the results of
  `gt_admixture()` and `gt_snmf()`) is now a factor whose levels are the
  group labels in the order of the groups, rather than characters, which
  sorted alphabetically; `gt_admix_reorder_q()` and `autoplot()` thus order
  the groups as the grouped functions do
* fix `tidy()`, `augment()` and `autoplot()` of a `q_matrix` with grouped
  `data` turning the groups into characters, so that `autoplot()` ordered them
  alphabetically rather than by the levels of a factor (as documented) or in
  numeric order; and ignoring the groups of `data` when the `q_matrix` had the
  ids of the individuals but no groups (e.g. from `gt_admixture()` or
  `gt_snmf()` on an ungrouped `gen_tibble`), so that `autoplot()` could not
  annotate the groups
* fix `autoplot()` of `q_matrix` and `gt_admix` objects with
  `reorder_within_groups = TRUE` sorting the groups even with
  `arrange_by_group = FALSE`
* fix `heatmap_pairwise()` and `autoplot()` of tidy pairwise outputs sorting
  the individuals or groups alphabetically, rather than keeping their order (as
  documented, and as for the same data as a matrix)
* fix `autoplot()` of `q_matrix` and `gt_admix` objects failing ("`breaks` and
  `labels` must have the same length") when the first group has a single
  individual, and putting the label of a group of a single individual on the
  last individual of the previous group
* fix `gen_tibble()` on a bigSNP .rds file overwriting that .rds and ignoring
  `backingfile`; the backing file is now copied to `backingfile`, so that the
  files of the bigSNP are never modified (`gt_impute_simple()` and
  `gt_impute_xgboost()` used to write into its .bk file, and `gen_tibble()`
  deleted both files when it stopped on loci with no alleles)
* `gen_tibble()` now works on read-only bigSNP .rds and .bk files (it used to
  fail when overwriting the .rds); the copy of the backing file is writable,
  so that the gen_tibble can be imputed
* fix `filter_high_relatedness()` (and `qc_report_indiv()` with
  `kings_threshold`) failing when some kinship values are missing (e.g. pairs of
  individuals with no loci genotyped in both); such pairs are treated as not
  related, with a warning
* `gen_tibble()` now errors on truncated or corrupted packedancestry .geno
  files, rather than reading garbage genotypes (or reading out of bounds)
* document that the `n_cores` argument of `gen_tibble()` is currently ignored
* speed up reading genotypes with the "cpp" VCF parser of `gen_tibble()`
* fix `gen_tibble()` truncating user-supplied `backingfile` names at the last
  dot (e.g. "pop_a.filtered" became "pop_a"), and naming the backing files of
  "x.vcf.gz" "x.vcf.bk" rather than "x.bk"
* fix `gen_tibble()` failing on VCF and PED files with a trailing blank line,
  and keeping a carriage return in the last sample name of VCF files with CRLF
  line endings
* fix `gen_tibble()` failing on files with upper case extensions (e.g. .BED,
  .PED, .GENO, .VCF.GZ)
* fix `gen_tibble()` on PLINK bed files turning the family, individual and
  parental IDs of the .fam file into numbers when they look like numbers (e.g.
  "001" became "1", and "1.10" became "1.1"); they are now kept as they are, as
  when reading ped files
* fix `gen_tibble()` on PLINK bed files setting quantitative phenotypes (with
  values other than 1, 2, 0 and -9 in the .fam file) to NA; as in PLINK, they
  are now kept as numbers, with -9 as missing
* fix `gen_tibble()` on VCF files leaving a backing file behind when the loci
  fail validation (e.g. duplicated loci), which caused a retry to use a new
  versioned file name
* fix `gen_tibble()` on VCF files setting the ploidy of individuals with a
  missing genotype (".") on the first variant to 1; their ploidy is now taken
  from their first non-missing genotype (searched in the first 1,000 variants
  by the "cpp" parser), or is the largest ploidy of the other individuals if
  there is none
* fix `gen_tibble()` on VCF files silently corrupting genotypes with more
  alleles than the ploidy estimated from the first variant (e.g. a tetraploid
  genotype in a diploid file) when they are not on the last line; this is now an
  error
* fix `gen_tibble()` rejecting "0" alleles (the PLINK missing value) when custom
  `missing_alleles` are given
* fix `gen_tibble()` on ped files reading a missing genotype "0 0" (the PLINK
  default, also written by `gt_as_plink()`) as an allele "0" when custom
  `missing_alleles` are given, which silently corrupted the genotypes of the
  locus (or stopped with an error for loci with only missing genotypes); "0"
  is now always missing in ped files
* fix `gen_tibble()` on VCF files with the "cpp" parser stopping with an error
  on missing alleles "." (e.g. the ALT of monomorphic loci) when custom
  `missing_alleles` are given; as with the "vcfR" parser, "." is now always
  missing
* fix `gt_as_plink(type = "raw")` writing a missing `allele_ref` as "NA"
  rather than "0" in the header of raw files; missing alleles are now "0" in
  both allele columns, as in bim files
* `gen_tibble()` on a matrix now checks that a vector of `ploidy` has one value
  per individual, and that dosages are non-negative integers not greater than
  each individual's ploidy (previously invalid values were silently corrupted)
* fix `gt_load()` on gen_tibbles saved in the old bigSNP format losing the
  individual ploidies of mixed-ploidy and pseudohaploid data, and setting too
  few ploidies for subsetted gen_tibbles
* fix `gt_admixture()` with `seed = NULL` (the default) giving identical runs
  when `n_runs > 1`, as ADMIXTURE used its default seed (43) for every run;
  runs now use seeds 43, 44, etc.
* fix `gt_snmf()` and `gt_admixture()` failing when the path of the input file
  contains "geno" (or the file name contains ".bed") other than in the extension
* fix `gt_snmf()` failing with the default `seed` when `n_runs > 1`
* fix `gt_snmf()` returning the results of earlier runs when adding runs to an
  existing project (`project = "continue"` or `"force"`)
* fix `gt_snmf()` returning cross-entropies (`$cv`) in the wrong order when
  `n_runs > 1` and several values of `k` are used
* `gt_pcadapt()` now errors if the individuals of `x` do not match (in identity
  and order) those used for `pca`
* fix `gt_impute_xgboost()` assigning loci to the wrong chromosome (and so
  imputing them from the wrong neighbours) when the loci table is not in the
  order of the backing file
* fix `gt_impute_simple()` with `n_cores > 1` resetting the option
  `bigstatsr.check.parallel.blas` to `TRUE` rather than to its previous value
* fix `filter_high_relatedness()` (and the `kings_threshold` option of
  `qc_report_indiv()`) treating strongly negative kinship coefficients as high
  relatedness
* fix error messages and documentation of `autoplot()` for `gt_admix` and
  `gt_pca` objects referring to non-existent arguments (`repeat`, `pc`)
* fix `autoplot()` of `gt_admix` objects with `type = "cv"` failing when the
  object has no `algorithm` element
* fix `gt_admix_reorder_q()` turning the Q matrices for K = 1 into vectors
* fix `reorder_within_groups` in `autoplot()` for `q_matrix` and `gt_admix`
  objects, which did not sort individuals by their dominant component
* fix `arrange_by_indiv` in `autoplot()` for `q_matrix` and `gt_admix`
  objects, which did not sort individuals by their id within each group
* fix `autoplot()` for `q_matrix` and `gt_admix` objects with numeric groups
  (e.g. from `gt_admixture()` or `gt_snmf()` on a `gen_tibble` grouped by a
  numeric column) drawing no group separators or labels
* fix group annotations of admixture barplots showing extra separators and
  labels for unused levels of the grouping factor
* `select_loci_if()` no longer loads the whole genotype matrix into memory to
  check the length of the selection
* fix `gt_as_genind()` failing when loci names contain '.' (which are now
  replaced by '_', with a warning)
* `gt_as_genlight()` now stores the alleles (`loc.all`), chromosome and position
  of each locus, so that the result can be converted back with
  `gt_from_genlight()`
* fix `gt_from_genlight()` rejecting `genlight` objects without a `pop` slot,
  which is documented as optional
* fix `rbind_dry_run()` failing when `ref` or `target` are paths to PLINK .bim
  files, as documented
* `rbind(use_position = TRUE)` now keeps the loci names of the reference
  dataset, instead of renaming loci to "chromosome_position" and adding a
  `name_old` column
* fix `rbind()` silently accepting `as_is = TRUE` with `flip_strand = TRUE` when
  `quiet = TRUE`
* fix `select()`, `rename()` and `relocate()` on a grouped `gen_tibble`
  returning a plain `grouped_df`
* fix `[` on a grouped `gen_tibble` returning a corrupt grouped object when the
  grouping variables are dropped
* fix `ungroup()` dropping the `sf` class and geometry column attribute of
  grouped `gen_tibble`s, and returning malformed classes when only some grouping
  variables are removed
* fix `group_by()` without grouping variables returning a corrupt grouped
  `gen_tibble`
* fix `mutate()` on an ungrouped `gen_tibble` returning a corrupt `grouped_df`
  when the genotypes column is dropped
* fix `windows_indiv_roh()` with `heterozygosity = TRUE` treating homozygotes
  for the alternate allele as heterozygotes
* fix `windows_stats_generic()` (and thus all `windows_*` functions) leaving out
  the last loci of a chromosome when the window range is not a multiple of
  `step_size`
* fix spurious warnings from `windows_pairwise_pop_fst()` and
  `windows_nwise_pop_pbs()` with `type = "tidy"` for some numbers of populations
* fix `gt_as_vcf()` writing a space-padded day in the `##fileDate` header line
  on the first nine days of a month
* fix `gen_tibble()` with `parser = "vcfR"` failing on VCF files containing
  biallelic indels; these are now dropped, as with the "cpp" parser
* `show_loci<-` now errors if the replacement table does not have a valid
  `big_index` column, rather than silently creating a broken `gen_tibble`
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
* `gen_tibble()` on PLINK ped files now stores the individual information in
  the same columns, with the same types, as for bed files (see "Details" in
  `?gen_tibble`): the parents are `paternal_ID` and `maternal_ID` (rather than
  `father` and `mother`); `sex` is a factor with levels "female" and "male"
  (rather than PLINK's codes 2 and 1); `phenotype` is a factor with levels
  "control" and "case" (rather than the codes 1 and 2, with -9 for missing
  values), or a number for a quantitative phenotype (with NA rather than -9
  for missing values); `population` is only created if the family IDs differ
  from the individual IDs, and parents, sex and phenotype only if they have
  some non-missing values. `gt_as_plink()` thus writes back the parents and
  sex read from ped files, which used to be written as unknown. `gen_tibble`s
  read from ped files by earlier versions keep the old columns: rename and
  recode them (e.g. with
  `rename(paternal_ID = father, maternal_ID = mother)`) to write them with
  `gt_as_plink()`
* fix `gt_as_plink(type = "bed")` not writing the population (as family ID),
  parents, sex and phenotype to the .fam file, and `gt_as_plink(type =
  "ped")`/`type = "raw"` ignoring the grouping variable as family ID
* fix `gt_as_plink(type = "ped")` and `gt_as_plink(type = "raw")` dropping the
  parents (`paternal_ID` and `maternal_ID`, as created when reading PLINK
  files)
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
* `gt_admixture()` with `conda_env = "auto"` (the default) uses a copy of
  ADMIXTURE on the PATH when conda is not installed, as documented: it
  failed with the error of `reticulate::conda_list()` ("Unable to find
  conda binary"). The check for ADMIXTURE on the PATH also works on Windows.
* fix `rbind()` and `rbind_dry_run()` failing ("NAs are not allowed in
  subscripted assignments") on loci with a missing `allele_ref`. A missing
  allele is now inferred from the other dataset whichever column it is in
  (previously only `allele_alt` was), and the report of `rbind_dry_run()`
  gains a `missing_allele_col` column giving the column that was filled. With
  `flip_strand = TRUE`, loci whose `allele_ref` is missing in both datasets
  are swapped so that the observed allele becomes `allele_ref`; these are
  flagged in a new `to_swap` column of the `ref` report.
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
* add `jackknife_*` functions (`jackknife_pairwise_pop_fst()`,
  `jackknife_nwise_pop_pbs()`, `jackknife_pop_tajimas_d()` and
  `jackknife_stats_generic()`) to estimate genome-wide statistics with
  standard errors from a weighted block jackknife, with blocks defined in
  SNPs, centiMorgans or base pairs (blocks in base pairs cover the same
  positions as the windows of the `windows_*` functions with `step_size`
  equal to `window_size`).
* `pop_tajimas_d()` now works with polyploid and pseudohaploid data.
* `pop_tajimas_d()`, `windows_pop_tajimas_d()` and `jackknife_pop_tajimas_d()`
  now account for missing genotypes by using the number of non-missing allele
  copies at each locus as its sample size, rather than the total number of
  allele copies in the sample. This follows Bailey et al. (2025), and results
  match those of `pixy` (version 2.2.4). Results will differ from previous
  versions when data contain missing genotypes.

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
