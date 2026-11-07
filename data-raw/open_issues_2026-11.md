# tidypopgen: issues still to solve (November 2026)

These are the findings of the October 2026 code review that were confirmed
but not fixed (the "Findings not fixed" table of
[code_review_2026-10.md](code_review_2026-10.md)), checked against `dev` in
November 2026. Numbers are those of that table. Line numbers are for `dev`
at the commit that adds this file.

## Still to solve

| # | Location | Issue | Suggestion |
|---|---|---|---|
| 5 | `R/rbind_gen_tibble.R:121-125` | `as_is` is not implemented: loci are always harmonised, although the message printed with `as_is = TRUE` says they will not be. | Decide what `as_is` should do when the loci differ (e.g. require identical loci tables, and error otherwise), and fix the message at the same time. |
| 6 | `R/filter_high_relatedness.R:140` | `mn2 <- mean(matrix2[-j, ])` (copied from `caret::findCorrelation()`) is not "the average relatedness of individual j", as the comment says. | Both versions are greedy heuristics that leave no pair above the threshold, and neither was clearly better in simulations. Either fix the comment, or use `mean(matrix2[j, ])` and document the change. |
| 11 | `R/autoplot_gt_dapc.R:122` (also `autoplot.gt_pca(type = "loadings", k = c(1, 2))`, via bigstatsr) | `geom_hex()` needs `hexbin`, which is not declared in `DESCRIPTION`. | Add `hexbin` to Suggests and call `rlang::check_installed("hexbin")`, or plot points. |
| 12 | `data/distruct_colours` | The palette repeats colours (positions 4/35, 12/54, 24/46, 59/60), so with K >= 35 some clusters share a colour, and K > 60 errors. | Replace the duplicated colours, and give a clear error (or recycle with a warning) for K > 60. |
| 13 | VCF parsers (`src/vcf_parser_fast.cpp:202`) | Ploidy is taken from the first non-missing genotype, so haploid calls (e.g. chrX in males) on the first variant set ploidy 1. The cpp parser only looks for that genotype in the first 1,000 variants; individuals missing on all of them get the largest ploidy of the others. | Documented limitation, but the docs of `gen_tibble()` (`R/gen_tibble.R:6-14`) say the ploidy comes from the first genotype that is not '.', without the 1,000 variants limit: update them. |
| 14 | `R/gen_tibble.R:529` (`new_vctrs_bigsnp()`) | `fbm_md5sum` is NA for VCF inputs, as the `.rds` is hashed before it is saved. | The attribute is not used anywhere: either drop it, or save the `.rds` before hashing it. |
| 15 | `R/gt_pca_autoSVD.R:141` | `infos_chr + abs(min(infos_chr) + 1)` is wrong when the smallest code is -1 or less. | No observable effect (bigsnpr only uses `infos.chr` to group SNPs), but `infos_chr - min(infos_chr) + 1` is what was meant. |

## Open from the November 2026 speed and ROH work

These come from the work on `speed_extra` and `roh` (see
[speed_review_2026-11.md](speed_review_2026-11.md)). They are decisions or
ideas, not bugs.

| Topic | Question or idea |
|---|---|
| `ignore_genetic_dist` default | `is_loci_table_ordered()` and `gt_order_loci()` check positions only by default (`ignore_genetic_dist = TRUE`). Missing positions, and missing genetic distances when they are checked, are now informative errors (`7dad1c9f`). Should checking genetic distances be the default, with positions only as the opt-in? |
| `genetic_dist` column | It takes 8 bytes per locus (0.32 GB at 40M loci) and is all zeros for VCF inputs. It is needed for cM blocks and windows, `gt_extract_f2()` and PLINK exports, but could become optional (filled with 0 when exported). |
| Runs of homozygosity across chromosomes | `windows_indiv_roh()` keeps the algorithm of detectRUNS: windows slide along the whole loci table, so near a boundary they mix the last loci of one chromosome and the first of the next (runs never span two chromosomes). PLINK resets windows at each chromosome; an option to do so would knowingly depart from detectRUNS. |
| Duplicated loci names | `validate_loci()` now finds duplicated positions in one pass when loci are sorted (`0b3df790`), but duplicated character names still take `anyDuplicated()`: 3.2 s and 0.5 GB for 40M names (integer names are cheap). |
| AVX2 popcount | Further idea 11 of the speed review: compiling the pairwise row kernel for AVX2 was 15% slower; a hand-written Harley-Seal popcount with run-time dispatch remains possible, with CRAN portability risks. |
| `read_pedfile()` | The genotype matrix is still held in memory (now as integers, never copied); streaming it to the FBM would need the allele coding of each locus, known only once all individuals are read. |
| `NEWS.md` | The speed commits on `speed_extra` have no `NEWS.md` entries (nothing visible changes); decide whether to add a general entry. |

### Not reproduced

- **#9** (`R/cbind_gen_tibble.R`, "`cbind()` on a grouped gen_tibble returns
  an invalid grouped object"): on `dev`, `cbind(grouped_gt, df)` returns a
  valid grouped gen_tibble (`dplyr:::validate_grouped_df()` passes, and
  `summarise()` works). Two related cases do not give an invalid grouped
  object either, but may be worth a look:
  - `cbind(grouped_gt, data.frame(population = ...))` (a column with the name
    of the grouping variable) renames both columns (`population...2` and
    `population...4`) and returns an ungrouped gen_tibble.
  - `cbind(df, gt)` (the gen_tibble not first) is dispatched to
    `cbind.data.frame()`, and returns a plain data.frame.

### Minor points found while checking

- `R/snp_king.R`: the KING matrix is stored in a variable called
  `ibs_counts_matrix`; rename it.
- `pairwise_ibs(type = "adjusted_counts")` matches PLINK's
  `--distance flat-missing`. PLINK's default correction for missing genotypes
  (without `flat-missing`) is not available; add it only if users ask for it.

## Resolved since the review

| # | Issue | Resolution |
|---|---|---|
| 1 | `gen_tibble()` on a bigSNP `.rds` overwrote that `.rds` and ignored `backingfile`. | Fixed in `2f12715a` (test in `c8e36257`): the `.bk` file is copied to `backingfile`. This doubles disk use for the bigSNP. |
| 2 | `pairwise_ibs(type = "adjusted_counts")` was a similarity in loci, documented as PLINK's "distance". | Fixed in the commit that adds this file: it is now PLINK's allele count distance (`--distance flat-missing`), tested against a PLINK 1.9 `.dist` file (`inst/extdata/related/test_plinkDist.dist`). No other function used it. |
| 3 | The pairwise kernels treated imputed genotypes as missing. | Fixed in `f26f046b`: genotypes are decoded through `code256`. |
| 4 | WC84 drops loci with no data in one population, whereas `hierfstat::wc()` keeps them with sigma_a = 0. | Kept. The difference is documented in the details of `pairwise_pop_fst()`, and tested, in the commit that adds this file. Nei87 already matches `hierfstat::pairwise.neifst()`. |
| 7 | `gt_extract_f2()` counted groups with `max(group_indices())`. | Fixed on `dev` (#218): groups are counted with `.group_info()`, and empty groups are dropped. |
| 8 | `pop_fst()` and `pop_fis()` counted groups with `factor(group_indices())`. | Fixed on `dev` (#218): empty groups give NaN. |
| 10 | `gt_as_hierfstat()` without a `population` column failed with an unclear "differing number of rows" error. | Fixed in `8949cf81` (on `speed_extra`, merged into `roh`): it stops with a clear error. (`c25bf1e2`, which builds the data frame column by column, had made it return a corrupt data frame instead.) |
