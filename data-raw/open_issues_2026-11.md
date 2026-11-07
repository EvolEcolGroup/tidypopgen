# tidypopgen: issues still to solve (November 2026)

These are the findings of the October 2026 code review that were confirmed
but not fixed (the "Findings not fixed" table of
[code_review_2026-10.md](code_review_2026-10.md)), checked against `dev` in
November 2026. Numbers are those of that table. Line numbers are for `dev`
at the commit that adds this file.

## Still to solve

All the confirmed findings of the October 2026 review are now resolved (see
[Resolved since the review](#resolved-since-the-review)); the points below
remain.

### Open from the November 2026 speed and ROH work

These come from the work on `speed_extra` and `roh` (see
[speed_review_2026-11.md](speed_review_2026-11.md)). They are decisions or
ideas, not bugs.

| Topic | Question or idea |
|---|---|
| `genetic_dist` column | It takes 8 bytes per locus (0.32 GB at 40M loci) and is all zeros for VCF inputs. It is needed for cM blocks and windows, `gt_extract_f2()` and PLINK exports, but could become optional (filled with 0 when exported). |
| Duplicated loci names | `validate_loci()` now finds duplicated positions in one pass when loci are sorted (`0b3df790`), but duplicated character names still take `anyDuplicated()`: 3.2 s and 0.5 GB for 40M names (integer names are cheap). |
| AVX2 popcount | Further idea 11 of the speed review: compiling the pairwise row kernel for AVX2 was 15% slower; a hand-written Harley-Seal popcount with run-time dispatch remains possible, with CRAN portability risks. |
| `read_pedfile()` | The genotype matrix is still held in memory (now as integers, never copied); streaming it to the FBM would need the allele coding of each locus, known only once all individuals are read. |

Resolved: windows of `windows_indiv_roh()` no longer cross chromosomes. They
restart at each chromosome, as in PLINK, so the runs are those of detectRUNS
run on each chromosome separately (`0248be2e`).

### Not reproduced

- **#9** (`R/cbind_gen_tibble.R`, "`cbind()` on a grouped gen_tibble returns
  an invalid grouped object"): on `dev`, `cbind(grouped_gt, df)` returns a
  valid grouped gen_tibble (`dplyr:::validate_grouped_df()` passes, and
  `summarise()` works). Two related cases do not give an invalid grouped
  object either, but may be worth a look:
  - `cbind(grouped_gt, data.frame(population = ...))` (a column with the name
    of the grouping variable) renames both columns (`population...2` and
    `population...4`), with a message, and returns an ungrouped gen_tibble.
    This is the behaviour of `dplyr::bind_cols()`, which dplyr uses for
    `cbind()` on grouped data frames (`cbind.grouped_df()`), and is kept.
    An ungrouped gen_tibble used to get two columns with the same name; it
    is now renamed in the same way, and a column called `id` or `genotypes`
    is an error (`7627bfc5`).
  - `cbind(df, gt)` (the gen_tibble not first) is dispatched to
    `cbind.data.frame()`, and returns a plain data.frame.

### Minor points found while checking

- `R/snp_king.R`: the KING matrix was stored in a variable called
  `ibs_counts_matrix`; it is now `king_matrix` (`caa23187`).
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
| 5 | `rbind()` did not implement `as_is` (loci were always harmonised), although with `as_is = TRUE` it said that they would not be. | Deprecated in `dd527bc5`: `as_is = TRUE` gives a warning and has no effect, and the misleading message is removed. With a `NEWS.md` entry. |
| 6 | `filter_high_relatedness()`: the comment said that the individual with the higher mean relatedness is removed, but (as in `caret::findCorrelation()`) the mean of the first individual is compared with the mean of the remaining matrix without the second. | The comments and the details of the help now describe the greedy rule (`78826707`); the code is unchanged. |
| 7 | `gt_extract_f2()` counted groups with `max(group_indices())`. | Fixed on `dev` (#218): groups are counted with `.group_info()`, and empty groups are dropped. |
| 8 | `pop_fst()` and `pop_fis()` counted groups with `factor(group_indices())`. | Fixed on `dev` (#218): empty groups give NaN. |
| 10 | `gt_as_hierfstat()` without a `population` column failed with an unclear "differing number of rows" error. | Fixed in `8949cf81` (on `speed_extra`, merged into `roh`): it stops with a clear error. (`c25bf1e2`, which builds the data frame column by column, had made it return a corrupt data frame instead.) |
| 11 | `geom_hex()` (loadings of several axes in `autoplot()` of `gt_dapc` and `gt_pca`) needs `hexbin`, which was not declared. | `hexbin` is in Suggests, and both functions call `rlang::check_installed("hexbin")` before plotting several axes (`60cfe4b8`). |
| 12 | The distruct palette repeats 4 of its 60 colours, so clusters shared colours for K >= 35, and K > 60 failed. | `distruct_colours` holds the 56 distinct colours, and `scale_fill_distruct()` recycles them with a warning when there are more groups (`d4e156df`; needs ggplot2 >= 3.5.0). With a `NEWS.md` entry. |
| 13 | The documentation of `gen_tibble()` did not say that the cpp parser only looks for the ploidy of individuals missing on the first variant in the next 1,000 variants. | Documented (`6d387fed`), with the behaviour of the vcfR parser (the whole file) and the fallback (the largest ploidy of the others, or 2). Haploid calls on the first variant remain a documented limitation. |
| 14 | `fbm_md5sum` was NA for VCF inputs, as the `.rds` was hashed before it was saved. | Both VCF readers save the `.rds` before the gen_tibble is created, and `gen_tibble()` then saves only the gen_tibble for every input, so the FBM is no longer saved twice (`b2e614d9`). Saving and hashing the `.rds` takes about 60 ms whatever the size of the data. |
| 15 | `gt_pca_autoSVD()` shifted chromosome codes with `infos_chr + abs(min(infos_chr) + 1)`, wrong when the smallest code is -1 or less. | Fixed to `infos_chr - min(infos_chr) + 1` (`c22fa88d`); no observable effect, as codes are never negative. |
