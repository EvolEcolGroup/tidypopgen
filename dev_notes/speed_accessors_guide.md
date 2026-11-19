# Writing fast genotype kernels in tidypopgen

This guide explains how to use the raw-byte accessors added on the `speed`
branch (November 2026) when writing new functions that read genotypes from the
file-backed matrix (FBM). It complements
[speed_review_2026-11.md](speed_review_2026-11.md), which describes the
changes and their measured effect.

## Contents

- [Why not use bigstatsr's accessors in the inner loop](#why-not-use-bigstatsrs-accessors-in-the-inner-loop)
- [Which tool to use](#which-tool-to-use)
- [Principles that apply to every kernel](#principles-that-apply-to-every-kernel)
- [`fbm_byte_acc.h` reference](#fbm_byte_acch-reference)
- [Worked example: heterozygotes by locus and group](#worked-example-heterozygotes-by-locus-and-group)
- [Calling a kernel from R](#calling-a-kernel-from-r)
- [Per-individual statistics: `gt_ind_counts.h`](#per-individual-statistics-gt_ind_countsh)
- [Pairwise statistics: `pairwise_bitpacked.h`](#pairwise-statistics-pairwise_bitpackedh)
- [Copying and exporting genotypes](#copying-and-exporting-genotypes)
- [Checklist for a new kernel](#checklist-for-a-new-kernel)
- [Pitfalls](#pitfalls)

## Why not use bigstatsr's accessors in the inner loop

Genotypes are stored in an `FBM.code256`: one byte per genotype, individuals
in rows and loci in columns, column-major. The byte is decoded through the
`code256` vector (256 doubles): with `bigsnpr::CODE_012`, bytes 0, 1 and 2
decode to dosages and byte 3 to `NA`. After `gt_impute_simple()`, bytes 4 to 6
hold imputed calls, and `gt_set_imputed()` swaps `code256` so that the same
bytes decode either to `NA` or to the imputed value.

The original kernels used bigstatsr's `SubBMCode256Acc` and accumulated into
Rcpp matrices inside the OpenMP loop. On HGDP (1,043 × 643,733) the
allele-frequency kernel took 6.1 s, while just reading every byte once took
0.40 s. The kernels review traced the 15× gap to:

1. `ploidy(i)`-style accesses to Rcpp vectors are bounds-checked and call
   `Rf_xlength` for every genotype (about 13% of the time).
2. Adding into an Rcpp `NumericMatrix` inside the loop (strided writes, and
   possible aliasing that stops the compiler keeping values in registers)
   doubles the cost.
3. The two-level indexing (`rowInd` and `colInd`, copied to `size_t`
   vectors), the decode to `double` and the branch on missing values cost
   another 3.4×.

The new kernels read the raw bytes and turn `code256` into small integer
lookup tables at call time, so each genotype costs one load and one integer
addition. They are **not** a replacement for bigstatsr in general: use
bigstatsr/bigsnpr for linear algebra (`big_SVD`, `big_prodMat`,
`snp_clumping`, ...) and for anything that is not an inner loop over
genotypes.

## Which tool to use

| Task | Tool | Examples |
|---|---|---|
| A statistic per locus, over all individuals or by group | `src/fbm_byte_acc.h` (`FBMByteAcc`, `GroupedRows`, `lut_sums()` / `LutSums`) | `loci_alt_freq()`, `loci_missingness()`, `loci_pi()`, `loci_hwe()`, `qc_report_loci()` |
| Alternate and valid allele counts by group, with ploidy weights | `alt_valid_counts()` in `fbm_byte_acc.h` | `loci_alt_freq()`, grouped and polyploid |
| Population summaries (diploid/pseudohaploid) per block of loci, followed by a second phase (pairs, windows) | `GroupedSummariesDipPseudo` in `fbm_byte_acc.h` | `pop_global_stats()`, `pairwise_pop_fst()`, `windows_pairwise_pop_fst()` |
| A statistic per individual (row-wise) | the pattern in `src/gt_ind_counts.h` | `indiv_het_obs()`, `indiv_missingness()` |
| A statistic per pair of individuals | `src/pairwise_bitpacked.h` | `pairwise_ibs()`, `pairwise_king()`, `pairwise_allele_sharing()` |
| Copying genotypes into a new FBM (subset, reorder, swap alleles, stack) | `gt_copy_fbm()` (`src/gt_copy_fbm.cpp`) | `gt_update_backingfile()`, `rbind()` |
| Writing genotypes to a text format | `gt_write_text_by_locus()` / `gt_write_text_by_indiv()` (`src/gt_write_text.cpp`) | `gt_as_vcf()`, `gt_as_geno_lea()`, `gt_as_plink()` |

## Principles that apply to every kernel

- **Derive everything from `code256` at call time.** Never hard-code byte
  values (other than in the bit-packed pairwise kernel, see below). Build each
  lookup table entry from `acc.code(b)`, so imputed data, `gt_set_imputed()`
  and dosage codes keep working. A genotype is missing when
  `!(acc.code(b) > -1)`, as in the original kernels (this is true for `NA`).
- **Exact integer sums.** Pack counts as integers. When a contribution is
  fractional (e.g. pseudohaploid dosages weighted by 0.5), store twice the
  value and divide at the end. Integer sums are exact, so the results are
  bit-identical to the double sums of the old kernels.
- **A fallback for non-integer codes.** If `acc.has_int_codes(max)` is false
  (dosages), use a double path (`dosage_sums()`) that sums in the original
  order of individuals within each group, so results stay identical.
- **No R inside parallel regions.** Do all Rcpp work first: build the
  accessor, check arguments, convert group ids, allocate outputs and take raw
  pointers (`double* out = mat.begin()`). Inside OpenMP regions and inside
  the `finish()` callbacks, use only plain C++ objects. No `stop()`, no Rcpp
  vectors, no R API.
- **Preallocate the output and fill it in place.** Do not wrap kernels in
  `bigstatsr::big_apply()` with `a.combine = "rbind"` or `"c"`: combining
  blocks needs twice the memory of the result at its peak. `lut_sums()`
  processes loci in blocks and checks for user interrupts between blocks, so
  `block_size` only controls how often that happens.
- **Loci-dimension outputs keep integer names.** Return unnamed vectors or
  matrices from C++, and let R attach names with `matrix_int_names()` /
  `format_grouped_output()`. Never put loci names into `dimnames` or call
  `as.character()` on them (40M names cost ~3 GB as strings).

## `fbm_byte_acc.h` reference

All types live in `src/fbm_byte_acc.h`. Include it with
`#include "fbm_byte_acc.h"`; it includes `bigstatsr/BMCodeAcc.h`, which
brings in Rcpp.

### `FBMByteAcc acc(BM, rowInd, colInd)`

A read-only view of the raw bytes of a subset of an FBM.

- `BM` is the FBM (`Environment`); `rowInd` and `colInd` are the 1-based R
  indices (from `.gt_fbm_rows()` and `.gt_fbm_cols()`). Both are
  bounds-checked.
- Rows are stored as 0-based `size_t`. `colInd` is **not** copied: the
  accessor keeps a pointer into the R vector, so `colInd` must outlive the
  accessor (always true inside an exported function).
- Methods:
  - `nrow()`, `ncol()`: the number of selected individuals and loci;
  - `rows()`: the 0-based FBM rows of the selected individuals;
  - `col(j)`: pointer to the first byte (FBM row 0) of the j-th selected
    locus, so the genotype of FBM row `r` is `acc.col(j)[r]`;
  - `code(b)`: `code256[b]` (a plain copy, safe to use anywhere);
  - `has_int_codes(max)`: true if every non-missing code is an integer in
    `[0, max]`, i.e. integer lookup tables give exact results.

### `check_group_ids(groupIds, n, ngroups)`

Converts 0-based group ids from R (`dplyr::group_indices(.x) - 1`) to
`std::vector<size_t>`, checking the length and range.

### `GroupedRows(rows, group_ids, lut_ids, max_len)`

Orders individuals by group and then by lookup table (a stable sort, so the
original order is kept within a group), and splits them into runs
(`RowRun{start, end, group, lut, contiguous}`) of at most `max_len`
individuals that share a group and a lookup table.

- `lut_ids[i]` selects which table individual `i` uses (e.g. one table per
  ploidy class). Use `std::vector<size_t>(n, 0)` for a single table.
- Runs of consecutive FBM rows are flagged `contiguous` and read straight
  from the column pointer.
- Members: `order` (original position of the k-th sorted individual), `rows`
  (sorted 0-based rows) and `runs`.

### Lookup tables and packing

A lookup table is a `std::vector<uint64_t>` with 256 entries per table (table
`t` occupies entries `256 * t` to `256 * t + 255`). Each entry packs
`NFIELDS` counters of `W = 64 / NFIELDS` bits:

```
entry = f0 | (f1 << W) | (f2 << 2W) | ...
```

`NFIELDS` can be 1, 2, 3 or 4 (64-, 32-, 21- or 16-bit counters). Summing the
entries of a run adds every field at once.

### `max_run_length<NFIELDS>(luts)`

The longest run for which no packed counter can overflow: the field mask
divided by the largest field value of any entry. Always pass it as `max_len`
to `GroupedRows`. Each single entry must also fit within its field. With this
rule the kernels are safe for any number of individuals; long groups are
just split into several runs.

### `LutSums<NFIELDS>(acc, grouped, luts, ngroups)` (pull interface)

`sums.locus(j, tot)` fills `tot[g * NFIELDS + k]` with the 64-bit total of
field `k` for group `g` at locus `j`.

- It is `const` and thread-safe, as long as each thread has its own `tot`
  buffer of `ngroups * NFIELDS` values.
- It keeps references to `acc`, `grouped` and `luts`, which must outlive it.
- Use it when you process blocks of loci in your own loop and need a second
  phase per block (for example, all pairs of populations, or windows).

### `lut_sums<NFIELDS>(acc, grouped, luts, ngroups, ncores, block_size, finish)` (push interface)

Loops over all loci with OpenMP, in blocks of `block_size`, and calls
`finish(j, tot)` for each locus with that locus's totals. Interrupts are
checked between blocks, outside the parallel region. This is the usual
choice for a statistic per locus.

### `dosage_sums(acc, grouped, weight, copies, ngroups, ncores, block_size, finish)`

The double fallback for codes that are not small integers. For each locus it
sums `x * weight[i]` and `copies[i]` over the valid genotypes of each group,
individual by individual in the original order, and calls `finish(j, tot)`
with the two sums of group `g` in `tot[2g]` and `tot[2g + 1]` (doubles).
Build `grouped` with all `lut_ids` equal and `max_len = SIZE_MAX`.

### `alt_valid_counts(acc, group_ids, ngroups, weight, copies, ncores, block_size, finish)`

Ready-made weighted alternate allele counts and numbers of valid allele
copies by group: `finish(j, g, alt, valid)` is called for each locus and
group. It uses integer tables (one per copies/weight class, storing 2 × alt)
when codes are integers up to 255, copies are integers up to 255 and weights
are 1 or 0.5, and falls back on `dosage_sums()` otherwise. Use it for
anything built on allele frequencies.

### `GroupedSummariesDipPseudo`

Per-group alternate allele sums, valid allele counts and 2 × heterozygote
counts for diploid and pseudohaploid data, via `locus(j, alt_sum, n_alleles,
het2, buffer)`. It uses `LutSums<3>` with one table per ploidy, and a double
fallback otherwise. Each thread needs its own buffer of `buffer_size()`
`uint64_t` values. It is used by the population statistics and pairwise Fst
kernels, which process loci in blocks of 1,024 and then compute statistics
per pair or per window.

### `genotype_count_lut(acc)` (in `src/hwe.cpp`)

Packs the counts of dosage 0, 1, 2 and missing genotypes into four 16-bit
fields (`LutSums<4>`). Used by `loci_hwe()` and `qc_report_loci()`; declare
it in your file to reuse it, as `src/loci_qc_cpp.cpp` does.

## Worked example: heterozygotes by locus and group

This kernel counts heterozygous and valid genotypes at each locus for each
group of diploid individuals. It was compiled and checked against
`show_genotypes()` on the example data (identical results, including
reordered individuals and loci). In the package, put it in its own file in
`src/` and run `Rcpp::compileAttributes()`.

```cpp
#include "fbm_byte_acc.h"

// Count heterozygous and valid (non-missing) genotypes at each locus, for
// each group of individuals (diploid data). Returns a list of two loci x
// groups matrices.
// [[Rcpp::export]]
List grouped_het_count_cpp(Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd,
                           const IntegerVector& groupIds,
                           int ngroups,
                           int ncores,
                           int block_size = 100000) {
  // 1. all Rcpp work happens before any parallel region
  FBMByteAcc acc(BM, rowInd, colInd);
  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);

  // 2. one lookup table built from code256 at call time: field 0 counts
  // heterozygotes, field 1 counts valid genotypes (32-bit fields)
  std::vector<uint64_t> lut(256, 0);
  for (int b = 0; b < 256; b++) {
    double x = acc.code(b);
    if (x > -1) { // not missing
      uint64_t het = (x == 1) ? 1 : 0;
      lut[b] = het | (uint64_t(1) << 32);
    }
  }

  // 3. runs of individuals by group, short enough to never overflow
  GroupedRows grouped(acc.rows(), group_ids, std::vector<size_t>(n, 0),
                      max_run_length<2>(lut));

  // 4. preallocate the outputs, and write through raw pointers
  NumericMatrix n_het(m, ngroups), n_valid(m, ngroups);
  double* p_het = n_het.begin();
  double* p_valid = n_valid.begin();

  lut_sums<2>(acc, grouped, lut, ngroups, ncores, block_size,
              [&](size_t j, const uint64_t* tot) {
                for (int g = 0; g < ngroups; g++) {
                  p_het[j + g * m] = static_cast<double>(tot[2 * g]);
                  p_valid[j + g * m] = static_cast<double>(tot[2 * g + 1]);
                }
              });
  return List::create(_["n_het"] = n_het, _["n_valid"] = n_valid);
}
```

Points to note:

- For an ungrouped statistic, pass `std::vector<size_t>(n, 0)` as the group
  ids and `ngroups = 1` (see `missingness_cpp()` in
  `src/grouped_missingness_cpp.cpp`).
- For polyploid or pseudohaploid data, build one table per ploidy class
  (`lut_ids`), as `gt_pi_polyploid()` and `loci_qc_cpp()` do.
- If the statistic needs codes other than 0, 1 and 2 (e.g. dosages), check
  `acc.has_int_codes()` and provide a double fallback.

## Calling a kernel from R

Follow the existing `loci_*` wrappers (e.g. `R/loci_missingness.R`):

```r
geno_fbm <- .gt_get_fbm(.x)
res <- grouped_het_count_cpp(
  BM = geno_fbm,
  rowInd = .gt_fbm_rows(.x),            # FBM rows of the individuals
  colInd = show_loci(.x)$big_index,     # FBM columns of the loci
  groupIds = dplyr::group_indices(.x) - 1,
  ngroups = dplyr::n_groups(.x),
  ncores = n_cores,
  block_size = block_size
)
# attach loci and group names without creating character row names
res <- format_grouped_output(
  out_mat = res$n_het,
  group_ids = dplyr::group_keys(.x) %>% dplyr::pull(1),
  loci_names = loci_names(.x),
  type = type
)
```

- Pass `n_cores` (default `bigstatsr::nb_cores()`) straight to the kernel;
  do not also parallelise in R.
- Keep `block_size` as an argument for consistency with the other `loci_*`
  functions; it now only sets how often interrupts are checked.
- `show_loci(.x)$big_index` and `loci_names(.x)` do not copy the loci table.
- Ploidy per individual is `indiv_ploidy(.x)`; pseudohaploid individuals have
  ploidy 1 but store their allele as a 0/2 dosage.

## Per-individual statistics: `gt_ind_counts.h`

Row-wise statistics (one value per individual, summed over loci) cannot use
`lut_sums()`, which sums over individuals. `src/gt_ind_counts.h` (used by
`indiv_het_obs()` and `indiv_missingness()`) shows the pattern:

- Map each byte through a 256-entry table (built from `code256`, one table
  per ploidy class) to a packed 64-bit increment (heterozygotes in the low 32
  bits, missing in the high 32 bits).
- Process a few loci at a time (`N_COLS`), so that each individual's counter
  is read and written once per group of loci rather than once per genotype.
- Template flags remove unneeded loads in the common case (all FBM rows in
  their original order, a single table).

It predates `fbm_byte_acc.h` and keeps its own small helpers. It could be
moved onto `FBMByteAcc` for consistency, with about the same speed. For a
new row-wise statistic, copy this structure and change the table.

## Pairwise statistics: `pairwise_bitpacked.h`

For statistics summed over loci for every pair of diploid individuals,
`src/pairwise_bitpacked.h` packs genotypes into bit planes (64 loci per
64-bit word) and counts with AND/XOR and popcount, as PLINK does. A new
pairwise statistic is a small struct (see `PwIbs`, `PwKing` and
`PwAlleleSharing`) that says which planes it needs and how to combine them
for one word of a pair.

- The pair loop is compiled twice on x86 (with and without the `popcnt`
  target) and chosen at run time with `__builtin_cpu_supports()`. On other
  architectures the generic loop is used, and the builtin is already a
  native instruction.
- As in the original kernels, this header reads raw bytes directly: 0, 1 and
  2 are valid genotypes and any other byte (3 for missing, 4 to 6 for imputed
  calls) is treated as missing, whatever `code256` says. Keep this in mind
  if a new statistic should use imputed data.

## Copying and exporting genotypes

- `gt_copy_fbm(BM_from, BM_to, rowInd, colInd, col_offset, row_offset, swap,
  swap_lut)` copies the raw bytes of selected rows and columns of one FBM
  into another, optionally at a row offset (to stack two datasets) and
  swapping alleles through a 256-entry table at flagged loci. Raw bytes are
  copied as they are, so imputed values are preserved. Use it instead of
  decoding genotypes into R and writing them back.
- `gt_write_text_by_locus()` (one line per locus: VCF, geno) and
  `gt_write_text_by_indiv()` (one line per individual: ped, raw) write text
  from the raw bytes through a table of 256 strings that R builds from
  `code256`. To add a new text format, build the table in R with the same
  recoding as the format needs, and reuse these functions.

## Checklist for a new kernel

1. Write a reference implementation in R (`show_genotypes()` on small data)
   for the tests.
2. Build the accessor and tables from `code256`, with a fallback for
   non-integer codes if the statistic accepts dosages.
3. Use `max_run_length()` for `GroupedRows`; never hard-code a run length.
4. Keep all Rcpp and R API calls outside the parallel regions and the
   `finish()` callbacks.
5. Preallocate outputs; do not use `big_apply()` with `rbind`/`c` to combine
   blocks.
6. Test against the reference, in the existing test file for the function,
   with small inline data: missing genotypes; reordered and filtered
   individuals (`dplyr::slice()`) and loci (`select_loci()`); grouped data;
   a single individual and a single locus; pseudohaploid and polyploid data
   where supported; imputed data (`gt_impute_simple()` +
   `gt_set_imputed()`); and the dosage fallback (set
   `attr(x$genotypes, "fbm")$code256 <- bigsnpr::CODE_DOSAGE`, then restore
   it). Keep tests small and fast; no timing tests.
7. Benchmark at -O2, never from a plain `load_all()`, which compiles the C++
   at -O0 (2–3× slower). Run `pkgbuild::clean_dll(path)`, then
   `pkgbuild::compile_dll(path, debug = FALSE)`, then
   `pkgload::load_all(path, compile = FALSE)`, and after any C++ edit run
   `compile_dll(path, debug = FALSE)` again before loading. Compare against
   the previous commit installed into its own library
   (`git archive <ref> | tar -x -C <dir>`, then
   `R CMD INSTALL --library=<lib> <dir>`), loaded in a separate R process.
   Also compare against the memory floor: the time to read every byte once
   (about 0.4 s for HGDP on a laptop). Ungrouped kernels should be close to
   it.

## Pitfalls

- **`colInd` lifetime.** `FBMByteAcc` keeps a pointer into `colInd`; do not
  pass a temporary vector that is destroyed before the accessor is used.
- **Overflow.** A table entry with a field larger than 1 (e.g. 2 × dosage)
  shortens the safe run length; `max_run_length()` accounts for it, but only
  if every table is included in `luts`.
- **Thread-local buffers.** With `LutSums::locus()` in your own OpenMP loop,
  allocate one `tot` buffer per thread (inside the parallel region).
- **Summation order.** If outputs must stay bit-identical to an existing
  function, keep the same order of floating-point operations (integer sums
  are always exact). The pop statistics kernels reproduce `rowMeans()` and
  `colMeans()` in long double for this reason.
- **Interrupts.** Never call `checkUserInterrupt()` inside a parallel
  region; `lut_sums()` does it between blocks.
- **Memory bandwidth.** On a laptop, ungrouped kernels saturate memory
  bandwidth with one or two threads. More threads help mainly for grouped
  kernels with many groups and for the per-pair or per-window phases.
