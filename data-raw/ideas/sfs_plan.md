# Plan: site frequency spectra (`pop_sfs()`, `pairwise_pop_sfs()`, `nwise_pop_sfs()`)

A family of functions that compute one-dimensional, two-dimensional
(joint) and multi-dimensional site frequency spectra (SFS) from a grouped
`gen_tibble`. Missing data are handled the way easySFS does it: each
population is projected down to a fixed number of allele copies with the
hypergeometric projection of dadi. The aim is that, for the same data,
projections and options, the spectra match easySFS's output.

## Status (branch `sfs`)

**Decisions taken:**

- **Choosing the sample size.** `n_proj = NULL` (the default) projects each
  population to the number of allele copies that maximises its number of
  segregating sites, which is the rationale of easySFS and the dadi manual.
  Ties go to the largest projection. `pop_sfs_preview()` shows the full curve
  and flags this choice in a `best` column, and a message lists the chosen
  values.
- **Return type.** A tidy tibble by default, like the other functions;
  `type = "array"` returns `sfs` arrays.
- **Cell order.** The order of the populations is respected by default:
  rows are sorted by the first population, then the second, and so on, so
  the last population changes fastest, as in dadi and easySFS.
  `fastest = "first"` makes the first population change fastest. This is set
  in `pairwise_pop_sfs()`, `nwise_pop_sfs()` and `tidy.sfs()`.

**Implemented:**

- stage 1 and stage 2 of the engine (`src/sfs_cpp.cpp`, `R/sfs_helpers.R`);
- `pop_sfs()`, `pairwise_pop_sfs()`, `nwise_pop_sfs()` and
  `pop_sfs_preview()`;
- the `sfs` class with `print()` and `tidy()`;
- the reference tests (`helper_sfs_reference.R`, a port of easySFS and
  dadi);
- `data-raw/create_easysfs_test_data.R` and `test_easysfs.R`, which is
  skipped until the easySFS fixtures are generated.
- `write_sfs()`, which writes one spectrum in dadi or fastsimcoal2 format,
  with `fastest` and with easySFS's `--total-length` (`total_length`) and
  `--dtype int` (`round`);
- `ancestral` (in `pop_sfs()`, `pairwise_pop_sfs()`, `nwise_pop_sfs()`,
  `pop_sfs_preview()` and `write_sfs()` for a gen_tibble): the ancestral
  allele of each locus, as a loci table column or a vector, compared ignoring
  case. Loci
  whose ancestral allele is the alternate allele are flipped in the counting
  pass, and loci without a usable ancestral allele are dropped (counted in
  `n_loci_unpolarised`), as dadi does with an outgroup allele. It is ignored,
  with a warning, for folded spectra, which polarisation does not change;
- `autoplot()` methods for the tidy outputs of `pop_sfs()` (class
  `pop_sfs_tbl`, bar plots with one panel per population) and
  `pairwise_pop_sfs()` (class `pairwise_sfs_tbl`, heatmaps with one panel per
  pair), for `sfs` objects of one or two populations, and for the output of
  `pop_sfs_preview()` (class `sfs_preview_tbl`, with the default projection
  marked). Masked cells are hidden by default. Heatmaps use a single blue
  ramp on a log scale that starts at `min_sites` (a thousandth of the
  largest cell by default, like dadi's `vmin`), because the projection spreads
  each locus over many cells with tiny fractions of a site;
- the block bootstrap: `bootstrap_pop_sfs()`, `bootstrap_pairwise_pop_sfs()`,
  `bootstrap_nwise_pop_sfs()`, and the `n_bootstrap` argument of
  `write_sfs()` for a gen_tibble. They work as follows:
  - the blocks are those of the jackknife (`.jackknife_block_ids()`);
  - the counting pass keys each configuration by block as well, so the
    genotypes are read once;
  - the projection sorts the rows by block, completes the spectrum of each
    block in turn, and adds it to every replicate, weighted by the number of
    times the block was drawn (`rmultinom()`, so `set.seed()` makes it
    reproducible);
  - replicates use the projections chosen with all loci;
  - `write_sfs()` projects the replicates in chunks of at most 1e7 cells, to
    limit memory;
- `write_sfs()` for a grouped gen_tibble (first written as
  `gt_as_easysfs()`, then split into `sfs_as_dadi()` and `sfs_as_simcoal()`,
  and finally made methods of `write_sfs()`: `write_sfs.sfs()` for a
  spectrum and `write_sfs.gen_tbl()` for a gen_tibble). It computes, in one
  pass, all the spectra that easySFS writes (1D, each pair, all populations)
  and writes them in `format` ("dadi" or "fastsimcoal2"), with easySFS's
  file names,
  directly in `out_dir`. Writing them to `dir/dadi` and `dir/fastsimcoal2`
  gives easySFS's layout. This was checked
  against emulated easySFS outputs, written by a literal R port of easySFS's
  writing code; `test_easysfs.R` compares every file once the real fixtures
  exist.

**Benchmark.** The data were 1000 diploids x 100k loci, 10 groups of 100,
with 5% missing genotypes. The branch was installed at -O2 in its own
library, as `CLAUDE.md` requires.

| Function | 1 core | 4 cores |
|---|---|---|
| grouped `loci_alt_freq()` (one pass, for reference) | 0.14 s | 0.06 s |
| `pop_sfs_preview()` | 0.13 s | 0.07 s |
| `pop_sfs(n_proj = 160)` | 0.35 s | 0.27 s |
| `pairwise_pop_sfs()`, 45 pairs, `n_proj = 40` | 2.1 s | 1.4 s |
| `nwise_pop_sfs()`, 3 populations, `n_proj = 40` | 0.35 s | 0.33 s |

- Bootstrap timings: 100 replicates, 1000 blocks of 100 SNPs, -O2 build.
  Without bootstrap, the timings are unchanged.

  | Function | 1 core | 4 cores |
  |---|---|---|
  | `bootstrap_pop_sfs(n_proj = 160)` | 0.85 s | 0.68 s |
  | `bootstrap_pairwise_pop_sfs()`, 45 pairs, `n_proj = 40` | 11.2 s | 4.8 s |
  | `bootstrap_nwise_pop_sfs()`, 3 populations, `n_proj = 40` | 6.1 s | 4.1 s |

  The projection dominates. Within a block, almost every locus has its own
  configuration (about 99k rows per spectrum), so the factorisation saves
  little. Each block spectrum is then added to every replicate. The block
  spectra are buffered (32 at a time) and added by tiles of cells, which cut
  the `nwise` time from about 9 s.
- With missing data, each pair has about 60k distinct configurations.
  Counting them takes about 1 s, and projecting them about 1 s.
- In the projection, building the probability tables and copying the
  configurations into C++ is serial, so it is the next thing to speed up if
  needed.

**Still to do:**

- generate the easySFS fixtures (run `data-raw/create_easysfs_test_data.R`);
- verify the cell order of fsc MSFS files against the fsc manual. In fsc's
  2D files the first population changes fastest along each row, which
  suggests that fsc may expect `fastest = "first"` for MSFS files too.
  easySFS writes them in dadi's order (`fastest = "last"`, the default);
- check whether fsc expects `_DSFS.obs` rather than `_MSFS.obs` for an
  unfolded multidimensional spectrum (easySFS always writes `_MSFS.obs`);
- an HGDP-sized benchmark.

Contents:

1. [What easySFS does](#1-what-easysfs-does-the-logic-to-reproduce)
2. [User-facing functions](#2-user-facing-functions)
3. [The `sfs` object](#3-the-sfs-object)
4. [Computation engine](#4-computation-engine)
5. [Ploidy, imputation, polarisation and edge cases](#5-ploidy-imputation-polarisation-and-edge-cases)
6. [Output files (dadi, fastsimcoal2)](#6-output-files-dadi-fastsimcoal2)
7. [Tests](#7-tests)
8. [Files to add or change](#8-files-to-add-or-change)
9. [Phases](#9-phases)
10. [Decisions for the maintainers](#10-decisions-for-the-maintainers)

---

## 1. What easySFS does (the logic to reproduce)

This section is based on the current `easySFS.py` and its bundled
`Spectrum.py`, a copy of dadi's `Spectrum` class
(github.com/isaacovercast/easySFS, master branch).

**Counting alleles (`make_datadict`).** For each site and population,
easySFS counts the reference and alternate allele copies among the
individuals that have a genotype call:

- diploids: `0/0` adds 2 REF, `1/1` adds 2 ALT, and `0/1` adds 1 of each;
- haploids (`--ploidy 1`): `0` and `1` add one copy each;
- for ploidy above 2, each allele in the genotype string is counted.

Missing genotypes (`./.`) add nothing. A site therefore has its own number of
called copies, `n`, in each population, and `k` derived (ALT) copies.

**Polarisation.** The REF allele is always taken as ancestral
(`outgroup_allele = REF`). With `--unfolded`, the derived count is the ALT
count. Without it (the default), the spectrum is folded.

**Projection (`Spectrum._from_count_dict` and `_cached_projection`).**

- Sites are first grouped into configurations: a tuple of `(n, k)` for each
  population, with a count of sites per tuple (`count_data_dict`).
- To project population `i` to `p_i` copies, a site with `n_i >= p_i`
  contributes the hypergeometric probabilities
  `P(j) = C(p,j) C(n-p,k-j) / C(n,k)` for `j = 0..p`. These are equal to
  `dhyper(j, k, n - k, p)`; I checked this numerically against dadi's log-gamma
  formula for all `n <= 60`, with differences below 1e-13.
- A site with `n_i < p_i` in any population of the spectrum gets a vector of
  zeros, so it is dropped. **This is how missing data are handled.**
- The contribution of a site to a joint SFS is the outer product of the
  projection vectors of its populations, multiplied by the site count. The
  spectrum is therefore made of fractional counts.
- Monomorphic sites, and sites that project to monomorphic, are added to the
  corner cells. Those cells stay in the data but are masked.

**Folding (`Spectrum.fold`).** A joint spectrum is folded on the *total*
derived count across all its populations, not one axis at a time:

- cells whose total is above `floor(sum(p) / 2)` are added to the reversed
  array and then set to 0;
- cells whose total is exactly `sum(p) / 2` are averaged with their mirror
  cell (`-0.5 * a + 0.5 * reverse(a)`);
- the mask covers the original mask, its reverse, and the folded-out cells.

**Outputs.** With `--proj p1,p2,...` (one value per population, in population
order or in the order given by `--order`), easySFS writes:

- a 1D SFS for every population;
- a 2D SFS for every pair, in `itertools.combinations` order (the same order
  as `utils::combn()`);
- one multi-dimensional SFS with all populations.

Each is written in dadi format and in fastsimcoal2 format; see
[section 6](#6-output-files-dadi-fastsimcoal2).

Each spectrum only uses the sites that pass the projection for *its own*
populations. The marginals of a 2D SFS are therefore not, in general, equal to
the 1D spectra.

**Preview (`--preview`).** For each population on its own, easySFS
computes `S()` of the projected 1D spectrum for every `p` from 2 to
`ploidy * n_individuals`. `S()` is the sum of the cells that are not corners,
that is the number of segregating sites. It prints `S()` rounded, and the
user picks the `p` that keeps the most segregating sites without losing too
much sample size. In 1D, `S()` is the same whether the spectrum is folded or
not, because REF is always a valid outgroup.

**`--total-length L`.** Before writing, cell 0 is set to `L` minus the
sum of all the other cells. It is an error if `L` is not greater than the
current sum of the spectrum.

**`--dtype int`.** The cells are rounded with `np.rint` (round half to even,
like R's `round()`).

**Out of scope for the SFS functions (done by selecting loci instead).**

- By default, easySFS samples one random SNP per `CHROM` (one per RAD locus).
  `-a` keeps all SNPs, and `--window-bp` / `--window-snp` thin by window.
- `--GQ` is parsed but never used.
- `--resample-vcf` (hidden) resamples SNPs with replacement.

In `tidypopgen`, thinning is done with `select_loci()`, `loci_ld_clump()`, or
`slice_sample()` per chromosome on `show_loci()`. The documentation should
show this, as RAD users will expect it.

## 2. User-facing functions

The functions follow the `pop_` / `pairwise_` / `nwise_` naming, so pkgdown
already places them in the "Population" and "Pairwise" sections. They take a
grouped `gen_tibble` and return `sfs` objects
([section 3](#3-the-sfs-object)).

As implemented on the `sfs` branch (`ancestral` and `total_length` are not
implemented yet):

```r
pop_sfs(.x, n_proj = NULL, folded = TRUE, type = c("tidy", "array"),
        n_cores = bigstatsr::nb_cores(),
        block_size = bigstatsr::block_size(nrow(.x), 1), ...)

pairwise_pop_sfs(.x, n_proj = NULL, folded = TRUE,
                 type = c("tidy", "array"), fastest = c("last", "first"),
                 n_cores = ..., block_size = ...)

nwise_pop_sfs(.x, n_proj = NULL, folded = TRUE, type = c("tidy", "array"),
              fastest = c("last", "first"), max_cells = 1e7,
              n_cores = ..., block_size = ...)

pop_sfs_preview(.x, n_cores = ..., block_size = ...)
```

The tidy outputs have these columns:

- `pop_sfs()`: `group` (or one column per grouping variable),
  `allele_count`, `n_sites` and `masked`;
- `pairwise_pop_sfs()`: `<var>_1`, `<var>_2`, `allele_count_1`,
  `allele_count_2`, `n_sites` and `masked`;
- `nwise_pop_sfs()`: one column per population (named by its label),
  `n_sites` and `masked`.

The tidy outputs carry the attributes `n_proj` and `folded`.

The plan as first written is below.

What each function returns:

- **`pop_sfs()`**:
  - grouped input: a named list of 1D `sfs`, one per group;
  - ungrouped input: a single `sfs` for all individuals.
  
  It has `tbl_df`, `grouped_df` and `vctrs_bigSNP` methods, like
  `pop_tajimas_d()`.
- **`pairwise_pop_sfs()`**: a list of 2D `sfs` objects, one per pair of
  groups, in `utils::combn()` order and named with `.pair_labels()` (for
  example `pop1.pop2`). The first population of the pair is `dim 1`. Stops
  if there are fewer than 2 groups.
- **`nwise_pop_sfs()`**: a single multi-dimensional `sfs` of all the groups,
  in group order (easySFS's MSFS). Users choose the populations with
  `filter()`, and their order with the factor levels of the grouping
  variable, which is also how easySFS's `--order` maps to `tidypopgen`. It
  stops before allocating when `prod(n_proj + 1) > max_cells`, as the array
  grows exponentially with the number of populations.
  - `nwise_pop_pbs()` uses all 3-way combinations. All `n`-way combinations
    of spectra are not needed for dadi or fsc, and would explode in size, so
    they are left out.
- **`pop_sfs_preview()`**: a tibble with one row per group and projection:
  - the group columns, from `.group_columns()`;
  - `n_proj`, from 2 to the group's maximum number of allele copies;
  - `n_seg_sites`, `S()` as a double, not rounded;
  - `n_loci`, the number of loci with at least `n_proj` called copies.
  
  It has an `autoplot()` method that draws `n_seg_sites` against `n_proj`,
  with one panel per group.

**Shared arguments:**

- **`n_proj`** (`NULL`, a single integer, or one per group):
  - a single value is recycled;
  - a named vector is matched to the group labels; an unnamed vector is taken
    in group order;
  - `NULL` means each group's maximum number of allele copies (the sum of the
    ploidy of its individuals): no projection, so loci with missing data in a
    group are dropped. A message gives the number of loci dropped and points
    to `pop_sfs_preview()`;
  - each value must be between 1 and the group's maximum;
  - the same `n_proj` per group is used for every spectrum, as in easySFS.
- **`folded`**: `TRUE` by default, as in easySFS.
- **`ancestral`**: `NULL` means that `allele_ref` is ancestral (easySFS).
  Otherwise, give the name of a loci table column, or a character vector with
  one allele per locus ([section 5](#5-ploidy-imputation-polarisation-and-edge-cases)).
- **`total_length`**: `--total-length`, applied after folding to each
  spectrum.
- **`n_cores`, `block_size`**: as in `pop_tajimas_d()`.

## 3. The `sfs` object

An `sfs` is an S3 class over a numeric array. It has one dimension per
population, of extent `n_proj + 1`.

**Attributes:**

| Attribute | Content |
|---|---|
| `dimnames` | `d<i>_<j>`-style labels |
| `pop_ids` | the group labels |
| `n_proj` | the projection of each population |
| `folded` | whether the spectrum is folded |
| `mask` | a logical array with the corners and folded-out cells |
| `n_loci` | the loci that contributed, i.e. those with `n >= n_proj` in all its populations |
| `n_loci_dropped` | the loci without enough called copies |

- The sum of the array equals `n_loci` (projection preserves mass) until
  `total_length` is applied. This is a useful invariant for the tests.
- Keeping the masked cells in the array, and the mask beside it, gives
  exactly what dadi writes, with no extra information.
- A 1D `sfs` is a 1D array, not a bare vector, so all the code paths are the
  same.

**Methods:**

- `print()`;
- `tidy()`, a long tibble with one column of derived counts per population,
  `count` and `masked`;
- `autoplot()`, a bar plot in 1D and a `geom_tile` heatmap in 2D; it stops
  for higher dimensions;
- `write_sfs()`, see [section 6](#6-output-files-dadi-fastsimcoal2).

## 4. Computation engine

All three functions call one internal engine. The engine computes any set of
spectra in a single pass over the genotypes.

```r
.sfs_compute(.x, spectra, n_proj, folded, ancestral, total_length,
             n_cores, block_size)
# spectra: a list of integer vectors of group indices, e.g.
#   pop:      as.list(seq_len(G))
#   pairwise: split(t(utils::combn(G, 2)), ...)
#   nwise:    list(seq_len(G))
```

### Stage 1 (C++, one pass): configuration counts

This is dadi's `count_data_dict`, run on the FBM.

The new file `src/sfs_config_counts.cpp` reuses `alt_valid_counts()` from
`fbm_byte_acc.h`, with the same `weight` and `copies` as
`tajimas_d_window_sums_cpp()`:

- pseudohaploids (ploidy 1 with `is_pseudohaploid`) have weight 0.5;
- every individual has `copies = ploidy`.

This gives, for every locus and group, `n` (valid copies) and `k` (derived
copies), as exact integer sums from the lookup tables.

**Per block.** For each block of loci, the `finish()` callback stores `n`
and `k` per group in block buffers, as `b_pi` and `b_valid` do in the
Tajima kernel. If the locus is flipped by `ancestral`, it stores
`k = n - alt`. For loci that `ancestral` cannot polarise, it stores a "drop"
flag, used only when `folded = FALSE`.

**`block_done(j0, j1)`.** This adds the loci of the block to each requested
spectrum:

- skip the locus if any population of the spectrum has `n < n_proj`, and
  increment `n_loci_dropped`;
- otherwise, encode each population's `(n, k)` as a triangular index
  `c = n (n + 1) / 2 + k`, combine the indices of the spectrum into one key,
  and increment its count.

The loop over spectra inside `block_done` runs in parallel
(`#pragma omp parallel for` over spectra). Each spectrum's table belongs to
one thread, so there are no locks and no merge step. Counts are integers, so
the results do not depend on `n_cores` or `block_size`, which the package
requires.

**Containers:**

- 1D: a dense triangular array of `int64` counts per group (`N (N + 1) / 2`
  cells; 16 MB for `N = 2000` copies).
- 2D / nD: `std::unordered_map<uint64_t, int64_t>`, with a mixed-radix key
  built from the per-population config indices. If the key cannot fit in 64
  bits, fall back to a key of type `std::vector<uint32_t>`; in practice the
  64-bit key always fits.

**Output.** For each spectrum, the distinct configurations (an integer matrix
of `n` and `k` per population), their counts, and `n_loci_dropped`, sorted by
key so the next stage is deterministic.

Other options and why they were rejected:

- **Return the `n` and `k` matrices to R** (the same as
  `grouped_alt_freq_dip_pseudo_cpp(as_counts = TRUE)`). This allocates
  `m x 2G` doubles, 1.6 GB for 10M loci and 10 groups, which goes against
  the kernel guide (`data-raw/speed_accessors_guide.md`).
- **Project each locus inside the kernel.** This costs
  `m x prod(n_proj + 1)` instead of `distinct configs x prod(n_proj + 1)`.
  Without missing data there are at most `N + 1` configurations per
  population, and repeats are very common even with missing data, so the
  configuration counts are much faster for 2D and nD.

### Stage 2 (C++): projection and folding

This is `src/sfs_project.cpp`, the equivalent of dadi's `_from_count_dict` and
`fold`.

**Projection.**

- For each population, precompute the hypergeometric vectors for every `n` in
  `[n_proj, N]` and `k` in `[0, n]`. Use `R::dhyper()`, or log-space
  `lchoose()` like dadi; they agree to about 1e-13.
- For each configuration, in key order, add
  `count x (outer product of its vectors)` to the dense output array. The
  array is stored in R's column-major order, so `dim 1` changes fastest.
- The cost is `distinct configs x prod(n_proj + 1)`, the same as dadi.

**Folding.** If `folded`, fold on the total derived count across populations,
using dadi's algorithm exactly, including the averaging of ambiguous cells
and the mask. Do this in C++ or in R; R is fine, as the array is small
compared with the genotypes.

**Total length.** Apply `total_length` last, as easySFS does.

**Preview.** `pop_sfs_preview()` only runs stage 1, with 1D spectra and the
smallest projection (2) for each group, so no locus with at least 2 copies is
dropped. It returns the dense `(n, k)` histogram `H_g` per group. It then
computes, in R, for every `p`:

```
S_g(p) = sum over n >= p and k of H_g[n, k] * (1 - dhyper(0, k, n - k, p)
                                                - dhyper(p, k, n - k, p))
```

There is one pass over the genotypes for all the projections, where easySFS
recomputes the full spectrum for each `p`.

### Speed and memory

Stage 1 is a grouped `alt_valid_counts()` pass plus small per-block work, so
it should run close to `loci_alt_freq()` on grouped data. The extra cost for
pairwise is G(G-1)/2 hash updates per locus, spread across threads by
spectrum. Benchmark it on HGDP-sized data, following `CLAUDE.md`:

- never time a plain `load_all()`;
- use `clean_dll()`, then `compile_dll(debug = FALSE)`, then
  `load_all(compile = FALSE)`;
- check with `objdump` that no -O0 code is left.

## 5. Ploidy, imputation, polarisation and edge cases

**Ploidy.** `alt_valid_counts()` handles diploid, polyploid, mixed-ploidy
and pseudohaploid data. Each group's maximum projection is the sum of the
ploidy of its individuals, so a pseudohaploid counts as one copy. This is the
same as easySFS for uniform ploidy (`--ploidy 1` for haploids and
pseudohaploids). Mixed ploidy within a population is supported here but not by
easySFS.

**Codes.** The SFS needs integer allele counts:

- if `acc.has_int_codes()` is false (dosage codes, e.g.
  `bigsnpr::CODE_DOSAGE`), stop with a clear error;
- with imputed integer codes (`gt_set_imputed(TRUE)`), the imputed genotypes
  are used like everywhere else in the package. The docs should say that the
  projection makes imputation unnecessary, and that `gt_set_imputed(FALSE)`
  is usually the right choice before computing an SFS.

**Polarisation.**

- With `folded = FALSE` and `ancestral = NULL`, `allele_ref` is ancestral,
  as in easySFS.
- The docs must warn that operations that swap alleles change what REF means
  (for example `rbind()` matching against a reference, or flipping
  alleles).
- With an `ancestral` column or vector, loci where it equals `allele_alt` are
  flipped (`k = n - alt`). Loci where it is `NA` or matches neither allele are
  dropped from unfolded spectra, as dadi does with an AA outgroup. They are
  counted in `n_loci_dropped`. This is a per-locus integer vector passed to
  stage 1, so it costs nothing when unused.

**Groups.**

- Empty groups (only from `.drop = FALSE`) are dropped with a message, as
  `gt_extract_f2()` does: an SFS is meant for external inference programs.
- With groups from several grouping variables, labels are pasted with `.`
  (`.group_info()`).

**Loci.** Monomorphic loci count towards the corners, which are masked. A
gen_tibble with no loci, or with no locus passing the projection, gives an
all-zero spectrum with `n_loci = 0`, and a warning.

## 6. Output files (dadi, fastsimcoal2)

`write_sfs(x, file, format = c("dadi", "fsc"), round = FALSE)` writes one
spectrum. The formats are those written by easySFS:

**dadi:**

- line 1: the dimensions (`n_proj + 1`), then `folded` or `unfolded`, then
  the population labels in quotes;
- line 2: the cells, in C order (last population fastest), with `%.16g`;
- line 3: the mask, as 0/1.

Keep in mind that R arrays are column-major, so `aperm()` or `rev(dim)` is
needed when flattening.

**fsc, 1D:**

- `1 observation`;
- the headers `d0_0` to `d0_p`, separated by tabs;
- the dadi data line, followed by an empty line (easySFS writes an extra
  `"\n"`);
- file name: `<pop>_MAFpop0.obs`, or `_DAFpop0.obs` when unfolded.

**fsc, 2D:**

- file name: `<prefix>_jointMAFpop<j>_<i>.obs`, where `i < j` are the 0-based
  population indices;
- the rows are population `j` (`d<j>_*`) and the columns population `i`
  (`d<i>_*`);
- the header line starts with a tab, and each row is the row label, a tab,
  then values separated by spaces.

**fsc, MSFS:**

- file name: `<prefix>_MSFS.obs`;
- `1 observations. No. of demes and sample sizes are on next line.`;
- the number of populations, a tab, then the projections separated by
  spaces;
- the dadi data line.

**Ordering to verify.** easySFS writes the MSFS cells in dadi's order, with
the last population fastest. Check this against the fastsimcoal2 manual
before release. If they disagree, follow fsc and document the difference from
easySFS.

**`round = TRUE`** gives `--dtype int` (R's `round()`, half to even, the same
as `np.rint`).

`gt_as_easysfs(.x, n_proj, out_dir, prefix, folded = TRUE, total_length = NULL,
round = FALSE)` is a drop-in replacement for an easySFS run. It calls the
three engines in one pass and writes the same directory layout and file
names:

- `dadi/<pop>-<p>.sfs`;
- `dadi/<pop1>-<pop2>.sfs`;
- `dadi/<pop1>-...-<popK>.sfs`;
- `fastsimcoal2/...`.

Its name follows the `gt_as_*` export family (pkgdown "Exporting" section).

## 7. Tests

1. **R reference implementation** (`tests/testthat/helper_sfs_reference.R`).
   This is a literal port of `make_datadict`, `count_data_dict`,
   `_from_count_dict` and `fold`, working from `show_genotypes()`. Compare
   `pop_sfs()`, `pairwise_pop_sfs()` and `nwise_pop_sfs()` with it:
   - on small inline data with ragged missingness;
   - with projections from 1 to the maximum, folded and unfolded;
   - for 1, 2 and 3 populations;
   - with reordered or filtered individuals and loci;
   - with two grouping variables.
2. **easySFS fixtures.** `data-raw/create_easysfs_test_data.R`, modelled on
   `data-raw/create_pixy_test_data.R`:
   - create a conda environment (python, numpy, scipy, pandas) and clone
     easySFS at a **pinned commit**;
   - simulate a biallelic VCF: 3 populations of unequal size, several CHROMs,
     ragged missingness, some sites monomorphic in a population;
   - run `easySFS.py -a -y -f`, with `-a` so that it does not sample one SNP
     per CHROM at random, for:
     - folded;
     - `--unfolded`;
     - `--total-length`;
     - `--dtype int`;
     - `--order`;
     - `--preview` (save the output);
     - a haploid VCF with `--ploidy 1`;
   - save the outputs in `tests/testthat/testdata/easysfs`.
   
   `test_easysfs.R` then checks:
   - the cells, to a tolerance;
   - the masks;
   - the preview `S`, after rounding;
   - with `write_sfs()` / `gt_as_easysfs()`, that the files are identical,
     byte for byte where the formatting allows.
3. **Invariants:**
   - with no missing data and `n_proj = NULL`, the unfolded 1D SFS equals
     `tabulate()` of the ALT counts;
   - the sum of the spectrum equals `n_loci`;
   - with no missing data, the marginals of the unfolded 2D spectra equal the
     unfolded 1D spectra (not when folded, because joint folding uses the
     total count);
   - folding keeps the total;
   - the results are the same for `n_cores = 1` and `2` and for any
     `block_size`;
   - with no missing data, Tajima's D computed from the unprojected 1D SFS
     equals `pop_tajimas_d()`.
4. **Data types and edge cases** (the checklist in
   `speed_accessors_guide.md`):
   - pseudohaploid, polyploid and mixed-ploidy data;
   - imputed data, which uses the imputed codes;
   - dosage codes, which give an error;
   - a single individual or locus;
   - empty groups;
   - `pairwise_pop_sfs()` with fewer than 2 groups (error);
   - the `max_cells` guard of `nwise_pop_sfs()`;
   - invalid `n_proj` (too large, wrong length, unknown names).

## 8. Files to add or change

| File | Content |
|---|---|
| `R/pop_sfs.R` | `pop_sfs()` (methods), `pop_sfs_preview()`, `autoplot` of the preview |
| `R/pairwise_pop_sfs.R` | `pairwise_pop_sfs()` |
| `R/nwise_pop_sfs.R` | `nwise_pop_sfs()` |
| `R/sfs_helpers.R` | `.sfs_compute()`, argument checks (`n_proj`, `ancestral`), `.sfs_fold()`, `new_sfs()`, `print`/`tidy`/`autoplot` methods |
| `R/write_sfs.R` | `write_sfs()`, `gt_as_easysfs()` |
| `src/sfs_config_counts.cpp` | stage 1 kernel |
| `src/sfs_project.cpp` | stage 2 projection (and fold, if in C++) |
| `tests/testthat/helper_sfs_reference.R` | R port of the easySFS/dadi logic |
| `tests/testthat/test_pop_sfs.R`, `test_pairwise_pop_sfs.R`, `test_nwise_pop_sfs.R`, `test_easysfs.R` | tests |
| `data-raw/create_easysfs_test_data.R`, `tests/testthat/testdata/easysfs/` | fixtures |
| `_pkgdown.yml` | `write_sfs`, `tidy.sfs` and `autoplot.sfs` in a new "Site frequency spectrum" section, or under "Population" |
| `NEWS.md`, `NAMESPACE`, `RcppExports` | the usual |
| `vignettes/articles/` | a short article: preview, choose projections, compute, export to fsc/dadi, with a note on thinning RAD data |

## 9. Phases

1. **Core.** Stage 1 and stage 2, `pop_sfs()`, `pairwise_pop_sfs()`,
   `nwise_pop_sfs()`, `pop_sfs_preview()`, the `sfs` class (`print`, `tidy`),
   and the reference and invariant tests.
2. **easySFS parity and export.** `total_length`, `write_sfs()`,
   `gt_as_easysfs()`, the easySFS fixtures and tests, the MSFS order check,
   and `autoplot()` methods.
3. **Extras.**
   - `ancestral` polarisation.
   - A block bootstrap of spectra for fsc/dadi confidence intervals. Stage 1
     would keep configuration counts per block, using the `block_id` of
     `.jackknife_block_ids()`. Replicates are then weighted sums of
     per-block spectra, without reading the genotypes again; this replaces
     easySFS's `--resample-vcf`.
   - Statistics from an `sfs` (theta_pi, theta_W, Tajima's D, 2D Fst).

## 10. Decisions for the maintainers

Decisions 1, 2 and 4 have been taken; see [Status](#status-branch-sfs).

1. **The default `n_proj`.** Decided: the projection that maximises the
   number of segregating sites.
2. **The return type.** Decided: a tidy tibble by default.
3. **One SNP per RAD locus.** Still open. The plan leaves this to loci
   selection and the documentation, rather than an argument of the SFS
   functions.
4. **The MSFS cell order.** Decided: the order of the populations by default,
   with the `fastest` option.
