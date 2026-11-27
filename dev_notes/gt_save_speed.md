# Saving and loading .gt files with many loci

How long `gt_save()`, `gen_tibble()` (which saves the *.gt* file) and
`gt_load()` take with tens of millions of loci, and why: the two allele
columns of the loci table. A prototype that saves the alleles as codes made
loading 4.2x and saving 2.9x faster on 91M loci, but it was **not adopted**
(Andrea, 20 Nov 2026: not worth the extra complexity; with this many loci,
`names_as_int = TRUE` is needed anyway). It is kept for a possible future as
a patch, `dev_notes/gt_save_speed/alleles_as_codes.patch` (commit
`332f5405` of the deleted local branch `gt_save_speed`; on 20 Nov it
applied cleanly to `bed_reader`, with `git am`). It touches `save_gt_file()`
(with the new `compact_loci_alleles()` and `restore_loci_alleles()`) in
`R/gt_save.R`, `gt_load()` in `R/gt_load.R`, `encode_alleles_cpp()` and
`decode_alleles_cpp()` in `src/strings_in_cpp.cpp`, the docs of `gt_save()`
and a test. Work of 20 Nov 2026 on manica01.

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| Keep the current format (character alleles, gzip level 1); do not merge the prototype | save the alleles as codes (the prototype) | the gains (17 -> 4 s to load, 21 -> 7 s to save on 91M loci) are not worth the extra code and the change of format; with tens of millions of loci, `names_as_int = TRUE` is needed anyway, as millions of strings slow down everything else too | 2026-11-20 | Andrea |
| (prototype) Save the alleles as codes into their distinct strings | no compression; zstd (R 4.5 has `zstdfile()`); alleles as factor columns | on 91M loci, uncompressed saves took 11.7 s (gzip level 1: 18.4) and loads 13.6 s (16.3): the cost is per string, not the compression. zstd is optional in R builds, so files would not be portable. Factor columns would be read by older versions as factors, and e.g. `ifelse()` on them gives codes: wrong results without an error | 2026-11-20 | Claude |
| (prototype) Remove the allele columns from the saved table, and keep the codes in an attribute of the table | keep the columns with codes in them | older versions then stop with an error (no `allele_ref` column) rather than use codes as alleles | 2026-11-20 | Claude, prototype |
| (prototype) Distinct alleles by the address of their CHARSXP | `match()`/`factor()` on the strings | R keeps one CHARSXP per string and encoding, so decoding gives back the very same strings (two encodings of the same text stay two levels); a pass of ~1 ns per locus, with no strings compared | 2026-11-20 | Claude |
| (prototype) Codes as raw with up to 255 distinct alleles, integer beyond | always integer | a quarter of the size to compress; datasets with indels can have millions of distinct alleles | 2026-11-20 | Claude |

## How it works

**Current.** `save_gt_file()` saves the gen_tibble with `saveRDS()` through
`gzfile(compression = 1)`; `gt_load()` reads it with `readRDS()`. Each
element of a character column is written (and read, creating or finding its
string in R's cache) one by one: with 91M loci the two allele columns take
two thirds of the time.

**Prototype (the patch).** `save_gt_file()` calls `compact_loci_alleles()` on the copy it saves: if the
loci table has character `allele_ref` and `allele_alt`, `encode_alleles_cpp()`
gives a code per locus for each (1-based, into one vector of levels shared by
both; NA is a level), the two columns are dropped and the attribute
`compact_alleles` of the table holds `columns` (the names of all columns,
in order), `allele_ref`, `allele_alt` (the codes) and `levels`. The
gen_tibble in memory is not modified (the saved table shares its other
columns). `gt_load()` calls `restore_loci_alleles()` right after `readRDS()`:
without the attribute (files of older versions) nothing changes; otherwise
the columns are decoded (`decode_alleles_cpp()`, which checks every code) and
put back in their place.

Tests (`test_gen_tibble_save_load.R`, "the alleles are saved as codes, and
load back identical"): a table with NA alleles and one with 300 distinct
alleles (integer codes); the saved table has no allele columns; the
gen_tibble in memory is unchanged; `gt_load()` gives an identical loci
table. The old-format tests (`old bigsnp gen_tibbles`) load files of older
versions. The test fails on `bed_reader`. The full suite passes (740
tests, on an -O2 install of `332f5405`).

## Numbers

manica01, 20 Nov 2026, the gen_tibble of the baboon data (486 individuals x
91,162,139 loci, `names_as_int = TRUE`), the same table saved in the old
format (gzip level 1, character alleles) and in the new one; one P-core
(`taskset -c 2`), -O2 installs (`bed_reader` at `4c7cbdeb` and
`gt_save_speed`), 3 runs each in A B B A A B order, quiet machine. Seconds,
median (range):

| format | `gt_load()` | saving the *.gt* (`save_gt_file()`) | file |
|---|---|---|---|
| old | 17.21 (17.14-17.26) | 20.89 (20.88-21.01) | 305 MB |
| new | 4.13 (4.12-4.18) | 7.25 (7.20-7.28) | 225 MB |

- The loci tables loaded from the two files are `identical()`.
- On a synthetic table of the same size (`bed_reader/scripts/save_split.R`):
  serialising to memory alone takes 7.5 s; uncompressed, saving takes 11.7 s
  and loading 13.6 s; with factor alleles, 6.9 and 3.3 s (+0.8 s for
  `as.character()`).
- Expected on all 166M loci of the baboon data: ~13 s to save and ~7.5 s to
  load (from 38 and 31 s).

To reproduce: `dev_notes/bed_reader/scripts/gt_setup.R <lib> <bed>
<prefix>` makes the gen_tibble and both files; `gt_time.R <lib> <gt file>
<out dir>` times `gt_load()` and `save_gt_file()` (logs in
`dev_notes/bed_reader/results/2026-11-20/logs_gt_time.txt`).

## Open questions

- **If it is ever revived:** *.gt* files saved with it can not be loaded by
  older versions (they stop with an error when they use the alleles), so it
  is best done with a new minor version. Apply the patch, then NEWS (the
  item on saving the *.gt* file faster), the version, the full suite, and
  re-time with `gt_setup.R` and `gt_time.R`.
- With character locus names (not `names_as_int`), each of millions of
  unique names is still saved and loaded one by one (and loading creates a
  string for each); codes do not help there.
