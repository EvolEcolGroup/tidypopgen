# Plan: integer chromosomes, and sidecar files for loci names

Two additions on top of `vcf_fast`:

1. `chromosomes_as_int`: store the chromosome of each locus as an integer
   code, as `names_as_int` already does for locus names.
2. Sidecar files: when names or chromosomes are stored as integers, the
   readers write the original names, chromosomes and positions to plain text
   files next to the backing file. New functions read them back, for all
   loci or for a few (e.g. the hits of a selection scan).

Status: **phases 1 to 3 are done** on the branch `vcf_names_tsv`; phase 4
is not started. See [loci_sidecar_handover.md](loci_sidecar_handover.md) for
what was built, how it differs from this plan, the benchmarks, and
explorations of phase 4. A position mismatch is an error, with
`check_positions = FALSE` as the way out (section 8).

## Why

The DiscoSnpRad data of the VCFArrow benchmark has 18.6 M variants on about
13.8 M "chromosomes" (one per RAD path). With `names_as_int = TRUE` the
locus names are gone from memory, but the chromosome factor still has 13.8 M
levels. Each level is an R string of about 80 bytes, so the levels take about
1 GB. R also has to create and sort those strings when the VCF is read, and
write them to the `.gt` file. With integer codes, the chromosome column is
74 MB of integers.

`names_as_int` also loses the VCF IDs for good. After a scan, there is no way
to find out which variants the hits are. The sidecars keep that information
on disk, where it costs no memory, and make the integer options lossless.

## 1. Integer chromosomes

### What the codes mean

Chosen: the rule of `cast_chromosome_to_int()`, which
`gt_as_plink(chromosomes_as_int = TRUE)`, `gt_order_loci()`,
`loci_ld_clump()` and `gt_pca_autoSVD()` already use.

- A name that is a number, optionally after a `chr`/`chromosome` prefix
  (`"7"`, `"chr7"`, `"Chr07"`), keeps that number as its code.
- Other names (`"chrX"`, `"scaffold_17"`, `"SNP_higher_path_42"`) are
  numbered after the largest number, in order of first appearance.

So chromosome 10 is coded 10 whatever order the file has, and Stacks catalog
loci keep their catalog number. Two changes to the R rule are needed to make
codes reversible:

- If two names give the same number (`"1"` and `"chr1"`), the second one is
  treated as an "other" name.
- A number above `.Machine$integer.max` is treated as an "other" name.

The alternative is simpler: number every chromosome in order of first
appearance, as `names_as_int` does for loci. Its drawback is that a file
sorted as 1, 10, 11, ..., 2 gives chromosome "10" the code 2, which is easy
to misread.

`"."` (missing) stays `NA`, as it does now.

### Where it is available

The same readers as `names_as_int`: `gen_tibble()` for VCF with
`parser = "cpp"` (phase 2), `gt_merge_vcfs()` and `gt_concat_vcfs()`
(phase 3), and PLINK bed / bigSNP rds (phase 4, optional). The default is
`FALSE` everywhere, including in merge and concat: chromosome names cost
little for most data.

### C++ (VCF parser)

The parser already gives each locus the index of its chromosome in order of
first appearance, with the names held once in a `ChromosomeTable`. With
`chromosomes_as_int`, at the end of the parse it will:

- work out the code of each of the K chromosomes;
- remap the index vector in place;
- write the chromosome sidecar straight from the `ChromosomeTable`;
- return no `chromosome_names`, so R never creates the K strings.

The code rule will live in one C++ function, used by the parser and, through
a wrapper on a character vector, by merge and concat (which already have the
names in R).

Missing IDs are filled as `<chromosome>_<position>` in R today, which needs
the chromosome name. With integer chromosomes, this filling moves into C++,
where CHROM is at hand. The C++ code will follow `fill_missing_loci_names()`
exactly, including `"NA_<pos>"` for a missing chromosome.

### Keeping the column an integer

`validate_loci()` turns any chromosome column that is not a factor into a
factor. With 13.8 M codes, that would recreate every string. Tests build
matrix loci with `chromosome = 1L` and expect a factor, so the default can't
simply change. Instead:

- `validate_loci(..., chromosomes_as_int = FALSE)`: the readers pass `TRUE`.
- `show_loci<-` keeps integer chromosomes if the current table has them.
  `gt_update_backingfile()` assigns `show_loci()` this way, so without this
  the codes would turn back into a factor.

Places that call `as.factor()` or `levels()` on the chromosome column need an
integer path:

| Function | File | Change |
|---|---|---|
| `any_unsorted_by_chromosome()` | `is_loci_table_ordered.R` | pass codes to the C++ helper with `level_group = seq_len(max)`, or `match(x, unique(x))` when the codes are sparse or include 0 |
| `chromosomes_not_adjacent()` | `is_loci_table_ordered.R` | same |
| duplicate check | `validate_loci()` | same |
| `windows_indiv_roh()`, `gt_roh_*()` | ROH files | keep integer codes; output has the same type as the input |
| `vcf_write_header()` | `gt_as_vcf.R` | check that it handles integer chromosomes |
| `gt_extract_f2()`, `gt_as_genlight()` | | check, then test |

The window functions, the jackknife helpers, `find_duplicated_loci()` and
the `cast_chromosome_to_int()` callers already work with integers.

### `rbind()` and `rbind_dry_run()`

Matching by position compares chromosome labels. Integer codes from two
datasets are not comparable: code 23 could be `chrX` in one file and
`scaffold_1` in the other.

`chromosome_codes()` in `rbind_dry_run.R` will take the labels of an integer
column from its chromosome sidecar, read only for the codes in use. Any mix
of factor and integer then matches by name. The merged object keeps the
reference's representation, and its sidecars are copied (see section 2). If
a sidecar is missing, it stops with a clear error.

## 2. The sidecar files

### Format (chosen)

Plain, uncompressed, tab-separated text:

- a header line;
- one line per item;
- the first column is the integer key, strictly increasing;
- no quoting (VCF fields cannot contain tabs or newlines);
- `"."` for a missing value;
- `\n` line endings on every platform.

Locus sidecar, `<backingfile>.loci.tsv`, written when `names_as_int = TRUE`:

```
int_name	name	chromosome	position
1	rs3094315	1	752566
2	rs12124819	1	776546
3	1_832918	1	832918
```

`int_name` is the value in `show_loci(x)$name`. `name` is what
`names_as_int = FALSE` would have given, with a missing ID filled. The
chromosome is the original name, so the file can be read on its own, like a
`.map` file.

Chromosome sidecar, `<backingfile>.chromosomes.tsv`, written when
`chromosomes_as_int = TRUE`:

```
int_chromosome	chromosome
1	chr1
2	chr2
23	chrX
24	scaffold_17
```

The files are readable as is: `data.table::fread(path, quote = "")`,
`readr::read_tsv(path, quote = "")`, `cut`, `grep` or `awk`. Joining
`show_loci(x)` on `name = int_name` gives the original names in any tool.

**Why two files?** A locus key only exists when names are integers. With
integer chromosomes and character names, only the chromosome table can map
the codes back. The chromosome table is also what turns the codes in a
windows or ROH result back into names. Getting it from the locus file would
mean scanning the whole file every time.

**Why the explicit key column?** The line number would do for the locus
file, but the chromosome codes have gaps (1..22, then 23 for `chrX`). With a
key column, both files have the same layout, one C++ reader serves both, and
a file can be cut down to some loci (see Lifecycle). The cost is about 9
bytes per line.

### Indexing

To fetch a few keys, a sparse index records the key and byte offset of
every 1024th line, about 18 k entries (200 KB) for 18.6 M loci. A lookup
binary-searches this index, then reads at most 1024 lines (about 50 KB). A
hundred lookups take milliseconds. Sorted keys are read in one forward pass,
so restoring every locus streams the file once.

The index is built by one `memchr` scan of the file, which takes about
0.3 s per GB when the file is in the page cache and a few seconds from disk.
It is kept for the session in a package environment, keyed on the
normalised path, size and mtime. The same scan checks that the keys strictly
increase, which catches truncated or edited files. Adding no `.idx` file
keeps the file set simple. If the first lookup turns out too slow, writing
the index at creation is a later option.

### Alternatives considered

| Format | For | Against |
|---|---|---|
| **Plain TSV + sparse index (chosen)** | Simple; readable by any tool; fast to write (about 1 s per GB); random access | Size: on the VCFArrow data at 2.3 M variants, 44 bytes per locus (103 MB) and 32 bytes per chromosome (55 MB for 1.7 M chromosomes). That is more than the 130 MB `.bk` with 56 samples; with hundreds of samples it is small beside the `.bk` |
| bgzip TSV + virtual-offset index | About 5x smaller; still `zcat`-able; we already read BGZF | Compression time on write (parallel, or libdeflate); a BGZF writer and virtual offsets to maintain; `fread` needs R.utils for `.gz`. A good upgrade later, with the same layout |
| tabix (`.tbi`) | Standard tooling; region queries | Indexes chromosome:position, not our keys; needs sorted input; we'd write `.tbi` ourselves |
| Fixed-width TSV | Seek by arithmetic, no index | Width unknown until the end of the parse; wasteful with long IDs |
| Binary columns, fst, Parquet, SQLite | Fastest lookups; SQLite also looks up by name | Not human-readable, several files, or a new dependency |

### Lifecycle

**Location.** Each sidecar sits next to its backing file:

```
foo.bk
foo.rds
foo.gt
foo.loci.tsv
foo.chromosomes.tsv
```

**Finding it.** The path comes from the FBM's backing file, so the
gen_tibble needs no new attribute. Several `.gt` files that share a `.bk`
share its sidecars. After moving all the files together,
`gt_load(reattach_to = )` finds the sidecars too.

**Files are never edited.** A function that writes a new backing file also
writes new sidecars:

- `gt_update_backingfile()`, and so `gt_order_loci()`;
- `rbind()`, from the reference's sidecars.

The new sidecars keep only the keys that are still in the loci table, in
streaming C++ with a bitset of keys, so memory is N/8 bytes. Integer names
keep their values, so the keys stay valid. If the old sidecar is missing,
the function warns and the new gen_tibble has none.

**Checking.** Only the file name ties a sidecar to its data, and every
locus sidecar has keys 1..n, so a wrong file still has a line for key 123.
If the wrong file sits next to the `.bk`, a lookup would return the name of
a different variant without any sign of a problem. This can happen when:

- files from two runs are copied by hand into one folder;
- the files are moved without their sidecars, and an older sidecar of the
  same name is already at the new place;
- a bug in our subset code writes misaligned lines.

The position on each sidecar line is the check. When names are integers,
`gt_lookup_loci()` and `gt_restore_loci()` compare the sidecar's position of
each locus with the `position` column of the table. Two files almost never
agree on the positions of the same keys, so a mismatch means that the
sidecar belongs to other data.

There is one false alarm. If the user changes positions on purpose (e.g.
`show_loci(x)$position <- lifted` after a liftover to a new assembly), each
variant keeps its key and the sidecar still gives the right name, but the
positions no longer match. Liftover is rare here: the package has no
liftover function, and de novo RAD data has no reference to lift.

When only the chromosomes are integers, there is nothing to compare: the
chromosome sidecar has no positions. The check is then that every code is
in the file.

Section 8 has the options for a mismatch.

**Cleanup.** The readers already delete the `.bk` if they fail. They will
delete the sidecars as well.

**Who writes what:**

| Producer | `.loci.tsv` | `.chromosomes.tsv` |
|---|---|---|
| `gen_tibble()` VCF cpp | streamed by the parser in the in-order metadata thread, buffered | written from the `ChromosomeTable` at the end |
| `gt_concat_vcfs()` | a temporary file per file in the scan pass, joined with offset keys | written from R |
| `gt_merge_vcfs()` | see below | written from R |
| `gt_update_backingfile()`, `rbind()` | subset | subset |
| bed / rds (phase 4) | streamed from the `.bim` in C++ | |

`gt_merge_vcfs()` cannot keep a temporary file per input: the SGDP merge has
hundreds of files, each with every site. Instead:

1. Each scan writes the file's loci to one temporary file.
2. After `merge_sites_add_file()`, the lines of `new_loci` (the sites first
   found in this file) are appended to a running "site names" file, keyed by
   site id. The per-file temporary file is then deleted.
3. At the end, the kept sites are gathered from the running file in merged
   order and rekeyed 1..n.

This matches the current rule that a locus takes its name from the first
file where it is found. Temporary disk use stays about the size of the final
sidecar plus one file's loci.

## 3. User functions

```r
# bring back: the same gen_tibble, with character names and factor
# chromosomes in its loci table (the result of reading with the integer
# options off). Meant for subsets, or for data small enough to hold the
# strings
gt_restore_loci(.x, names = TRUE, chromosomes = TRUE)

# fetch: give back both the original names and the original chromosomes
# for a subset of the loci table, without changing .x
gt_lookup_loci(.x, loci)

# translate chromosome codes (e.g. the chromosome column of
# windows_pairwise_pop_fst() or gt_roh_window()) into names
gt_lookup_chromosomes(.x, chromosomes)
```

### `gt_lookup_loci()`: both names and chromosomes from a subset

`loci` is a subset of the loci table: any rows of `show_loci(.x)`, in any
order, with any extra columns (e.g. the statistics of a scan). The function
returns the same table, with the same rows in the same order and every
other column untouched, where:

- an integer `name` column becomes the original names, read from the locus
  sidecar with the integer names as keys;
- an integer `chromosome` column becomes the original chromosome names, read
  from the chromosome sidecar with the codes as keys.

Both columns are character in the result. A column that already holds names
(character or factor) is left as it is. So the same call works whichever
integer options the data was read with, and returns data read without them
unchanged. A sidecar is only needed for a column that has to be translated.

`loci` can also be row numbers of `show_loci(.x)`, as given by `which()`.
This is the same as passing `show_loci(.x)[loci, ]`.

The keys are taken from the table, not from row numbers, so a subset stays
valid after `.x` is filtered or reordered, as long as its loci are still in
`.x`. A key that is not in the sidecar (a locus that was dropped from `.x`,
or a table from other data) is an error that lists the first few such keys.

For example, after a selection scan:

```r
x <- gen_tibble(vcf, names_as_int = TRUE, chromosomes_as_int = TRUE)
hits <- show_loci(x) %>%
  mutate(p = p_values) %>% # one p-value per locus
  filter(p < 1e-8)
gt_lookup_loci(x, hits)
#> # A tibble: 12 x 8
#>   big_index name    chromosome              position ... p
#>       <int> <chr>   <chr>                      <int>     <dbl>
#> 1      4021 9999561 SNP_higher_path_9999561       77 3.1e-09
#> ...

# the same, from row numbers
gt_lookup_loci(x, which(p_values < 1e-8))
```

Internally, `gt_restore_loci()` takes the same path for names, on the whole
table. For chromosomes it builds a factor from the codes instead, so that no
string is created per locus.

Each function errors clearly if it needs a sidecar that does not exist (data
read before this feature, or a sidecar that was deleted). The doc pages of
`gen_tibble()`, `gt_save()` and `gt_load()` will describe the sidecar files,
and `gt_save()` will list them in its messages.

Later (phase 4):

- `gt_find_loci(.x, names)`: the reverse lookup. It streams the locus
  sidecar with a hash set of the query names, so that loci can be selected
  by rsID when names are integers.
- `gt_as_vcf(restore = TRUE)`: writes the original IDs and CHROM from the
  sidecars, chunk by chunk, so memory stays bounded. The `##contig` lines
  are already written in chunks. `gt_as_plink()` would get the same option.

## 4. Internals

`src/loci_sidecar.cpp` will hold:

- `SidecarWriter`: a buffered `FILE*` in binary mode, used by the parser.
- `sidecar_index(path, block)`: the scan that builds the sparse index and
  checks the keys.
- `sidecar_fetch(path, keys, index, columns)`: the requested columns, for
  unsorted keys, with `NA` for missing keys. Strings are made with
  `Rf_mkCharLen`, as the parser does.
- `sidecar_subset(path_in, path_out, keys)`.
- `sidecar_gather(path_in, keys_in_order, path_out)` for merge (phase 3).
- `sidecar_find(path, column, values)` (phase 4).

All offsets are 64-bit, so files over 2 GB work on Windows too.

`R/loci_sidecar.R` will hold the paths, the session cache, the position
check, and the three exported functions. Restoring chromosomes builds the
factor from the codes, with levels sorted by `order_as_levels()`, as
`chromosome_index_to_factor()` does. No string is created per locus.

## 5. Phases (one commit each, with code, docs, RcppExports, tests and NEWS)

Work on the branch `vcf_names_tsv`, made from `vcf_fast` (which is PR #226,
under review).

1. **Locus sidecar for `names_as_int`, VCF.**
   - The C++ toolkit and the writer in `vcf_to_bk()`.
   - `gt_restore_loci()` (names) and `gt_lookup_loci()`.
   - Subset copies in `gt_update_backingfile()` and `rbind()`.
2. **`chromosomes_as_int`, VCF.**
   - The codes and the chromosome sidecar in C++.
   - `validate_loci()` and `show_loci<-`.
   - The integer paths in the table of section 1.
   - Chromosome matching in `rbind()`.
   - Restoring chromosomes, and `gt_lookup_chromosomes()`.
3. **Merge and concat.**
   - `chromosomes_as_int` and both sidecars, with the temporary files
     described in section 2.
   - Update the big data article and the merge/concat benchmark notes.
4. **Optional.**
   - `gt_as_vcf()` / `gt_as_plink()` restoring names.
   - bed/rds sidecars and `chromosomes_as_int`.
   - `gt_find_loci()`.

## 6. Tests (few and fast)

- **Round trip.** For the anolis VCF and a small VCF with missing IDs and a
  `"."` chromosome:
  `gen_tibble(vcf, names_as_int = TRUE, chromosomes_as_int = TRUE) %>% gt_restore_loci() %>% show_loci()`
  is identical to `show_loci(gen_tibble(vcf))`. The same holds with
  `n_cores = 2`, and for merge and concat.
- **Lookups.** On subsets and reordered loci, after `gt_order_loci()` and
  after `rbind()`. Use a tiny internal block size to cross block edges.
  `gt_lookup_loci()` on a filtered, reordered table with an extra column
  gives back the same rows and columns, with names and chromosomes
  translated. It returns a table without integer columns unchanged, and the
  row-number form equals the table form.
- **Errors.** A sidecar with an edited position, a missing sidecar, and
  unsorted keys.
- **Codes.** Numbers keep their value. `chrX` is coded after the largest
  number. `"1"` and `"chr1"` together, and an overflowing number.
- **Integer chromosome column.** `is_loci_table_ordered()`, one windows
  function, and `rbind()` of factor x integer and integer x integer.

## 7. Performance checks

All timings use -O2 builds, as `CLAUDE.md` requires:

- `pkgbuild::clean_dll()`;
- `compile_dll(debug = FALSE)`;
- `load_all(compile = FALSE)`;
- the baseline installed into its own library.

At 2.3 M variants of the VCFArrow data (`benchmarks/vcfarrow`, with a new
`read_vcf_int_chr` step), measure:

- the time and peak memory of reading, with and without
  `chromosomes_as_int`;
- the cost of writing the sidecars (expected under 5% of the read);
- one lookup of 100 loci, cold and warm;
- `gt_restore_loci()` on all loci.

## 8. Decisions

Agreed:

1. Integer chromosomes: numbers keep their value, and other names are
   numbered after the largest number.
2. Sidecar format: plain TSV with an explicit key column and a session index.
3. Sidecars are always written when an integer option is on. There is no
   argument to turn them off.
4. File names `foo.loci.tsv` and `foo.chromosomes.tsv`. Function names
   `gt_restore_loci()`, `gt_lookup_loci()` and `gt_lookup_chromosomes()`.
   `gt_lookup_loci()` takes a subset of the loci table and gives back both
   names and chromosomes (section 3).
5. Branch `vcf_names_tsv`, made from `vcf_fast`.
6. A sidecar's positions that do not match the table are an error, with an
   off switch (the first option below, chosen by the user). The options
   were (see "Checking" in section 2):

- **Error, with an off switch (recommended).** `gt_lookup_loci()` and
  `gt_restore_loci()` stop with an error saying that the sidecar does not
  belong to this data. An argument `check_positions = TRUE` can be set to
  `FALSE` after a deliberate change of positions, and the error message
  says so. A wrong file can never give wrong names silently.
- **Warning.** Return the names anyway and warn how many positions differ.
  A liftover would not need an argument, but a warning is easy to miss,
  and the names from a wrong file would then be used.
