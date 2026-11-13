# Draft issue for LEA: the seed of `R_sNMF()` is read and written as a `long long`

**Status:** draft, not posted (2026-11-13).

## Why we need it

`gt_snmf()` wraps `LEA::snmf()`, and `test_gt_snmf.R` runs whenever LEA is
installed. Valgrind then reports, in LEA's C code:

- 5 errors `Invalid write of size 8` at `R_sNMF.c:64`;
- 10 errors "uninitialised value" in `print_summary_snmf()`.

These are the only valgrind errors that come from a package's own code, on
R-hub's `valgrind` platform. Both runs report them in the same numbers:

- `nov_fixes` at `68e3a315`, run 43477155061;
- the 1.0.0 submission, run 42381223691.

The other errors come from R itself (the regex library), from libc, from
OpenMP's thread-local storage (the "possibly lost" blocks in
`sfs_config_counts_cpp`), and from the udunits, fontconfig and glib
libraries.

The bug is in LEA, whatever tidypopgen passes to it, so there is nothing to
fix in tidypopgen. If CRAN's valgrind check flags it, point to this issue.

## Where to post it

- LEA's GitHub repository (`bcm-uga/LEA`) has its issues turned off.
- Post the issue on the Bioconductor support site
  (<https://support.bioconductor.org>, tag `LEA`), or email the maintainer
  (Olivier François, `cre` in LEA's DESCRIPTION).
- A pull request on `bcm-uga/LEA` with the fix below is also possible.

The text below is ready to paste.

---

**Title:** `R_sNMF()` and `R_LFMM()` declare the seed as `long long *`, but R
passes an integer (an 8-byte write into a 4-byte integer)

In LEA 3.24.0, and in the current devel (3.25.2), `snmf()` passes the seed to
C as an R integer:

```r
# R/sNMF.R
resC = .C("R_sNMF",
    ...
    s = as.integer(s),
    ...
```

But `R_sNMF()` declares it as a 64-bit integer. It reads 8 bytes from it,
and writes 8 bytes back:

```c
// src/R_sNMF.c
void R_sNMF(..., int *R_iteration, long long *R_seed, ...)
...
        param->seed = *R_seed;   /* lines 44 and 53: 8 bytes read from a 4-byte int */
...
        *R_seed = param->seed;   /* line 64: 8 bytes written into a 4-byte int */
```

`.C()` passes an integer vector as `int *`, so these accesses go 4 bytes past
the end of the data of the vector.

Under valgrind, each call of `snmf()` reports the write. This was on R-hub's
`valgrind` platform, R-devel 2026-11-12 r90881, with LEA 3.24.0:

```
Invalid write of size 8
   at 0x...: R_sNMF (R_sNMF.c:64)
   by 0x...: do_dotCode (dotcode.c:2035)
 Address 0x... is 4,984 bytes inside a block of size 7,960 alloc'd
```

The upper 32 bits of `param->seed` are not initialised either. Valgrind
reports it when the summary of the options is printed:

```
Use of uninitialised value of size 8
   at 0x...: _itoa_word (in /usr/lib64/libc.so.6)
   ...
   by 0x...: print_summary_snmf (print_snmf.c:100)
   by 0x...: sNMF (sNMF.c:77)
   by 0x...: R_sNMF (R_sNMF.c:60)
```

The printed seed shows the stray bits. In a single call of `snmf()` with
`seed = 1`:

- `create.dataset` prints `-s (seed random init) 1`, because
  `R_createDataSet()` takes an `int *`;
- the sNMF summary prints `-s (seed random init) 12884901889`, that is
  3 * 2^32 + 1. In another call it printed `4294967297`, that is 2^32 + 1.

As far as we can tell, the results are not affected:

- under R, `init_random()` does not seed `rand()`;
- `drand()` uses `unif_rand()`, which `set.seed(s)` in `snmf()` seeds;
- on little-endian machines, `resC$s` reads back the lower 4 bytes, which
  hold `s`.

But the write goes past the end of an R vector, which can corrupt R's heap.
It also makes the valgrind checks fail for every package that runs `snmf()`
in its tests or examples.

`R_LFMM()` has the same mismatch: `src/R_LFMM.c` declares `long long *R_seed`
on line 33 and reads it on line 54, while `lfmm()` passes
`s = as.integer(s)`. It only reads the seed, though.

**A possible fix** is what `R_createDataSet()` already does: take an `int *`,
and convert.

```c
void R_sNMF(..., int *R_iteration, int *R_seed, ...)
...
        param->seed = (long long) *R_seed;
...
        *R_seed = (int) param->seed;
```

Do the same in `R_LFMM()`: `int *R_seed`, and
`param->seed = (long long) *R_seed;`. Update the prototypes in `R_sNMF.h`
and `R_LFMM.h` to match.

We found this while checking tidypopgen
(<https://github.com/EvolEcolGroup/tidypopgen>), whose `gt_snmf()` wraps
`snmf()`.
