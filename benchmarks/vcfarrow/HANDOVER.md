# Handover: run the VCFArrow benchmark on a desktop

Written on 2026-11-08 at the end of a Claude Code session on Andrea's laptop,
which was too small and slow for the full run. The task is to run this
harness on a desktop and report to the user. This file is self-contained.

Status: run on the desktop on 2026-11-09 (24 cores, 126 GB). The results
are in `README.md` ("Results (desktop ...)") and in `results/desktop_24c*`.

## What to do

```sh
cd <your tidypopgen clone>
git fetch origin
git checkout vcf_fast
cd benchmarks/vcfarrow
./setup_libs.sh     # installs v1.0.0 and vcf_fast into lib/
./run_all.sh        # 1.0.0 vs vcf_fast at 2.3 M, 4.6 M and 18.6 M variants
```

Check the output of `setup_libs.sh`. Each version must report "non-inlined
Rcpp indexing calls (0 = no -O0 code): 0", then "loaded 1.0.0" for v1.0.0
and "loaded 1.0.0.9000" for vcf_fast.

`run_all.sh`:
- downloads the VCFArrow subsample (two small files) into `data/`;
- generates the VCFs (3.4, 6.8 and 27 GB, so about 40 GB of disk is
  needed);
- runs every step for both versions, then writes `results/summary.md`.

Settings are environment variables (see the top of `run_all.sh`). On a
desktop:
- raise `THREADS` (the comparison used 4 and 15);
- raise `MEMORY_MAX` (12G per step by default; the comparison had 50 GB).

For example: `THREADS=8 MEMORY_MAX=32G ./run_all.sh`. For a first look, run
`SIZES="2322460" ./run_all.sh`.

Before timing, check that the machine is idle (`uptime`). Do not time
anything while other heavy jobs run.

## Requirements

- R (4.6.1 on the laptop), with:
  - the dependencies of tidypopgen 1.0.0;
  - `dplyr`, `data.table`, `ggplot2`, `adegenet`, `jsonlite`, `tibble`.
- python3, standard library only.
- curl.
- A systemd user session: `bench.py` runs each R step in
  `systemd-run --user --scope` with `MemoryMax` and `MemorySwapMax=0`.

## What the user wants

VCFArrow (https://github.com/legalLab/VCFArrow, `comparison/tool_comparison.md`)
benchmarked 17 functions with tidypopgen 0.4.4 on an 18.6 M-variant x 56-sample
DiscoSnpRad VCF. 0.4.4 is known to be slow, and its published numbers appear
in the summary. The user wants:
- 1.0.0 against `vcf_fast`, step by step;
- next to VCFArrow's published numbers;
- with an account of where tidypopgen is still slow, and why.

Report:
- per step: time and peak PSS for both versions, and the speed-up;
- whether the results agree (`summary.md` checks this);
- the estimates at 18.6 M variants next to the published tidypopgen 0.4.4
  and VCFArrow numbers.

## The branch `vcf_fast`

Its changes since 1.0.0 are in NEWS.md:
- `gt_as_vcf()` writes the header in linear time (it was quadratic in the
  number of chromosomes, and this data has one chromosome per RAD locus),
  and writes missing chromosomes as ".".
- `gt_as_genlight()` works in chunks.
- `select_loci()` resolves calls (`sample()`, `which()`) without tidyselect,
  and gains `.engine`.
- `.gt` files are saved at gzip level 1.
- The VCF reader:
  - reads in a single pass, in parallel with `n_cores` (BGZF blocks and
    genotypes);
  - looks up chromosome names without an allocation per chromosome, and
    builds the chromosome factor faster;
  - gives clearer errors.
- The last commit adds this folder, `benchmarks/`, which is not part of the
  package. Whether to keep it on the branch is for the user to decide before
  merging.

## Earlier results (laptop)

See the end of `README.md` and `results/laptop_test_200k_summary.md`. At
200 k variants, the two versions give identical results on all steps. Write
VCF takes over 2 min in 1.0.0 and 1.1 s in vcf_fast.

Things to check at the larger sizes:
- `to_genlight` memory: at 200 k variants vcf_fast peaked at 0.66 GB,
  against 0.56 GB for 1.0.0. The chunks are meant to cap memory on large
  data.
- `merge` (`rbind()`).
- `subsample_loci_bare`: a bare index vector misses the fast path, which is
  expected (tidyselect deprecates it).

## Known remaining costs (VCF reading, 2.3 M variants)

- R creating 4 M new strings for marker and chromosome names, about 3.4 s.
  This is intrinsic to R, and `names_as_int = TRUE` avoids it.
- Saving the `.gt` file, 2.8 s. The user does not want this step changed.
- libdeflate was set aside by the user, but is to be kept in mind. It would
  save about 3 s on 1 core. See `libdeflate_proto/README.md` and the
  comment above `VcfLineReader::read_bgzf()` in `src/vcf_parser_fast.cpp`.

## The user's working rules

- Never time -O0 builds: `load_all()` compiles the C++ at -O0. Time
  `R CMD INSTALL` builds (as `setup_libs.sh` does), or use
  `pkgbuild::clean_dll()`, `compile_dll(debug = FALSE)` and
  `load_all(compile = FALSE)`.
- Commits: one commit per fix, on a named branch, holding the code, docs,
  RcppExports, tests and NEWS together.
- Tests: few and fast, never performance tests.
- Match the existing style, and comment densely.
- Never push unless asked.
