# tidypopgen on the VCFArrow benchmark

VCFArrow (https://github.com/legalLab/VCFArrow, `comparison/tool_comparison.md`,
October 2026) compares itself with vcftools, vcfR and tidypopgen 0.4.4 on a
DiscoSnpRad VCF of *Phyllomedusa vaillantii* (18,579,683 variants x 56
samples, 26 GB). Neither the dataset nor the scripts are public; only the
results and a description of each call are. This harness reruns the
tidypopgen side on synthetic data of the same shape, comparing 1.0.0 with the
branch `vcf_fast` (0.4.4 is known to be much slower: its published numbers
are shown next to ours). This folder is not part of the package (it is in
`.Rbuildignore` and excluded from lintr).

## Running

```sh
cd benchmarks/vcfarrow
./setup_libs.sh                    # installs v1.0.0 and vcf_fast in lib/
./run_all.sh                       # full plan, several hours
SIZES="2322460" ./run_all.sh       # first look: 1/8 of the data
```

`setup_libs.sh` installs the committed refs of the repository this folder is
in, each into its own library with `R CMD INSTALL` (never `load_all()`, which
compiles the C++ at -O0), and checks for -O0 code. `run_all.sh` downloads the
VCFArrow subsample (two small files) into `data/` the first time, and
generates the VCF of each size there (3.4 GB, 6.8 GB and 27 GB).

Results: `results/vcfarrow_bench.tsv` (one row per version x size x step) and
`results/summary.md` (1.0.0 vs vcf_fast per size, results check, and
estimates at 18.6 M variants next to the published numbers). The log is in
`logs/run_all_*.log`; each step's R output and memory trace is in
`work/<version>_<n>_t<threads>/`. All of these are git-ignored.

Default plan (`run_all.sh`): 1.0.0 and vcf_fast at 2,322,460 (1/8),
4,644,921 (1/4) and 18,579,683 (all) variants. Settings (environment
variables, see `run_all.sh`): `SIZES`, `OLD` (v1.0.0), `NEW` (vcf_fast),
`OLD_MAX`, `THREADS` (4), `MEMORY_MAX` (12G per step, no swap; the
comparison had 50 GB), `TIMEOUT` (60 min at full size as in the comparison,
scaled down with size, at least 10 min), `STEPS` (`all` or `core`). Outputs
of writers are deleted after each step. Another ref can be compared with
`./setup_libs.sh <ref>` and `OLD=... NEW=<ref> ./run_all.sh`.

Compared with 1.0.0, vcf_fast has (see NEWS.md):
- `gt_as_vcf()` writing the header in linear time (it was quadratic in the
  number of chromosomes, so "write VCF" timed out on this data), and
  writing missing chromosomes as ".";
- `gt_as_genlight()` working in chunks;
- `select_loci()` resolving calls such as `sample()` directly;
- the .gt file saved at gzip level 1;
- a VCF reader that reads in a single pass, in parallel with `n_cores`,
  with faster chromosome handling.

## Files

- `setup_libs.sh`: installs refs with `R CMD INSTALL` into `lib/<label>/`
  (never `load_all()`, which compiles the C++ at -O0) and checks for -O0 code.
  For 0.4.4 use `v0.4.4=5ab2687c9` (CRAN 0.4.4; the `v0.4.4` tag is one
  commit earlier, with the same R/ and src/ but version 0.4.3.9006).
- `make_vcf.R`: synthetic VCF of any size, resampled from the 10,000 x 18
  subsample of the real data shipped with VCFArrow
  (`inst/extdata/vaillantii_discosnp_sub.vcf.gz` and `strata` there,
  downloaded into `data/` by `run_all.sh`): same FORMAT fields
  (GT:DP:PL:AD:HQ), line length (~1.46 kB with 56 samples), missingness
  (~33%), 7% insertions, four lineages (18/18/17 + 3 outgroup), and RAD paths
  as CHROM ("SNP_higher_path_<n>", 1-5 SNPs each, so ~74% of variants start
  a new chromosome: ~13.8 M chromosomes at full size), in descending
  lexicographic order as in the subsample.
- `steps.R`: the tidypopgen call for each function of the comparison (below)
  plus extra steps.
- `run_step.R`: runs steps in a fresh R process; every step but the read
  starts from `gt_load()` of the gen_tibble saved after the read.
- `bench.py`: one process per step, PSS/RSS of the process tree sampled
  every 0.5 s (the comparison reported PSS of all processes), time limit per
  step, memory limit through a systemd user scope.
- `summarise.R`: tables in Markdown.

## The 17 functions, as we read them

The comparison names only some calls (`gt[1:10, ]`, `dplyr::filter()`,
`loci_pi()`, `qc_report_indiv()`, `rbind()`); the others are our best guess.
All start from the unfiltered data. `n` = 4 threads.

| comparison                  | tidypopgen call |
|-----------------------------|-----------------|
| read VCF                    | `gen_tibble(vcf, parser = "cpp", n_cores = n)` |
| set groups                  | `mutate(population = strata$pop[match(id, strata$id)])` |
| extract samples             | `gt[1:10, ]` |
| extract groups              | `filter(gt, population %in% c("GS", "BS", "WA"))` |
| subsample loci (1000)       | `select_loci(gt, all_of(sort(sample.int(count_loci(gt), 1000))))` |
| filter MAF (0.05)           | `select_loci_if(gt, coalesce(loci_maf(genotypes) >= 0.05, FALSE))` |
| filter locus missing (0.2)  | `select_loci_if(gt, loci_missingness(genotypes) <= 0.2)` |
| filter sample missing (0.8) | `filter(gt, indiv_missingness(genotypes) <= 0.8)` |
| filter invariant            | `select_loci_if(gt, coalesce(loci_maf(genotypes) > 0, FALSE))` |
| merge two VCFs              | `rbind(ingroups, outgroup, backingfile = ...)` (two sample subsets) |
| sample statistics           | `qc_report_indiv(gt)` |
| theta / pi                  | `loci_pi(gt)` |
| missing data plot + table   | `indiv_missingness()`, CSV, `ggsave()` of a bar chart |
| write VCF                   | `gt_as_vcf(gt)` |
| -> genlight                 | `gt_as_genlight(gt)` |
| -> PLINK ped                | `gt_as_plink(gt, type = "ped")` |
| -> PLINK bed                | `gt_as_plink(gt, type = "bed")` |

Extra steps (not in the comparison): `subsample_loci_bare` (the same
selection without `all_of()`, which misses the 1.0.0 fast path),
`theta_pi_grouped` (pi per lineage, closer to `vcf_theta()`),
`filter_chain_update` (the VCFArrow vignette filter chain, then
`gt_update_backingfile()`, since VCFArrow filters always write a complete new
object) and `read_vcf_int` (`names_as_int = TRUE`, 1.0.0 only).

## Differences from the published runs

- Synthetic data of the same shape (see `make_vcf.R`), not the real VCF.
- One fresh R process per step, not one session for all functions, so the
  memory of a step is the loaded gen_tibble plus the step, without what
  earlier steps left behind.
- Another machine, a 12 GB limit per step by default, and 4 threads (the
  published 4 and 15 threads gave the same times).
- Both tidypopgen versions drop the insertions (non-ATCG alleles) when
  reading, so they keep ~93% of the variants.

## Results (desktop: 24 cores, 126 GB, 9 November 2026)

`results/desktop_24c.tsv` and `results/desktop_24c_summary.md`: 1.0.0
against vcf_fast (commit 87598491) at 2.3, 4.6 and 18.6 M variants, with
`THREADS=4 MEMORY_MAX=50G` and `OLD_STEPS` set to all steps but `write_vcf`
for 1.0.0 (quadratic in the number of chromosomes, over 2 min already at
200,000 variants, so it would only time out). The two versions give the same
results on all steps.

At 18.6 M variants (17.3 M biallelic SNPs), vcf_fast:
- reads the VCF in 36.5 s and 3.6 GB (1.0.0: 1.8 min and 5.6 GB; published
  0.4.4: 30 min; VCFArrow: 6.9 min), or 21.3 s with `names_as_int = TRUE`
  (1.0.0: 1.4 min). With 15 threads (`results/desktop_24c_read_t15.tsv`) it
  takes 35.4 s and 20.2 s: on a plain (not bgzipped) VCF, reading is limited
  by its serial part (finding the lines and the loci, looking up
  chromosomes, creating R strings), not by parsing the genotypes;
- writes the VCF in 47 s (published 0.4.4: over 60 min; VCFArrow: 8.5 min);
- converts to genlight in 1.1 min and 9.5 GB (1.0.0: 2.0 min and 19 GB);
- merges two sample sets with `rbind()` in 39 s (1.0.0: 53 s; VCFArrow:
  10.3 min);
- runs every other function of the comparison in under 20 s.

Where tidypopgen was still slow in this run:
- `subsample_loci_bare` (a bare index vector, without `all_of()`): 5.9 min
  and 8.9 GB in both versions, and quadratic in the number of loci (7 s, 26 s
  and 355 s at the three sizes), because `tidyselect` creates an R symbol for
  every locus name. The first run of vcf_fast at 4.6 M variants gave 32.9 s,
  but other work was running at the time; on a quiet machine both versions
  take 26 s (`results/desktop_24c_recheck_subsample_bare.tsv`). Fixed later
  that day on dev (`2826366e`), see below.
- `to_ped` (20 s, to write 3.9 GB of text) and `filter_chain_update` (17 s:
  the filters of the VCFArrow vignette, then `gt_update_backingfile()`,
  which writes a new copy of the data), profiled below.

## Rechecks on dev (desktop, 13 November 2026)

### `select_loci()` with a bare variable

`results/recheck_subsample_bare_dev.tsv`: dev at `216165b7` (installed by
`setup_libs.sh`, -O2), 4 threads, 3 runs per size (load average 1.4-2.1 on
24 cores). `select_loci(gt, idx)` is now resolved without `tidyselect`
(`2826366e`), and is linear in the number of loci:

| variants (loci) | `all_of(idx)` | bare `idx` (best of 3) | bare `idx` on 9 Nov |
|---|---|---|---|
| 2.3 M (2.16 M) | 0.004 s | 0.06 s | 7 s |
| 4.6 M (4.33 M) | 0.006 s | 0.09 s | 26 s |
| 18.6 M (17.3 M) | 0.003 s | 0.27 s | 355 s, 8.9 GB |

Peak PSS stays that of the loaded gen_tibble (2.96 GB at 18.6 M; at most
+0.2 GB). The bare form still costs more than `all_of()` because it scans
the loci names once, to check that no locus is named like the variable
(which `tidyselect` would select instead); a vector of names given by name
is resolved in the same way (0.09 s and 0.16 s at the two smaller sizes).

What still goes through `tidyselect` (and its data mask of every locus
name, never freed): a bare locus name (`select_loci(gt, rs1)`, or a
variable named like a locus), combinations (`c(all_of(idx), 1)`,
`starts_with("a") & !ends_with("b")`, binary `-`), `where()`, anything
`.engine = "tidyselect"`, and any selection the fast path cannot resolve
(e.g. missing or out-of-range positions, which `tidyselect` turns into its
errors). Measured on dev in a fresh process per call (`select_loci()` alone,
1 thread; a bare locus name, `c(all_of(idx), 1)` and `all_of(idx)` with
`.engine = "tidyselect"` all cost the same):

| loci | time (best of 3) | max RSS of the process | `all_of(idx)` |
|---|---|---|---|
| 2.16 M | 7.6-10.5 s | 1.26-1.31 GB | 0.00 s, 0.51 GB |
| 4.33 M | 27-28 s | 2.41-2.49 GB | 0.00 s, 0.89 GB |

So still superlinear (3.5 times longer for twice the loci, about 355 s at
17.3 M as measured on 9 November), and about 350 bytes per locus. Other
agents were compiling and testing on the machine during these runs (load
average 1.5-8 at 2.16 M, 9-18 at 4.33 M), but the 4.33 M times match the
26 s of the quiet run on 9 November.

### `filter_chain_update`

At 18.6 M variants (17.3 M loci, 13 M chromosomes), 4 threads, best of 3
runs, dev (`216165b7`) against the branch `select_filter_speed`, both
installed by `setup_libs.sh`. Load average 2-8 during the whole step
(`results/recheck_filter_chain.tsv`) and 9-12 during the parts (each part
after a full `gc()`):

| part | dev | branch |
|---|---|---|
| whole step (harness) | 18.4 s, 5.6 GB | 8.1 s, 5.2 GB |
| `filter(population != "OG")` | 0.00 s | 0.00 s |
| `select_loci_if(loci_maf() >= 0.03)` | 3.2 s | 3.1 s |
| - of which `loci_maf()` | 0.24 s | 0.23 s |
| `select_loci_if(loci_missingness() <= 0.2)` | 0.29 s | 0.30 s |
| `filter(indiv_missingness() <= 0.3)` | 0.13 s | 0.12 s |
| `gt_update_backingfile()` | 12.6 s | 4.0 s |
| - check of genetic distances | 8.3 s | skipped (3.1 s if run) |
| - `.gt` file (`gt_save_light()`) | 3.8 s | 3.7 s |

- `gt_update_backingfile()` checked that genetic distances are sorted within
  each chromosome, casting the 13 M chromosome names to integers (`gsub()`,
  `grepl()`, `unique()` and `match()` on 13 M strings), although they are
  all zero. The check is now skipped when the distances are sorted across
  the whole table, and `cast_chromosome_to_int()` takes 2.5 s rather than
  7.8 s (it is also used by `is_loci_table_ordered()`, `gt_order_loci()` and
  others). The copy of the genotypes takes about 0.3 s.
- What remains is the `.gt` file (saving the loci table, left unchanged on
  purpose) and garbage collection: the first `select_loci_if()` spends
  almost all its time in R's garbage collector, which has to scan the
  ~30 M live strings (loci and chromosome names) on each full collection;
  `loci_maf()` itself takes 0.24 s.
- MAF and missingness in a single pass (`qc_report_loci(x, hwe = FALSE)`,
  new on the branch) gives the same loci: the two kernels take 0.26 s +
  0.13 s against 0.25 s for the single pass on 4 threads (0.38 s + 0.22 s
  against 0.50 s on 1 thread), and the two `select_loci_if()` 1.6 s against
  0.65 s for `qc_report_loci()` and one `select_loci_if()` (the saving is
  mostly one garbage collection less). Load average 9 during these runs.

Light work (single-core, niced compiles and small tests) ran on the machine
during part of the run (vcf_fast at 4.6 M variants and the first steps of
1.0.0 at 18.6 M).

## Earlier results (laptop: 4 cores, 15 W, 15 GB)

- Smoke test, 20,000 variants, 2 threads: all steps run in 0.4.4 and 1.0.0
  with identical results; `gt_as_plink(type = "ped")` 21 s in 0.4.4, 0.09 s
  in 1.0.0; `gt_as_vcf()` 2.5 s in both, but quadratic in the number of
  chromosomes in 1.0.0 (fixed in vcf_fast: 186k loci, 141k chromosomes,
  271 s -> 1.3 s, output byte-identical).
- `results/laptop_test_200k.tsv` and `results/laptop_test_200k_summary.md`:
  200,000 variants, 4 threads, 2 min per step, 1.0.0 against the branch
  (then called gt_as_vcf_speed, commit 79605daf): identical results on all
  steps; write VCF > 2 min -> 1.1 s, read VCF 2.1 s -> 1.1 s.
- Reading the 2.3 M-variant VCF (`gen_tibble(parser = "cpp")`): 0.4.4 292 s;
  vcf_fast 11.1 s (1 core) and 10.1 s (4 cores), or 16.8 s and 12.4 s for
  the bgzipped file.
