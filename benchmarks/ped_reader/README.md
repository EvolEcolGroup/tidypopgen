# Benchmark: gen_tibble() on PLINK .ped files

This harness times `gen_tibble()` on a synthetic *.ped* file, for two
installed versions of tidypopgen: before (`dev` at `216165b7`, where
`read_pedfile()` read the lines in R and held the genotype matrix in memory)
and after the streaming reader (`ped_to_fbm()`, which writes the genotypes
straight into the backing file).

`bench_step.R` also saves the loci and individual tables and the md5 sum of
the backing file of each run: all runs below gave identical tables and
backing files.

## Running it

```sh
# install each version into its own library (never time a load_all() build,
# which is compiled at -O0)
git archive <ref> | tar -x -C <src_dir>
R CMD INSTALL --library=<lib> <src_dir>
# 2,000 individuals x 50,000 loci, 2% missing (a 382 MB .ped file)
Rscript make_data.R <lib> <data_dir> 2000 50000
# check that the machine is idle (uptime, top), then
./run_bench.sh <data_dir>/synthetic.ped <work_dir> 3 \
  old=<old_lib>:1 new=<new_lib>:1 new4=<new_lib>:4
```

## Results

November 2026, on a 24-core i9-13900 (shared machine, load average 2.5 to
4.5 from other jobs during the runs), best of 3 runs (1 run for the old
reader on the larger file). Peak RSS includes about 135 MB for R with
tidypopgen loaded, and the pages of the backing file mapped while it is
written (95 MB and 381 MB); the peak R heap is from `gc()`.

| file | reader | n_cores | time (s) | peak RSS (MB) | peak R heap (MB) |
|---|---|---|---|---|---|
| 2,000 x 50,000 (382 MB) | old | 1 | 44.7 | 943 | 804 |
| | new | 1 | 1.03 | 279 | 132 |
| | new | 4 | 0.49 | 279 | 132 |
| 8,000 x 50,000 (1.5 GB) | old | 1 | 184.6 | 2806 | 2334 |
| | new | 1 | 4.25 | 566 | 133 |
| | new | 4 | 1.57 | 566 | 133 |

The R heap of the new reader does not grow with the data (it holds the
loci and individual tables); its RSS only grows with the mapped pages of the
backing file, which the kernel can write back and drop.
