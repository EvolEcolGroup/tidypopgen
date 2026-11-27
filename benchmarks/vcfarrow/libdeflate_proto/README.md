# libdeflate vs zlib for bgzipped VCFs (prototype, November 2026)

Set aside for now: adopting it would add a dependency and a configure step.
This directory keeps the prototype in case we come back to it.

- `bgzf_proto.cpp`: Rcpp functions that decompress a BGZF file and count its
  bytes and lines. `zlib_stream()` reads it with gzread, as a plain gzip
  stream. `bgzf_parallel(file, n_threads, use_libdeflate)` reads batches of
  BGZF blocks and inflates them in parallel (OpenMP) with zlib raw inflate
  or with libdeflate. `crc_cost()` times the CRC32 of the blocks with zlib
  or libdeflate.
- `Makevars_proto`: flags for sourceCpp (`-lz -ldeflate`, OpenMP, -O2). It
  needs the system libdeflate (Debian/Ubuntu: `libdeflate-dev`).
- `time_proto.R`: `Rscript time_proto.R <file.vcf.gz>` times all of them.

Results on `data/discosnp_2322460_gz.vcf.gz` (373 MB compressed, 3.4 GB of
text). The laptop has an i7-10510U with 4 cores and a 15 W limit, so all-core
speed is well below single-core turbo, and the 4-thread times vary from run
to run.

| decompression (3.4 GB)        | 1 thread | 2 threads | 4 threads   |
|-------------------------------|----------|-----------|-------------|
| zlib, gzread (sequential)     | 5.6 s    |           |             |
| BGZF blocks, zlib             | 5.3 s    | 2.9-3.5 s | 1.7-2.6 s   |
| BGZF blocks, libdeflate       | 2.3 s    | 1.3-1.6 s | 0.85-1.3 s  |
| CRC32 of the blocks, zlib     | 1.07 s   |           | 0.35 s      |
| CRC32 of the blocks, libdeflate | 0.30 s |           | 0.24 s      |

With the branch `gt_as_vcf_speed` (commit 09a84bee), `gen_tibble()` reads
this file in 16.8 s on 1 core and 12.4 s on 4. Using libdeflate would save
about 3 s on 1 core and about 1 s on 4.

How to adopt it:
- The CRAN package `libdeflate` (Tyler Morgan-Wall, 1.25-0, MIT) always
  builds a static `libdeflate.a` with CMake. It installs it in
  `lib/<Sys.info()[["machine"]]>/`, with the header in `include/`.
- A downstream package locates it in a configure step with
  `Rscript -e 'cat(system.file("lib", Sys.info()[["machine"]], package = "libdeflate"))'`,
  then adds `-I<that>/../include` to CPPFLAGS and `-L<that> -ldeflate` to
  PKG_LIBS. Check the Windows paths before relying on this.
- Other options are the system library with a bundled fallback (as GMMAT
  and MAGEE do) or bundling it (as pgenlibr does).
- Only BGZF blocks can use it: libdeflate works on whole buffers. Plain
  gzip files would still be read with zlib.
- In tidypopgen, the place to change is `VcfLineReader::read_bgzf()` in
  `src/vcf_parser_fast.cpp`, which has a comment about this.
