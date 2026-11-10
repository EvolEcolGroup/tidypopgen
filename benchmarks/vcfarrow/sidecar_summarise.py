#!/usr/bin/env python3
"""Markdown tables of the sidecar benchmark (sidecar_bench.sh):
  sidecar_summarise.py [results/sidecar_cost.tsv] [results/sidecar_compress.tsv]
                       [results/sidecar_chromosomes.tsv]
Times are the mean of the repeats, with their range; peak memory is the peak
PSS of the R process during the step (sampled every 0.2 s)."""
import csv
import json
import statistics
import sys
from collections import defaultdict

cost_file = sys.argv[1] if len(sys.argv) > 1 else "results/sidecar_cost.tsv"
comp_file = sys.argv[2] if len(sys.argv) > 2 else \
    "results/sidecar_compress.tsv"
chr_file = sys.argv[3] if len(sys.argv) > 3 else \
    "results/sidecar_chromosomes.tsv"


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def mean_range(xs, fmt="{:.2f}", stat=statistics.mean):
    xs = [x for x in xs if x == x]
    if not xs:
        return "NA"
    m = stat(xs)
    if len(xs) == 1:
        return fmt.format(m)
    return (fmt.format(m) + " (" + fmt.format(min(xs)) + "-" +
            fmt.format(max(xs)) + ")")


def mb(x):
    return f"{x / 1e6:,.0f}"


cost = read_tsv(cost_file)
rows = defaultdict(list)
for r in cost:
    if r["status"] != "ok":
        print(f"<!-- not ok: {r['label']} {r['n_variants']} {r['step']} "
              f"{r['status']} {r['summary']} -->")
        continue
    rows[(r["label"], int(r["n_variants"]), r["step"])].append(r)
sizes = sorted({k[1] for k in rows})


def col(label, n, step, field):
    return [num(r[field]) for r in rows.get((label, n, step), [])]


def extra(label, n, step):
    return [json.loads(r["extra"]) for r in rows.get((label, n, step), [])]


# 1. reading
print("### 1. Reading with names_as_int = TRUE (n_cores = 4)\n")
print("Wall time: median (range) of the repeats; the difference is between "
      "the medians (and between the means). Peak memory: the largest RSS of "
      "the R process (getrusage, so including gt_save() after the read), "
      "median of the repeats.\n")
print("| variants | biallelic SNPs | reps | baseline (s) | sidecar (s) | "
      "difference (s) | difference (%) | difference of means (%) | "
      "CPU time difference (s) | names_as_int = FALSE (s) | "
      "max RSS baseline (GB) | max RSS sidecar (GB) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
med = statistics.median
for n in sizes:
    b = col("base", n, "read_int", "elapsed_s")
    s = col("sidecar", n, "read_int", "elapsed_s")
    c = col("sidecar", n, "read_chr", "elapsed_s")
    if not b or not s:
        continue
    cpu_b = [u + y for u, y in zip(col("base", n, "read_int", "user_s"),
                                   col("base", n, "read_int", "sys_s"))]
    cpu_s = [u + y for u, y in zip(col("sidecar", n, "read_int", "user_s"),
                                   col("sidecar", n, "read_int", "sys_s"))]
    nl = extra("sidecar", n, "read_int")[0]["n_loci"]
    d = med(s) - med(b)
    dm = statistics.mean(s) - statistics.mean(b)
    print(f"| {n:,} | {nl:,} | {len(b)}/{len(s)} | "
          f"{mean_range(b, stat=med)} | {mean_range(s, stat=med)} | "
          f"{d:+.2f} | {100 * d / med(b):+.1f}% | "
          f"{100 * dm / statistics.mean(b):+.1f}% | "
          f"{med(cpu_s) - med(cpu_b):+.2f} | {mean_range(c, stat=med)} | "
          f"{med(col('base', n, 'read_int', 'proc_maxrss_gb')):.3f} | "
          f"{med(col('sidecar', n, 'read_int', 'proc_maxrss_gb')):.3f} |")
print()
print("| variants | .bk (MB) | .rds (MB) | .gt (MB) | sidecar (MB) | "
      "bytes per locus | sidecar / .bk |")
print("|---|---|---|---|---|---|---|")
for n in sizes:
    e = extra("sidecar", n, "read_int")
    g = extra("sidecar", n, "gt_save")
    if not e:
        continue
    e = e[0]
    gt = g[0]["gt_bytes"] if g else float("nan")
    print(f"| {n:,} | {mb(e['bk_bytes'])} | {mb(e['rds_bytes'])} | "
          f"{mb(gt)} | {mb(e['sidecar_bytes'])} | "
          f"{e['sidecar_bytes'] / e['n_loci']:.1f} | "
          f"{100 * e['sidecar_bytes'] / e['bk_bytes']:.0f}% |")
print()

# 3. lookups, restore, update
print("### 3. Lookups, restore and update (sidecar version)\n")
print("| variants | step | time (s) | notes |")
print("|---|---|---|---|")
for n in sizes:
    if not rows.get(("sidecar", n, "lookup100_first")):
        continue
    first = rows[("sidecar", n, "lookup100_first")]
    warm_cache = [num(r["elapsed_s"]) for r in first
                  if r["action"] == "lookup"]
    cold = [num(r["elapsed_s"]) for r in first
            if r["action"] == "lookup_cold"]
    cached = [num(r["cached_frac"]) for r in first
              if r["action"] == "lookup_cold"]
    print(f"| {n:,} | first lookup of 100 rows (builds the index), sidecar "
          f"in page cache | {mean_range(warm_cache, '{:.3f}')} | |")
    print(f"| {n:,} | first lookup of 100 rows, sidecar evicted from the "
          f"page cache | {mean_range(cold, '{:.3f}')} | "
          f"{100 * max(cached):.0f}% of the file cached after eviction |")
    idx = extra("sidecar", n, "sidecar_index")
    if idx:
        print(f"| {n:,} | sidecar_index() alone (median of 3) | "
              f"{idx[0]['median']:.3f} | "
              f"{idx[0]['sidecar_bytes'] / 1e6:,.0f} MB, "
              f"{idx[0]['sidecar_bytes'] / 1e9 / idx[0]['median']:.1f} "
              f"GB/s |")
    same = col("sidecar", n, "lookup100_same", "elapsed_s")
    print(f"| {n:,} | the same 100 rows again | "
          f"{mean_range(same, '{:.4f}')} | |")
    for size in [100, 1000, 10000, 100000, 1000000]:
        for kind, label in [("lookup", "gt_lookup_loci()"),
                            ("fetch", "C++ fetch alone")]:
            e = extra("sidecar", n, f"{kind}{size}_warm")
            if not e:
                continue
            t = e[0]["times"]
            print(f"| {n:,} | warm {label}, {size:,} random rows | "
                  f"{statistics.median(t):.4f} | median of {len(t)} "
                  f"(range {min(t):.4f}-{max(t):.4f}); "
                  f"{1e6 * statistics.median(t) / size:.1f} us per row |")
print()
print("| variants | step | time (s) | peak PSS (GB) | PSS before (GB) | "
      "notes |")
print("|---|---|---|---|---|---|")
for n in sizes:
    for label, step, note in [
            ("sidecar", "restore", "gt_restore_loci(), all loci"),
            ("base", "select10", "select_loci() of 10% (baseline)"),
            ("sidecar", "select10", "select_loci() of 10% (sidecar)"),
            ("base", "update10", "gt_update_backingfile(), baseline"),
            ("sidecar", "update10", "gt_update_backingfile(), sidecar"),
            ("sidecar", "sidecar_copy10",
             "copy_loci_sidecar() alone (part of the above)")]:
        t = col(label, n, step, "elapsed_s")
        if not t:
            continue
        e = extra(label, n, step)
        sz = ""
        if e and e[0].get("sidecar_bytes") is not None:
            sz = f"; new sidecar {e[0]['sidecar_bytes'] / 1e6:,.1f} MB"
        if e and e[0].get("bk_bytes") is not None:
            sz += f", new .bk {e[0]['bk_bytes'] / 1e6:,.1f} MB"
        print(f"| {n:,} | {note} | {mean_range(t, '{:.3f}')} | "
              f"{mean_range(col(label, n, step, 'pss_peak_gb'))} | "
              f"{mean_range(col(label, n, step, 'pss_start_gb'))} | "
              f"{len(t)} reps{sz} |")
print()

# the parts of gt_restore_loci() (results/sidecar_restore_parts.tsv)
try:
    parts = read_tsv(cost_file.replace("sidecar_cost", "sidecar_restore_parts"))
except OSError:
    parts = []
if parts:
    print("Parts of gt_restore_loci() (sidecar_step.R actions restore, "
          "restore_unchecked and restore_parts; each action in a fresh R "
          "process):\n")
    print("| variants | action | step | time (s) | peak PSS (GB) |")
    print("|---|---|---|---|---|")
    for r in parts:
        if r["step"] in ("gt_load", "gt_save", "read_int"):
            continue
        print(f"| {int(r['n_variants']):,} | {r['action']} | {r['step']} | "
              f"{num(r['elapsed_s']):.3f} | {num(r['pss_peak_gb']):.2f} |")
    print()

# 2. compression
print("### 2. Compressing the sidecar\n")
comp = read_tsv(comp_file)
print("| variants | codec | size (MB) | ratio | compress 1 thread (s) | "
      "compress 4 threads (s) | decompress 1 thread (s) | "
      "decompress 4 threads (s) | decompress + scan lines, 1 thread (s) | "
      "random lookup (us) |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in comp:
    if r["codec"] == "plain":
        print(f"| {int(r['n_variants']):,} | plain TSV | "
              f"{mb(num(r['out_bytes']))} | 1 | | | | | "
              f"{num(r['decomp_scan_1t_s']):.3f} (scan only) | "
              f"{num(r['lookup_us']):.0f} (pread 1 MB), "
              f"{num(r['lookup64k_us']):.0f} (pread 64 KB) |")
        continue

    def f(x, d=2):
        v = num(r[x])
        return "" if v != v else f"{v:.{d}f}"
    print(f"| {int(r['n_variants']):,} | {r['codec']} | "
          f"{mb(num(r['out_bytes']))} | {f('ratio')} | {f('comp_1t_s')} | "
          f"{f('comp_nt_s')} | {f('decomp_1t_s')} | {f('decomp_nt_s')} | "
          f"{f('decomp_scan_1t_s')} | {f('lookup_us', 0)} |")
print()

# 4. chromosome sidecar
print("### 4. Chromosome sidecar\n")
print("| variants | biallelic SNPs | chromosomes | locus sidecar (MB) | "
      "chromosome sidecar (MB) | bytes per chromosome |")
print("|---|---|---|---|---|---|")
for r in read_tsv(chr_file):
    k = int(r["n_chromosomes"])
    b = int(r["chromosome_sidecar_bytes"])
    print(f"| {int(r['n_variants']):,} | {int(r['n_loci']):,} | {k:,} | "
          f"{mb(int(r['loci_sidecar_bytes']))} | {mb(b)} | {b / k:.1f} |")
