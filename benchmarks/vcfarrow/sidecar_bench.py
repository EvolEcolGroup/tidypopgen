#!/usr/bin/env python3
"""Run one action of sidecar_step.R (the cost of the locus sidecar of
names_as_int = TRUE) in a fresh R process for one installed tidypopgen,
sampling the memory of its process tree as bench.py does (PSS and RSS summed
over the R process and its descendants), in a systemd user scope with
--memory-max and no swap. Appends one row per timed step to --results.

  sidecar_bench.py --lib lib/sidecar --label sidecar --vcf data/x.vcf \
      --work work/sidecar_sidecar_x --results results/sidecar_cost.tsv \
      --action read_int [--threads 4] [--rep 1] [--memory-max 50G]
      [--warm-vcf] [--evict FILE]

--warm-vcf reads the VCF once before the step (so every read starts with the
VCF in the page cache, as in bench.py); --evict drops FILE (e.g. the
sidecar) from the page cache before the step (posix_fadvise DONTNEED, after a
sync), and reports how much of it was still cached (mincore).
"""
import argparse
import ctypes
import ctypes.util
import json
import mmap
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# (bench.py only runs its main() as a script)
from bench import tree_memory, window_stats, read_new_records, fmt  # noqa

COLUMNS = ["label", "commit", "version", "n_variants", "threads", "action",
           "rep", "step", "status", "elapsed_s", "user_s", "sys_s",
           "pss_start_gb", "pss_median_gb", "pss_peak_gb", "rss_peak_gb",
           "proc_maxrss_gb", "load1", "busy_cores", "cached_frac", "extra",
           "summary"]

libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)


def cached_fraction(path):
    """fraction of the pages of path in the page cache (mmap + mincore)"""
    size = os.path.getsize(path)
    if size == 0:
        return float("nan")
    libc.mmap.restype = ctypes.c_void_p
    libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int,
                          ctypes.c_int, ctypes.c_int, ctypes.c_long]
    libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t,
                             ctypes.POINTER(ctypes.c_ubyte)]
    fd = os.open(path, os.O_RDONLY)
    try:
        addr = libc.mmap(None, size, mmap.PROT_READ, mmap.MAP_SHARED, fd, 0)
        if addr in (None, ctypes.c_void_p(-1).value):
            return float("nan")
        n_pages = (size + mmap.PAGESIZE - 1) // mmap.PAGESIZE
        vec = (ctypes.c_ubyte * n_pages)()
        ok = libc.mincore(addr, size, vec) == 0
        libc.munmap(addr, size)
        if not ok:
            return float("nan")
        return sum(v & 1 for v in vec) / n_pages
    finally:
        os.close(fd)


def busy_cores(seconds=2.0):
    """the number of cores busy over the next `seconds` (from /proc/stat)"""
    def read():
        with open("/proc/stat") as f:
            v = [int(x) for x in f.readline().split()[1:9]]
        return sum(v), v[3] + v[4]  # all, idle + iowait
    t1, i1 = read()
    time.sleep(seconds)
    t2, i2 = read()
    return (1 - (i2 - i1) / max(t2 - t1, 1)) * os.cpu_count()


def wait_idle(max_busy, max_wait):
    """wait (up to max_wait s) until at most max_busy cores are busy (with
    other work: this runs between our measurements); returns the number of
    busy cores"""
    t0 = time.time()
    while True:
        b = busy_cores()
        if b <= max_busy or time.time() - t0 > max_wait:
            return b
        top = subprocess.run(["ps", "-eo", "pcpu,comm", "--sort=-pcpu"],
                             capture_output=True, text=True).stdout
        top = "; ".join(" ".join(x.split()) for x in top.splitlines()[1:4])
        print(f"    waiting: {b:.1f} cores busy ({top})", flush=True)
        time.sleep(15)


def evict(path):
    subprocess.run(["sync"])
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--action", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--rep", type=int, default=1)
    ap.add_argument("--interval", type=float, default=0.2)
    ap.add_argument("--memory-max", default="50G")
    ap.add_argument("--warm-vcf", action="store_true")
    ap.add_argument("--evict", default=None)
    ap.add_argument("--max-busy", type=float, default=3.0,
                    help="wait until at most this many cores are busy")
    ap.add_argument("--max-wait", type=float, default=900)
    args = ap.parse_args()
    args.lib = os.path.abspath(args.lib)
    args.vcf = os.path.abspath(args.vcf)
    args.work = os.path.abspath(args.work)
    os.makedirs(args.work, exist_ok=True)

    with open(os.path.join(args.lib, "COMMIT")) as f:
        commit = f.read().strip()[:8]
    nv_file = args.vcf.rsplit(".vcf", 1)[0] + ".n_variants"
    with open(nv_file) as f:
        n_variants = int(f.read().strip())

    # same starting point for every run: no dirty pages being written back
    # (from the outputs of the previous run), the VCF in the page cache
    subprocess.run(["sync"])
    if args.warm_vcf:
        subprocess.run(["cat", args.vcf], stdout=subprocess.DEVNULL)
    cached = float("nan")
    if args.evict:
        evict(args.evict)
        cached = cached_fraction(args.evict)
    busy = wait_idle(args.max_busy, args.max_wait)
    load1 = os.getloadavg()[0]

    results_jsonl = os.path.join(args.work, "steps.jsonl")
    n_before = 0
    if os.path.exists(results_jsonl):
        with open(results_jsonl) as f:
            n_before = sum(1 for _ in f)
    try:
        os.remove(os.path.join(args.work, "current_step"))
    except OSError:
        pass
    cmd = ["Rscript", os.path.join(HERE, "sidecar_step.R"), args.lib,
           args.work, args.vcf, str(args.threads), results_jsonl, args.action]
    if args.memory_max:
        cmd = ["systemd-run", "--user", "--scope", "--quiet",
               "-p", f"MemoryMax={args.memory_max}",
               "-p", "MemorySwapMax=0"] + cmd
    env = dict(os.environ)
    env["OPENBLAS_NUM_THREADS"] = "1"
    print(f"[{time.strftime('%H:%M:%S')}] {args.label} {n_variants} "
          f"{args.action} rep {args.rep} (load {load1:.2f}, "
          f"{busy:.1f} cores busy)", flush=True)
    log = open(os.path.join(args.work, f"log_{args.action}.txt"), "w")
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                            env=env, start_new_session=True)
    samples = []
    while True:
        pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
        if pid != 0:
            maxrss = ru.ru_maxrss
            break
        pss, rss, n = tree_memory(proc.pid)
        samples.append((time.time(), pss, rss, n))
        time.sleep(args.interval)
    log.close()
    records = read_new_records(results_jsonl, n_before)
    if not records or os.WEXITSTATUS(status) != 0:
        print(f"    R exited with status {status}: see "
              f"{args.work}/log_{args.action}.txt", flush=True)
    with open(os.path.join(args.work,
                           f"mem_{args.action}_{args.rep}.tsv"), "w") as f:
        f.write("time\tpss_kb\trss_kb\tn_proc\n")
        for s in samples:
            f.write("\t".join(str(v) for v in s) + "\n")

    new_file = not os.path.exists(args.results)
    with open(args.results, "a") as out:
        if new_file:
            out.write("\t".join(COLUMNS) + "\n")
        for r in records:
            start, med, peak, rss, _ = window_stats(
                samples, r["t_start"], r["t_end"])
            row = {
                "label": args.label, "commit": commit,
                "version": r.get("version"), "n_variants": n_variants,
                "threads": args.threads, "action": args.action,
                "rep": args.rep, "step": r["step"], "status": r["status"],
                "elapsed_s": r["elapsed"], "user_s": r["user"],
                "sys_s": r["sys"], "pss_start_gb": start,
                "pss_median_gb": med, "pss_peak_gb": peak,
                "rss_peak_gb": rss, "proc_maxrss_gb": maxrss / 1024 ** 2,
                "load1": load1, "busy_cores": busy, "cached_frac": cached,
                "extra": json.dumps(r.get("extra") or {},
                                    separators=(",", ":")),
                "summary": r["summary"]}
            out.write("\t".join(fmt(row[c]) for c in COLUMNS) + "\n")
            print(f"    {r['step']}: {r['status']}, {fmt(r['elapsed'])} s, "
                  f"peak PSS {fmt(peak)} GB; {r['summary']}", flush=True)
    return 0 if records else 1


if __name__ == "__main__":
    sys.exit(main())
