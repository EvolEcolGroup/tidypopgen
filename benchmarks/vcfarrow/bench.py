#!/usr/bin/env python3
"""Run the steps of steps.R for one tidypopgen library, one fresh R process
per step, sampling the memory (PSS and RSS, summed over the R process and all
its descendants, as in the VCFArrow comparison) every --interval seconds, and
killing a step that runs longer than --timeout seconds (60 min in the
comparison). With --memory-max (e.g. 12G; 50 GB in the comparison), each R
process runs in a systemd user scope with that MemoryMax and no swap, so a
step that needs more is killed and reported as out of memory.

  bench.py --lib lib/v1.0.0 --label v1.0.0 --vcf data/x.vcf \
           --work work/v1.0.0_x --results results/x.tsv [--threads 4]
           [--steps read_vcf,filter_maf,...] [--timeout 3600]
           [--memory-max 12G]

Writes one row per step to --results (appending), the raw step records to
<work>/steps.jsonl and the memory samples to <work>/mem_<step>.tsv.
"""
import argparse
import json
import re
import os
import signal
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# the 17 functions of the comparison, in the order of its tables
CORE = ["read_vcf", "extract_samples", "extract_groups", "subsample_loci",
        "filter_maf", "filter_locus_missing", "filter_sample_missing",
        "filter_invariant", "merge", "sample_stats", "theta_pi",
        "missing_plot", "write_vcf", "to_genlight", "to_ped", "to_bed"]
EXTRA = ["subsample_loci_bare", "theta_pi_grouped", "filter_chain_update",
         "read_vcf_int"]

COLUMNS = ["label", "version", "n_variants", "threads", "step", "status",
           "elapsed_s", "user_s", "sys_s", "load_s", "pss_start_gb",
           "pss_median_gb", "pss_peak_gb", "rss_peak_gb", "proc_maxrss_gb",
           "n_samples", "summary"]


def descendants(root):
    """pids of root and all its descendants"""
    children = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                stat = f.read()
            ppid = int(stat[stat.rindex(")") + 2:].split()[1])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(ppid, []).append(int(d))
    out, todo = [], [root]
    while todo:
        p = todo.pop()
        out.append(p)
        todo.extend(children.get(p, []))
    return out


def tree_memory(root):
    """(pss_kb, rss_kb, n_processes) summed over the process tree"""
    pss = rss = n = 0
    for p in descendants(root):
        try:
            with open(f"/proc/{p}/smaps_rollup") as f:
                for line in f:
                    if line.startswith("Pss:"):
                        pss += int(line.split()[1])
                    elif line.startswith("Rss:"):
                        rss += int(line.split()[1])
            n += 1
        except OSError:
            pass
    return pss, rss, n


def read_current_step(work):
    try:
        with open(os.path.join(work, "current_step")) as f:
            step, t = f.read().strip().split("\t")
        return step, float(t)
    except (OSError, ValueError):
        return None, None


def run_process(args, steps, results_jsonl):
    """run one R process for `steps`; returns (records, samples, maxrss_kb,
    timed_out_step)"""
    try:
        os.remove(os.path.join(args.work, "current_step"))
    except OSError:
        pass
    n_before = 0
    if os.path.exists(results_jsonl):
        with open(results_jsonl) as f:
            n_before = sum(1 for _ in f)
    cmd = ["Rscript", os.path.join(HERE, "run_step.R"), args.lib, args.work,
           args.vcf, str(args.threads), results_jsonl] + steps
    if args.memory_max:
        cmd = ["systemd-run", "--user", "--scope", "--quiet",
               "-p", f"MemoryMax={args.memory_max}",
               "-p", "MemorySwapMax=0"] + cmd
    log = open(os.path.join(args.work, f"log_{steps[0]}.txt"), "w")
    env = dict(os.environ)
    # no multi-threaded BLAS on top of the threads given to tidypopgen
    env["OPENBLAS_NUM_THREADS"] = "1"
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                            env=env, start_new_session=True)
    samples = []
    maxrss = None
    timed_out = None
    killed = False
    while True:
        pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
        if pid != 0:
            maxrss = ru.ru_maxrss
            proc.returncode = status
            killed = ((os.WIFSIGNALED(status) and os.WTERMSIG(status) == 9)
                      or (os.WIFEXITED(status)
                          and os.WEXITSTATUS(status) == 137))
            break
        pss, rss, n = tree_memory(proc.pid)
        samples.append((time.time(), pss, rss, n))
        cur, t_start = read_current_step(args.work)
        if cur is not None and time.time() - t_start > args.timeout:
            timed_out = (cur, t_start)
            os.killpg(proc.pid, signal.SIGKILL)
        time.sleep(args.interval)
    log.close()
    # (a step that died before recording anything may leave no results file)
    try:
        with open(results_jsonl) as f:
            records = [json.loads(line) for line in list(f)[n_before:]]
    except FileNotFoundError:
        records = []
    oom = killed and timed_out is None
    return records, samples, maxrss, timed_out, oom


def window_stats(samples, t0, t1):
    """PSS at start, median and peak of PSS and peak RSS (GB) in [t0, t1]"""
    gb = 1024 ** 2
    # drop samples taken after the process had gone
    samples = [s for s in samples if s[3] > 0]
    before = [s for s in samples if s[0] <= t0]
    inside = [s for s in samples if t0 <= s[0] <= t1]
    if not inside:
        # step shorter than the sampling interval: use the next sample
        after = [s for s in samples if s[0] > t1]
        inside = after[:1] or before[-1:]
    start = before[-1][1] / gb if before else float("nan")
    if not inside:
        return start, float("nan"), float("nan"), float("nan"), 0
    pss = [s[1] / gb for s in inside]
    return (start, statistics.median(pss), max(pss),
            max(s[2] for s in inside) / gb, len(inside))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--steps", default="all",
                    help="'all', 'core', or comma-separated step names")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--interval", type=float, default=0.5)
    ap.add_argument("--memory-max", default=None,
                    help="systemd MemoryMax for each R process, e.g. 12G")
    args = ap.parse_args()
    args.lib = os.path.abspath(args.lib)
    args.vcf = os.path.abspath(args.vcf)
    args.work = os.path.abspath(args.work)
    os.makedirs(args.work, exist_ok=True)

    if args.steps == "all":
        steps = CORE + EXTRA
    elif args.steps == "core":
        steps = CORE
    else:
        steps = args.steps.split(",")
    # read_vcf (with set_groups) first, if asked for or if there is no data
    if "read_vcf" not in steps and not os.path.exists(
            os.path.join(args.work, "gt.gt")):
        steps = ["read_vcf"] + steps

    # number of variants: from the sidecar written by make_vcf.R
    # (counting the lines of a 27 GB VCF takes minutes)
    sidecar = re.sub(r"\.vcf(\.gz)?$", ".n_variants", args.vcf)
    if os.path.exists(sidecar):
        with open(sidecar) as f:
            n_variants = int(f.read().strip())
    else:
        with open(args.vcf) as f:
            n_variants = sum(1 for line in f if not line.startswith("#"))
    with open("/proc/meminfo") as f:
        mem_total = int(f.readline().split()[1]) * 1024
    results_jsonl = os.path.join(args.work, "steps.jsonl")
    new_file = not os.path.exists(args.results)
    out = open(args.results, "a")
    if new_file:
        out.write("\t".join(COLUMNS) + "\n")

    for step in steps:
        if step == "read_vcf" and os.path.getsize(args.vcf) < 0.4 * mem_total:
            # same starting point for every version: the VCF in the page
            # cache, if it fits (at full size, 27 GB, it does not)
            subprocess.run(["cat", args.vcf], stdout=subprocess.DEVNULL)
        print(f"[{time.strftime('%H:%M:%S')}] {args.label}: {step}",
              flush=True)
        records, samples, maxrss, timed_out, oom = run_process(
            args, [step], results_jsonl)
        with open(os.path.join(args.work, f"mem_{step}.tsv"), "w") as f:
            f.write("time\tpss_kb\trss_kb\tn_proc\n")
            for s in samples:
                f.write("\t".join(str(v) for v in s) + "\n")
        load_s = next((r["elapsed"] for r in records
                       if r["step"] == "gt_load"), None)
        done = {r["step"] for r in records}
        if timed_out is not None and timed_out[0] not in done:
            records.append({
                "step": timed_out[0], "phase": "step",
                "t_start": timed_out[1], "t_end": time.time(),
                "elapsed": None, "user": None, "sys": None,
                "status": f"timeout (> {args.timeout / 60:.0f} min)",
                "summary": None, "version": None})
        elif step not in done:
            cur, t_start = read_current_step(args.work)
            started = cur == step
            if oom:
                status = f"out of memory (> {args.memory_max})"
            else:
                status = "crashed" + ("" if started else " before the step")
            records.append({
                "step": step, "phase": "step",
                "t_start": t_start if started else None,
                "t_end": time.time() if started else None,
                "elapsed": None, "user": None, "sys": None,
                "status": status,
                "summary": f"see {args.work}/log_{step}.txt",
                "version": None})
        version = next((r["version"] for r in records if r.get("version")),
                       "")
        for r in records:
            if r["phase"] != "step":
                continue
            if r["t_start"] is not None:
                t1 = r["t_end"] if r["t_end"] is not None else time.time()
                start, med, peak, rss, n = window_stats(
                    samples, r["t_start"], t1)
            else:
                start = med = peak = rss = float("nan")
                n = 0
            row = {
                "label": args.label, "version": version,
                "n_variants": n_variants, "threads": args.threads,
                "step": r["step"], "status": r["status"],
                "elapsed_s": r["elapsed"], "user_s": r["user"],
                "sys_s": r["sys"], "load_s": load_s,
                "pss_start_gb": start, "pss_median_gb": med,
                "pss_peak_gb": peak, "rss_peak_gb": rss,
                "proc_maxrss_gb": (maxrss / 1024 ** 2) if maxrss else None,
                "n_samples": n, "summary": r["summary"]}
            out.write("\t".join(fmt(row[c]) for c in COLUMNS) + "\n")
            out.flush()
            print(f"    {r['step']}: {r['status']}, {fmt(r['elapsed'])} s, "
                  f"peak PSS {fmt(peak)} GB; {r['summary']}", flush=True)
        if step == "read_vcf" and not os.path.exists(
                os.path.join(args.work, "gt.gt")):
            print("    no gen_tibble was saved: stopping", flush=True)
            break
    out.close()


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return "NA" if v != v else f"{v:.4g}"
    return str(v).replace("\t", " ").replace("\n", " ")


if __name__ == "__main__":
    sys.exit(main())
