#!/usr/bin/env bash
# Set up the gt_merge_vcfs() timings on the desktop (manica01): the three
# synthetic datasets, the four builds, bcftools 1.25 and plink2. Untimed: run
# it on the E-cores while the machine is in use, before run_desktop.sh.
#   bash dev_notes/vcf_merge/setup_desktop.sh
# Needs: R with data.table, bgzip and tabix (htslib), a C toolchain, internet.
set -euo pipefail
REPO=${TPG_REPO:-$HOME/git/tidypopgen}
export TPG_WORK=${TPG_WORK:-$HOME/claude_work/tidypopgen}
W=$TPG_WORK/merge_vcf
CPUS=${CPUS:-8-23}  # E-cores
mkdir -p "$W"/{data,libs,tools/src,results/meta,rtmp,out,scripts}
export TMPDIR=$W/rtmp
cp "$REPO"/dev_notes/vcf_merge/scripts/* "$W/scripts/"
chmod +x "$W"/scripts/*.sh

# 1. bcftools 1.25 from the release tarball (its bundled htslib, -O2)
if [ ! -x "$W/tools/bcftools/bin/bcftools" ]; then
  cd "$W/tools/src"
  curl -sSLO https://github.com/samtools/bcftools/releases/download/1.25/bcftools-1.25.tar.bz2
  tar -xjf bcftools-1.25.tar.bz2 && cd bcftools-1.25
  ./configure --prefix="$W/tools/bcftools" CFLAGS="-g -O2" --disable-libcurl > ../configure.log
  taskset -c "$CPUS" make -j8 > ../make.log && make install >> ../make.log
fi
B=$W/tools/bcftools/bin/bcftools
# bgzip and tabix: from the htslib 1.25 bundled with bcftools when they are not
# installed (manica01 had neither in Nov 2026; `sudo apt install tabix` would
# give htslib 1.19)
if ! command -v bgzip > /dev/null || ! command -v tabix > /dev/null; then
  if [ ! -x "$W/tools/bcftools/bin/bgzip" ]; then
    taskset -c "$CPUS" make -C "$W/tools/src/bcftools-1.25/htslib-1.25" -j8 bgzip tabix \
      >> "$W/tools/src/make.log"
    cp "$W"/tools/src/bcftools-1.25/htslib-1.25/{bgzip,tabix} "$W/tools/bcftools/bin/"
  fi
  export PATH=$W/tools/bcftools/bin:$PATH
fi

# 2. plink2: the AVX2 Linux build (the laptop used v2.0.0-b.1.4LM of 7 Nov 2026)
if [ ! -x "$TPG_WORK/plink_king/bin/plink2" ]; then
  echo "Download the AVX2 Linux build of plink2 from https://www.cog-genomics.org/plink/2.0/"
  echo "and put the binary at $TPG_WORK/plink_king/bin/plink2 (tools_timing.sh looks there)."
fi

# 3. the datasets (seeded: identical to the laptop's), bgzipped and indexed
cd "$W"
for sc in pops sgdp many; do
  [ -d "data/$sc" ] && ls data/$sc/*.vcf.gz > /dev/null 2>&1 && continue
  taskset -c "$CPUS" Rscript scripts/make_merge_data.R data "$sc"
  for f in data/$sc/*.vcf; do taskset -c "$CPUS" bgzip -f "$f"; done
done
# the sgdp headers declare their INFO fields (bcftools needs them), as on the laptop
INFO=$(for x in "DP Integer" "MQ Float" "FS Float" "QD Float" "CpG Integer" \
  "RepeatMasker String" "Map20 Integer" "TandemRepeat Integer" "SegDup Integer" \
  "ANN String"; do set -- $x
  printf '##INFO=<ID=%s,Number=1,Type=%s,Description="annotation">\\n' "$1" "$2"; done)
mkdir -p rtmp/hdr
# zcat | head ends zcat with SIGPIPE: without this, pipefail makes the check
# below always false and set -e stops the script
set +o pipefail
for f in data/sgdp/*.vcf.gz; do
  zcat "$f" | head -40 | grep -q '^##INFO' && continue
  b=$(basename "$f" .vcf.gz)
  zcat "$f" 2>/dev/null | head -40 | grep "^#" |
    awk -v info="$INFO" '/^##contig/ && !done {printf "%s", info; done=1} {print}' > rtmp/hdr/$b.txt
  "$B" reheader -h rtmp/hdr/$b.txt -o rtmp/hdr/$b.vcf.gz "$f" && mv rtmp/hdr/$b.vcf.gz "$f"
done
set -o pipefail
for f in data/*/*.vcf.gz; do [ -e "$f.tbi" ] || tabix -f -p vcf "$f"; done
du -sh data/*

# 4. the four versions, each installed at -O2 into its own library
declare -A SHAS=([base]=5e207d78 [p1]=d3b982d0 [p3]=bc6856b6 [p4]=6982994c)
for v in base p1 p3 p4; do
  [ -e "$W/libs/$v/COMMIT" ] && continue
  rm -rf "$W/src_copies/$v" && mkdir -p "$W/src_copies/$v" "$W/libs/$v"
  git -C "$REPO" archive "${SHAS[$v]}" | tar -x -C "$W/src_copies/$v"
  taskset -c "$CPUS" R CMD INSTALL --library="$W/libs/$v" "$W/src_copies/$v" > "$W/results/install_$v.log" 2>&1
  echo "${SHAS[$v]}" > "$W/libs/$v/COMMIT"
done
echo "setup done; check: each run's md5 in results/runs.tsv must be pops 87fe0252_1cb5b27f,"
echo "sgdp 2e2fdeae_e4f9c90e, many f286411f_d33f758d (as on the laptop)"
