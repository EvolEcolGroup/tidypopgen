#!/usr/bin/env bash
set -euo pipefail

# Usage: create_sgdp_subset_s02.sh [out_dir]
#
# Download the VCFs of the SGDP individuals listed in the ${pop}_subset.txt
# files in data_files/sgdp/, next to this script (written by
# create_sgdp_subset_s01_metadata.R), into <out_dir>/sgdp_vcfs/ (one VCF per
# individual). <out_dir> defaults to that folder (data-raw/data_files/sgdp/).
#
# SGDP distributes its per-individual VCFs as a single 57 GB tar archive. We
# do not download the archive: SGDP also publishes its listing (`tar -tv`
# output, in archive order), and since every file in a tar is a 512-byte
# header followed by the file contents padded to a multiple of 512 bytes, the
# listing tells us where each file starts. We then download just the files we
# need with HTTP range requests (about 10.8 GB for 50 individuals), into
# <out_dir>/sgdp_vcfs/. An interrupted download resumes when the script is run
# again.
#
# Each per-individual VCF only contains the sites where that individual
# carries a non-reference allele. The VCFs are then merged in R with
# tidypopgen::gt_merge_vcfs(variants = "union", missing_to_ref = TRUE), so that
# a site that is absent from an individual's VCF is set to homozygous
# reference (as with `bcftools merge --missing-to-ref`), rather than to
# missing (see vignettes/articles/big_data.Rmd.orig).
#
# Requires: curl, awk, and md5sum (or md5 on macOS).

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBSET_DIR=$SCRIPT_DIR/data_files/sgdp # the ${pop}_subset.txt files
OUT_DIR=${1:-$SUBSET_DIR}
mkdir -p "$OUT_DIR"
OUT_DIR=$(cd "$OUT_DIR" && pwd)

BASE_URL=https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants
TAR_URL=$BASE_URL/vcfs.variants.public_samples.279samples.tar
LIST_URL=$BASE_URL/vcfs.variants.public_samples.279samples.file.list
POPS=(afr oceania america west_eur east_asia) # expects ${pop}_subset.txt for each
DL_DIR=$OUT_DIR/sgdp_vcfs           # per-individual files go here

mkdir -p "$DL_DIR"

# 1. Get the listing of the archive and work out where each file starts:
#    "<offset of the header> <size> <name>"
[ -s "$DL_DIR/file.list" ] ||
  curl -fsSL --retry 5 -o "$DL_DIR/file.list" "$LIST_URL"
awk '{
  printf "%.0f %.0f %s\n", offset, $3, $6
  offset += 512 + int(($3 + 511) / 512) * 512
}' "$DL_DIR/file.list" > "$DL_DIR/offsets.txt"

file_size() { echo $(($(wc -c < "$1"))); }

md5_of() {
  if command -v md5sum > /dev/null; then
    md5sum "$1" | cut -d ' ' -f 1
  else
    md5 -q "$1"
  fi
}

# Download one file out of the archive into $DL_DIR
fetch() {
  local name=$1 dest=$DL_DIR/$1 line header size start part have found tries=0
  line=$(awk -v n="$name" '$3 == n' "$DL_DIR/offsets.txt")
  if [ -z "$line" ]; then
    echo "ERROR: $name is not in the SGDP archive" >&2
    return 1
  fi
  read -r header size _ <<< "$line"
  if [ -f "$dest" ] && [ "$(file_size "$dest")" -eq "$size" ]; then
    return 0
  fi

  # Check that the tar header at this offset names the file we expect, in case
  # the listing and the archive ever disagree. --max-filesize makes curl stop
  # straight away if the server ignores the range and sends the whole archive.
  found=$(curl -fsS --retry 5 --max-filesize 512 \
    -r "$header-$((header + 99))" "$TAR_URL" | tr -d '\0')
  if [ "$found" != "$name" ]; then
    echo "ERROR: expected $name at byte $header of the archive, found '$found'" >&2
    return 1
  fi

  start=$((header + 512))
  part=$dest.part
  touch "$part"
  while have=$(file_size "$part") && [ "$have" -lt "$size" ]; do
    tries=$((tries + 1))
    if [ "$tries" -gt 10 ]; then
      echo "ERROR: could not download $name" >&2
      return 1
    fi
    curl -fsS -r "$((start + have))-$((start + size - 1))" "$TAR_URL" >> "$part" ||
      sleep 5
  done
  if [ "$have" -ne "$size" ]; then
    echo "ERROR: $name is $have bytes, expected $size; deleted it" >&2
    rm -f "$part"
    return 1
  fi
  mv "$part" "$dest"
}

# Compare a downloaded file with the checksum SGDP ships alongside it
check_md5() {
  local expected actual
  expected=$(cut -d ' ' -f 1 "$DL_DIR/$1.md5sum")
  actual=$(md5_of "$DL_DIR/$1")
  if [ "$expected" != "$actual" ]; then
    echo "ERROR: checksum mismatch for $1; deleted it" >&2
    rm -f "$DL_DIR/$1"
    return 1
  fi
}

# 2. Download each individual's VCF, and check it (no index is needed to
#    merge them with gt_merge_vcfs())
n_vcfs=0
for pop in "${POPS[@]}"; do
  ids=$(sed '/^[[:space:]]*$/d' "$SUBSET_DIR/${pop}_subset.txt" | LC_ALL=C sort -u)
  echo "[$pop] downloading $(echo "$ids" | grep -c . || true) individuals"
  for id in $ids; do
    vcf=$id.annotated.nh2.variants.vcf.gz
    for f in "$vcf" "$vcf.md5sum"; do
      fetch "$f"
    done
    check_md5 "$vcf"
    n_vcfs=$((n_vcfs + 1))
  done
done

echo "Done: $n_vcfs VCFs in $DL_DIR/"
echo "Merge them with tidypopgen::gt_merge_vcfs() (see the big data article)."
