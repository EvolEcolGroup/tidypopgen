#!/usr/bin/env bash
set -euo pipefail

# Requires: bcftools, wget (or curl).
# All paths are relative to the folder containing this script (data-raw/)

cd "$(dirname "${BASH_SOURCE[0]}")"

TAR=vcfs.variants.public_samples.279samples.tar
URL=https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/$TAR
POPS=(afr oceania america west_eur)   # expects ${pop}_subset.txt for each
EXTRACT_DIR=extracted

# 0. Download the SGDP tar if it is not already here (large file!)
[ -f "$TAR" ] || wget -O "$TAR" "$URL"

mkdir -p "$EXTRACT_DIR"

# 1. Extract only the files needed across all four ID lists (once)
cat "${POPS[@]/%/_subset.txt}" | sed '/^[[:space:]]*$/d' | sort -u > all_ids.txt
tar -tf "$TAR" | grep -F -f all_ids.txt > wanted.txt
tar -xf "$TAR" -C "$EXTRACT_DIR" -T wanted.txt

# 2. Index every extracted VCF (required for bcftools merge)
find "$EXTRACT_DIR" -name '*.vcf.gz' | while read -r f; do
    [ -f "$f.tbi" ] || [ -f "$f.csi" ] || bcftools index -t "$f"
done

# 3. Build one merged VCF per population
for pop in "${POPS[@]}"; do
    list="${pop}_files.txt"
    : > "$list"

    while read -r id; do
        [ -z "${id// }" ] && continue
        match=$(find "$EXTRACT_DIR" -name "*${id}*.vcf.gz")
        if [ -z "$match" ]; then
            echo "WARNING [$pop]: no VCF found for ID '$id'" >&2
        else
            echo "$match" >> "$list"
        fi
    done < "${pop}_subset.txt"

    sort -u -o "$list" "$list"
    echo "[$pop] merging $(wc -l < "$list") files"

    bcftools merge -l "$list" -Oz -o "${pop}_subset.vcf.gz"
    bcftools index -t "${pop}_subset.vcf.gz"
done

bcftools merge "${POPS[@]/%/_subset.vcf.gz}" -Oz -o all_subset.vcf.gz
bcftools index -t all_subset.vcf.gz

echo "Done."
