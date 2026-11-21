Desktop (manica01), 2026-11-21; medians (min-max) in seconds.

## KING-robust on HGDP650 (PR #239)

| task | threads | dev | king_avx2 (AVX2) | king_avx2 (scalar) | AVX2 / scalar vs dev |
|---|---|---|---|---|---|
| full | 1 | 6.59 (6.59-6.60) | 3.63 (3.63-3.63) | 5.11 (5.10-5.12) | 1.82x / 1.29x |
| grouped | 1 | 0.35 (0.35-0.35) | 0.32 (0.32-0.32) | 0.36 (0.36-0.36) | 1.11x / 0.99x |
| full | 8 | 0.94 (0.94-0.95) | 0.62 (0.62-0.62) | 0.77 (0.76-0.77) | 1.52x / 1.22x |
| grouped | 8 | 0.06 (0.06-0.06) | 0.06 (0.06-0.16) | 0.06 (0.06-0.10) | 1.04x / 0.93x |

Outputs (check_proto.R):

- check_king.rds identical to dev 
- check_king_scalar.rds identical to dev 

## gt_merge_vcfs() (PR #240)

| scenario | threads | dev | vcf_merge_speed | speed-up | same .bk and sidecar |
|---|---|---|---|---|---|
| pops | 1 | 7.46 (7.44-7.52) | 6.18 (6.17-6.18) | 1.21x | yes |
| sgdp | 1 | 38.97 (38.94-39.31) | 24.27 (24.24-24.29) | 1.61x | yes |
| many | 1 | 28.92 (28.91-29.53) | 21.69 (21.68-21.72) | 1.33x | yes |
| pops | 8 | 3.11 (3.11-3.14) | 2.75 (2.75-2.77) | 1.13x | yes |
| sgdp | 8 | 16.77 (16.69-16.80) | 10.62 (10.59-10.65) | 1.58x | yes |
| many | 8 | 18.24 (18.21-18.28) | 11.53 (11.51-11.57) | 1.58x | yes |

## gen_tibble() on one VCF (PR #240)

| file | threads | dev | vcf_merge_speed | speed-up | same .bk |
|---|---|---|---|---|---|
| pop1.vcf.gz | 1 | 1.38 (1.38-1.38) | 1.39 (1.39-1.45) | 0.99x | yes |
| S001.vcf.gz | 1 | 0.77 (0.76-0.77) | 0.75 (0.75-0.75) | 1.02x | yes |
| pop1.vcf.gz | 8 | 0.64 (0.63-0.66) | 0.60 (0.59-0.62) | 1.06x | yes |
| S001.vcf.gz | 8 | 0.53 (0.53-0.54) | 0.48 (0.47-0.49) | 1.10x | yes |

Commits: dev=ce748638, vcf=2f4dc9b6 
Load (1 min) before the merge runs: 0.95 (0.58-2.10) 
