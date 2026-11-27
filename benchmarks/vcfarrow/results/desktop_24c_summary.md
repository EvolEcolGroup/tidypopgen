# tidypopgen on the VCFArrow benchmark

Synthetic DiscoSnpRad VCF (56 samples), one fresh R process per step, peak PSS of the process tree (data loaded + step).

## 2,322,460 variants, 4 threads

| step | v1.0.0 time | v1.0.0 peak GB | vcf_fast time | vcf_fast peak GB | speed_up | memory | same_result |
|---|---|---|---|---|---|---|---|
| read_vcf | 12.9 s | 0.78 | 4.4 s | 0.57 | 2.9x | 0.73x | yes |
| set_groups | 0.0 s | 0.79 | 0.0 s | 0.57 | 2.0x | 0.73x | yes |
| extract_samples | 0.0 s | 0.46 | 0.0 s | 0.46 | 1.0x | 1.00x | yes |
| extract_groups | 0.0 s | 0.46 | 0.0 s | 0.46 | 1.0x | 1.00x | yes |
| subsample_loci | 0.0 s | 0.46 | 0.0 s | 0.46 | 1.2x | 1.00x | yes |
| filter_maf | 0.2 s | 0.68 | 0.2 s | 0.46 | 1.1x | 0.68x | yes |
| filter_locus_missing | 0.1 s | 0.46 | 0.1 s | 0.46 | 1.1x | 1.00x | yes |
| filter_sample_missing | 0.1 s | 0.46 | 0.1 s | 0.46 | 1.0x | 1.00x | yes |
| filter_invariant | 0.2 s | 0.46 | 0.2 s | 0.46 | 1.1x | 1.00x | yes |
| merge | 7.3 s | 1.32 | 4.7 s | 1.32 | 1.5x | 1.00x | yes |
| sample_stats | 0.1 s | 0.46 | 0.1 s | 0.46 | 1.1x | 1.00x | yes |
| theta_pi | 0.0 s | 0.46 | 0.0 s | 0.46 | 1.0x | 1.00x | yes |
| missing_plot | 0.3 s | 0.67 | 0.3 s | 0.71 | 1.2x | 1.05x | yes |
| write_vcf | NA | NA | 5.8 s | 0.88 |  |  |  |
| to_genlight | 12.9 s | 2.71 | 9.1 s | 1.79 | 1.4x | 0.66x | yes |
| to_ped | 2.4 s | 0.91 | 2.3 s | 0.91 | 1.0x | 1.00x | yes |
| to_bed | 0.7 s | 0.73 | 0.7 s | 0.69 | 1.0x | 0.95x | yes |
| subsample_loci_bare | 7.3 s | 1.19 | 7.4 s | 1.19 | 1.0x | 1.00x | yes |
| theta_pi_grouped | 0.4 s | 0.71 | 0.4 s | 0.71 | 1.0x | 1.00x | yes |
| filter_chain_update | 2.6 s | 0.84 | 2.0 s | 0.84 | 1.3x | 1.00x | yes |
| read_vcf_int | 9.5 s | 0.87 | 2.4 s | 0.63 | 3.9x | 0.73x | yes |

## 4,644,921 variants, 4 threads

| step | v1.0.0 time | v1.0.0 peak GB | vcf_fast time | vcf_fast peak GB | speed_up | memory | same_result |
|---|---|---|---|---|---|---|---|
| read_vcf | 24.9 s | 1.46 | 10.4 s | 0.94 | 2.4x | 0.64x | yes |
| set_groups | 0.0 s | 1.44 | 0.0 s | 0.89 | 0.5x | 0.62x | yes |
| extract_samples | 0.0 s | 0.82 | 0.0 s | 0.82 | 1.0x | 1.00x | yes |
| extract_groups | 0.0 s | 0.82 | 0.0 s | 0.82 | 1.8x | 1.00x | yes |
| subsample_loci | 0.0 s | 0.84 | 0.0 s | 0.82 | 1.0x | 0.98x | yes |
| filter_maf | 0.9 s | 1.21 | 1.1 s | 1.21 | 0.8x | 1.00x | yes |
| filter_locus_missing | 0.1 s | 1.02 | 0.1 s | 0.82 | 0.9x | 0.80x | yes |
| filter_sample_missing | 0.1 s | 0.82 | 0.1 s | 0.82 | 1.0x | 1.00x | yes |
| filter_invariant | 0.9 s | 1.23 | 1.1 s | 1.23 | 0.8x | 1.00x | yes |
| merge | 13.0 s | 2.59 | 10.5 s | 2.59 | 1.2x | 1.00x | yes |
| sample_stats | 0.1 s | 0.85 | 0.1 s | 0.82 | 1.0x | 0.96x | yes |
| theta_pi | 0.1 s | 0.82 | 0.1 s | 0.82 | 1.0x | 1.00x | yes |
| missing_plot | 0.3 s | 0.82 | 0.3 s | 0.82 | 1.0x | 1.00x | yes |
| write_vcf | NA | NA | 11.9 s | 1.59 |  |  |  |
| to_genlight | 25.4 s | 5.06 | 20.5 s | 2.95 | 1.2x | 0.58x | yes |
| to_ped | 5.0 s | 1.62 | 5.3 s | 1.62 | 0.9x | 1.00x | yes |
| to_bed | 1.0 s | 1.12 | 1.1 s | 1.47 | 0.9x | 1.32x | yes |
| subsample_loci_bare | 25.3 s | 2.28 | 32.9 s | 2.29 | 0.8x | 1.00x | yes |
| theta_pi_grouped | 1.3 s | 1.30 | 1.5 s | 1.37 | 0.8x | 1.05x | yes |
| filter_chain_update | 5.9 s | 1.51 | 5.2 s | 1.56 | 1.1x | 1.03x | yes |
| read_vcf_int | 20.1 s | 1.62 | 5.9 s | 1.10 | 3.4x | 0.68x | yes |

## 18,579,683 variants, 4 threads

| step | v1.0.0 time | v1.0.0 peak GB | vcf_fast time | vcf_fast peak GB | speed_up | memory | same_result |
|---|---|---|---|---|---|---|---|
| read_vcf | 1.8 min | 5.59 | 36.5 s | 3.56 | 3.0x | 0.64x | yes |
| set_groups | 0.0 s | 5.36 | 0.0 s | 2.99 | 1.0x | 0.56x | yes |
| extract_samples | 0.0 s | 2.96 | 0.0 s | 2.96 | 1.0x | 1.00x | yes |
| extract_groups | 0.0 s | 2.96 | 0.0 s | 2.96 | 1.2x | 1.00x | yes |
| subsample_loci | 0.0 s | 2.96 | 0.0 s | 2.96 | 1.0x | 1.00x | yes |
| filter_maf | 3.2 s | 4.52 | 3.1 s | 4.52 | 1.0x | 1.00x | yes |
| filter_locus_missing | 0.9 s | 4.44 | 0.9 s | 4.44 | 1.1x | 1.00x | yes |
| filter_sample_missing | 0.5 s | 3.72 | 0.5 s | 3.98 | 1.0x | 1.07x | yes |
| filter_invariant | 3.1 s | 4.51 | 3.1 s | 4.74 | 1.0x | 1.05x | yes |
| merge | 53.1 s | 9.03 | 38.6 s | 9.03 | 1.4x | 1.00x | yes |
| sample_stats | 0.5 s | 3.83 | 0.5 s | 3.83 | 1.0x | 1.00x | yes |
| theta_pi | 0.3 s | 2.96 | 0.3 s | 2.96 | 1.0x | 1.00x | yes |
| missing_plot | 0.7 s | 3.86 | 0.7 s | 3.99 | 1.0x | 1.03x | yes |
| write_vcf | NA | NA | 47.1 s | 6.27 |  |  |  |
| to_genlight | 2.0 min | 18.96 | 1.1 min | 9.54 | 1.9x | 0.50x | yes |
| to_ped | 20.1 s | 5.46 | 19.8 s | 5.46 | 1.0x | 1.00x | yes |
| to_bed | 2.2 s | 4.72 | 2.4 s | 4.76 | 0.9x | 1.01x | yes |
| subsample_loci_bare | 5.9 min | 8.90 | 5.9 min | 8.90 | 1.0x | 1.00x | yes |
| theta_pi_grouped | 4.8 s | 4.88 | 4.7 s | 5.07 | 1.0x | 1.04x | yes |
| filter_chain_update | 22.5 s | 5.62 | 17.4 s | 5.62 | 1.3x | 1.00x | yes |
| read_vcf_int | 1.4 min | 6.24 | 21.3 s | 4.24 | 3.9x | 0.68x | yes |

## Estimated at 18.6 M variants (4 threads), with the published numbers

Measured on this machine (24 cores, 126 GB) at 2,322,460, 4,644,921, 18,579,683 variants; values at 18.6 M are measured where that size was run, otherwise extrapolated from the two largest sizes ("—": not enough sizes). The published runs used another machine, 50 GB, a 60-min limit, and one R session for all functions.

| function | published tidypopgen (4) | published tidypopgen mem | published VCFArrow (4) | v1.0.0 est. time | v1.0.0 est. peak GB | vcf_fast est. time | vcf_fast est. peak GB |
|---|---|---|---|---|---|---|---|
| read VCF | 30.0 min | 23.38/32.29 | 6.9 min | 1.8 min | 5.59 | 36.5 s | 3.56 |
| set groups | 0.0 s | 30.96/30.96 | 0.0 s | 0.0 s | 5.36 | 0.0 s | 2.99 |
| extract samples | 0.0 s | 30.96/30.96 | 55.6 s | 0.0 s | 2.96 | 0.0 s | 2.96 |
| extract groups | 0.0 s | 30.96/30.96 | 1.0 min | 0.0 s | 2.96 | 0.0 s | 2.96 |
| subsample loci (1000) | 6.1 min | 30.96/30.96 | 0.5 s | 0.0 s | 2.96 | 0.0 s | 2.96 |
| filter MAF (0.05) | 24.1 s | 30.96/38.25 | 33.7 s | 3.2 s | 4.52 | 3.1 s | 4.52 |
| filter locus missing (0.2) | 21.0 s | 30.96/38.25 | 19.0 s | 0.9 s | 4.44 | 0.9 s | 4.44 |
| filter sample missing (0.8) | 4.7 s | 30.96/30.96 | 3.6 min | 0.5 s | 3.72 | 0.5 s | 3.98 |
| filter invariant | 20.1 s | 30.96/38.25 | 41.8 s | 3.1 s | 4.51 | 3.1 s | 4.74 |
| merge two VCFs | 1.4 min | 32.65/38.25 | 10.3 min | 53.1 s | 9.03 | 38.6 s | 9.03 |
| sample statistics | 4.7 s | 30.96/33.82 | 1.0 min | 0.5 s | 3.83 | 0.5 s | 3.83 |
| theta / pi | 0.6 s | 30.96/30.96 | 1.1 min | 0.3 s | 2.96 | 0.3 s | 2.96 |
| missing data plot + table | 4.7 s | 30.96/30.96 | 46.0 s | 0.7 s | 3.86 | 0.7 s | 3.99 |
| write VCF | > 60 min | 30.96/30.96 | 8.5 min | — | — | 47.1 s | 6.27 |
| -> genlight | 3.3 min | 38.25/39.01 | 6.9 min | 2.0 min | 18.96 | 1.1 min | 9.54 |
| -> PLINK ped | > 60 min | 32.53/38.64 | 2.9 min | 20.1 s | 5.46 | 19.8 s | 5.46 |
| -> PLINK bed | 4.1 s | 30.84/31.10 | 1.7 min | 2.2 s | 4.72 | 2.4 s | 4.76 |

