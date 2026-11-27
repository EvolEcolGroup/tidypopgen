# tidypopgen on the VCFArrow benchmark

Synthetic DiscoSnpRad VCF (56 samples), one fresh R process per step, peak PSS of the process tree (data loaded + step).

## 200,000 variants, 4 threads

| step | v1.0.0 time | v1.0.0 peak GB | gt_as_vcf_speed time | gt_as_vcf_speed peak GB | speed_up | memory | same_result |
|---|---|---|---|---|---|---|---|
| read_vcf | 2.1 s | 0.19 | 1.1 s | 0.16 | 1.8x | 0.83x | yes |
| set_groups | 0.0 s | 0.19 | 0.0 s | 0.17 | 0.6x | 0.91x | yes |
| extract_samples | 0.0 s | 0.13 | 0.0 s | 0.15 | 1.0x | 1.09x | yes |
| extract_groups | 0.0 s | 0.13 | 0.0 s | 0.15 | 1.0x | 1.09x | yes |
| subsample_loci | 0.0 s | 0.15 | 0.0 s | 0.15 | 0.9x | 0.99x | yes |
| filter_maf | 0.0 s | 0.15 | 0.0 s | 0.15 | 1.0x | 1.01x | yes |
| filter_locus_missing | 0.0 s | 0.15 | 0.0 s | 0.15 | 0.9x | 1.00x | yes |
| filter_sample_missing | 0.0 s | 0.14 | 0.0 s | 0.15 | 1.3x | 1.03x | yes |
| filter_invariant | 0.0 s | 0.14 | 0.0 s | 0.15 | 1.0x | 1.04x | yes |
| merge | 1.6 s | 0.25 | 1.2 s | 0.23 | 1.3x | 0.92x | yes |
| sample_stats | 0.0 s | 0.15 | 0.0 s | 0.15 | 1.0x | 1.01x | yes |
| theta_pi | 0.0 s | 0.15 | 0.0 s | 0.15 | 1.0x | 0.99x | yes |
| missing_plot | 0.5 s | 0.23 | 0.5 s | 0.20 | 0.9x | 0.84x | yes |
| write_vcf | timeout (> 2 min) | 0.29 | 1.1 s | 0.20 |  | 0.68x |  |
| to_genlight | 3.9 s | 0.56 | 3.3 s | 0.66 | 1.2x | 1.17x | yes |
| to_ped | 0.4 s | 0.18 | 0.4 s | 0.22 | 0.9x | 1.21x | yes |
| to_bed | 1.0 s | 0.26 | 1.1 s | 0.22 | 1.0x | 0.87x | yes |
| subsample_loci_bare | 0.5 s | 0.19 | 0.5 s | 0.19 | 1.0x | 1.00x | yes |
| theta_pi_grouped | 0.0 s | 0.15 | 0.0 s | 0.14 | 1.0x | 0.94x | yes |
| filter_chain_update | 0.6 s | 0.19 | 0.6 s | 0.19 | 1.1x | 0.98x | yes |
| read_vcf_int | 1.6 s | 0.19 | 0.5 s | 0.16 | 3.0x | 0.86x | yes |

## Estimated at 18.6 M variants (4 threads), with the published numbers

Measured on this machine (4 cores, 15 GB) at 200,000 variants; values at 18.6 M are measured where that size was run, otherwise extrapolated from the two largest sizes ("—": not enough sizes). The published runs used another machine, 50 GB, a 60-min limit, and one R session for all functions.

| function | published tidypopgen (4) | published tidypopgen mem | published VCFArrow (4) | v1.0.0 est. time | v1.0.0 est. peak GB | gt_as_vcf_speed est. time | gt_as_vcf_speed est. peak GB |
|---|---|---|---|---|---|---|---|
| read VCF | 30.0 min | 23.38/32.29 | 6.9 min | — | — | — | — |
| set groups | 0.0 s | 30.96/30.96 | 0.0 s | — | — | — | — |
| extract samples | 0.0 s | 30.96/30.96 | 55.6 s | — | — | — | — |
| extract groups | 0.0 s | 30.96/30.96 | 1.0 min | — | — | — | — |
| subsample loci (1000) | 6.1 min | 30.96/30.96 | 0.5 s | — | — | — | — |
| filter MAF (0.05) | 24.1 s | 30.96/38.25 | 33.7 s | — | — | — | — |
| filter locus missing (0.2) | 21.0 s | 30.96/38.25 | 19.0 s | — | — | — | — |
| filter sample missing (0.8) | 4.7 s | 30.96/30.96 | 3.6 min | — | — | — | — |
| filter invariant | 20.1 s | 30.96/38.25 | 41.8 s | — | — | — | — |
| merge two VCFs | 1.4 min | 32.65/38.25 | 10.3 min | — | — | — | — |
| sample statistics | 4.7 s | 30.96/33.82 | 1.0 min | — | — | — | — |
| theta / pi | 0.6 s | 30.96/30.96 | 1.1 min | — | — | — | — |
| missing data plot + table | 4.7 s | 30.96/30.96 | 46.0 s | — | — | — | — |
| write VCF | > 60 min | 30.96/30.96 | 8.5 min | — | — | — | — |
| -> genlight | 3.3 min | 38.25/39.01 | 6.9 min | — | — | — | — |
| -> PLINK ped | > 60 min | 32.53/38.64 | 2.9 min | — | — | — | — |
| -> PLINK bed | 4.1 s | 30.84/31.10 | 1.7 min | — | — | — | — |

