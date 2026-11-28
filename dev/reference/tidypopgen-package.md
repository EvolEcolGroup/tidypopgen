# tidypopgen: Tidy Population Genetics

We provide a tidy grammar of population genetics, facilitating the
manipulation and analysis of data on biallelic single nucleotide
polymorphisms (SNPs). 'tidypopgen' scales to very large genetic datasets
by storing genotypes on disk, and performing operations on them in
chunks, without ever loading all data in memory. The full
functionalities of the package are described in Carter et al. (2025)
[doi:10.1111/2041-210x.70204](https://doi.org/10.1111/2041-210x.70204) .

## See also

Useful links:

- <https://github.com/EvolEcolGroup/tidypopgen>

- <https://evolecolgroup.github.io/tidypopgen/>

- Report bugs at <https://github.com/EvolEcolGroup/tidypopgen/issues>

## Author

**Maintainer**: Andrea Manica <am315@cam.ac.uk>
([ORCID](https://orcid.org/0000-0003-1895-450X)) \[copyright holder\]

Authors:

- Andrea Manica <am315@cam.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-1895-450X)) \[copyright holder\]

- Evie Carter

- Eirlys Tysall

Other contributors:

- Max Carter-Brown \[contributor\]

- Jason Hodgson \[contributor\]

- Christopher Chang (Author of PLINK 1.90: Hardy-Weinberg exact test
  code used in loci_hwe() and qc_report_loci(); runs of homozygosity
  algorithm followed by windows_indiv_roh(); LD pruning algorithm
  followed by loci_ld_prune()) \[contributor, copyright holder\]

- Shaun Purcell (Author of PLINK 1.90: Hardy-Weinberg exact test code
  used in loci_hwe() and qc_report_loci(); runs of homozygosity
  algorithm followed by windows_indiv_roh(); LD pruning algorithm
  followed by loci_ld_prune()) \[contributor, copyright holder\]

- Henrik Bengtsson (Author of countLines() in R.utils, modified for .vcf
  in count_vcf_variants()) \[contributor, copyright holder\]

- Jan Wigginton (Author of the original Hardy-Weinberg exact test code
  (SNPHWE), on which the PLINK 1.90 code used in loci_hwe() is based)
  \[contributor, copyright holder\]

- Florian Privé (Author of bigsnpr and bigstatsr, code from which is
  adapted in read_packedancestry_to_fbm(), gt_as_plink(),
  predict.gt_pca(), autoplot.gt_dapc() and loci_ld_clump())
  \[contributor, copyright holder\]

- Jérôme Goudet (Author of hierfstat, code from which is adapted in
  pop_global_stats(), pop_het_obs(), pop_het_exp() and
  indiv_inbreeding()) \[contributor, copyright holder\]

- Robert Maier (Author of admixtools, code from which is adapted in
  gt_extract_f2()) \[contributor, copyright holder\]

- Ryan Gutenkunst (Author of dadi, whose Spectrum class is translated to
  R in the tests of the site frequency spectra) \[contributor, copyright
  holder\]

- Isaac Overcast (Author of easySFS, translated to R in the tests of the
  site frequency spectra) \[contributor, copyright holder\]
