# tidypopgen: Tidy Population Genetics

We provide a tidy grammar of population genetics, facilitating the
manipulation and analysis of data on biallelic single nucleotide
polymorphisms (SNPs). 'tidypopgen' scales to very large genetic datasets
by storing genotypes on disk, and performing operations on them in
chunks, without ever loading all data in memory. The full
functionalities of the package are described in Carter et al. (2025)
[doi:10.1101/2025.06.06.658325](https://doi.org/10.1101/2025.06.06.658325)
.

## See also

Useful links:

- <https://github.com/EvolEcolGroup/tidypopgen>

- <https://evolecolgroup.github.io/tidypopgen/>

- Report bugs at <https://github.com/EvolEcolGroup/tidypopgen/issues>

## Author

**Maintainer**: Andrea Manica <am315@cam.ac.uk>
([ORCID](https://orcid.org/0000-0003-1895-450X)) \[copyright holder\]

Authors:

- Evie Carter

- Eirlys Tysall

Other contributors:

- Chang Christopher (Author of Hardy-Weinberg Equilibrium algorithm in
  PLINK 1.90, used in loci_hwe()) \[contributor\]

- Shaun Purcell (Author of Hardy-Weinberg Equilibrium algorithm in PLINK
  1.90, used in loci_hwe()) \[contributor\]

- Bengtsson Henrik (Author of countLines in R.utils, modified for .vcf
  in count_vcf_variants()) \[contributor\]
