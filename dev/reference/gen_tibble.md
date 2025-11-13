# Constructor for a `gen_tibble`

A `gen_tibble` stores genotypes for individuals in a tidy format.
DESCRIBE here the format

## Usage

``` r
gen_tibble(
  x,
  ...,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  quiet = FALSE
)

# S3 method for class 'character'
gen_tibble(
  x,
  ...,
  parser = c("cpp", "vcfR"),
  n_cores = 1,
  chunk_size = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  quiet = FALSE
)

# S3 method for class 'matrix'
gen_tibble(
  x,
  indiv_meta,
  loci,
  ...,
  ploidy = 2,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  quiet = FALSE
)
```

## Arguments

- x:

  can be:

  - a string giving the path to a PLINK BED or PED file. The associated
    BIM and FAM files for the BED, or MAP for PED are expected to be in
    the same directory and have the same file name.

  - a string giving the path to a RDS file storing a `bigSNP` object
    from the `bigsnpr` package (usually created with
    [`bigsnpr::snp_readBed()`](https://privefl.github.io/bigsnpr/reference/snp_readBed.html))

  - a string giving the path to a vcf file. Only biallelic SNPs will be
    considered.

  - a string giving the path to a *packedancestry* .geno file. The
    associated .ind and .snp files are expected to be in the same
    directory and share the same file name prefix.

  - a genotype matrix of dosages (0, 1, 2, NA) giving the dosage of the
    alternate allele.

- ...:

  if `x` is the name of a vcf file, additional arguments passed to
  [`vcfR::read.vcfR()`](https://rdrr.io/pkg/vcfR/man/io_vcfR.html).
  Otherwise, unused.

- valid_alleles:

  a vector of valid allele values; it defaults to 'A','T', 'C' and 'G'.

- missing_alleles:

  a vector of values in the BIM file/loci dataframe that indicate a
  missing value for the allele value (e.g. when we have a monomorphic
  locus with only one allele). It defaults to '0' and '.' (the same as
  PLINK 1.9).

- backingfile:

  the path, including the file name without extension, for backing files
  used to store the data (they will be given a .bk and .RDS
  automatically). This is not needed if `x` is already an .RDS file. If
  `x` is a .BED or a *VCF* file and `backingfile` is left NULL, the
  backing file will be saved in the same directory as the bed or vcf
  file, using the same file name but with a different file type (.bk
  rather than .bed or .vcf). If `x` is a genotype matrix and
  `backingfile` is NULL, then a temporary file will be created (but note
  that R will delete it at the end of the session!)

- allow_duplicates:

  logical. If TRUE, the tibble will allow duplicated loci (those with
  genomic coordinate (chromosome + position) or locus name appearing
  more than once). If FALSE, an error will be thrown if duplicated loci
  are found. These validations run before backing files are saved.
  Default is FALSE.

- quiet:

  provide information on the files used to store the data

- parser:

  the name of the parser used for *VCF*, either "cpp" to use a fast C++
  parser (the default), or "vcfR" to use the R package `vcfR`. The
  latter is slower but more robust; if "cpp" gives an error, try using
  "vcfR" in case your *VCF* has an unusual structure.

- n_cores:

  the number of cores to use for parallel processing

- chunk_size:

  the number of loci or individuals (depending on the format) processed
  at a time (currently used if `x` is a vcf with parser "vcfR")

- indiv_meta:

  a list, data.frame or tibble with compulsory columns 'id' and
  'population', plus any additional metadata of interest. This is only
  used if `x` is a genotype matrix. Otherwise this information is
  extracted directly from the files.

- loci:

  a data.frame or tibble, with compulsory columns 'name', 'chromosome',
  and 'position','genetic_dist', 'allele_ref' and 'allele_alt'. This is
  only used if `x` is a genotype matrix. Otherwise this information is
  extracted directly from the files.

- ploidy:

  the ploidy of the samples (either a single value, or a vector of
  values for mixed ploidy). Only used if creating a gen_tibble from a
  matrix of data; otherwise, ploidy is determined automatically from the
  data as they are read.

## Value

an object of the class `gen_tbl`.

## Details

- *VCF* files: the fast `cpp` parser is used by default. Both `cpp` and
  `vcfR` parsers attempt to establish ploidy from the first variant; if
  that variant is found in a sex chromosome (or mtDNA), the parser will
  fail with 'Error: a genotype has more than max_ploidy alleles...'. To
  successful import such a *VCF*, change the order of variants so that
  the first chromosome is an autosome using a tool such as `vcftools`.
  Currently, only biallelic SNPs are supported. If haploid variants
  (e.g. sex chromosomes) are included in the *VCF*, they are not
  transformed into homozygous calls. Instead, reference alleles will be
  coded as 0 and alternative alleles will be coded as 1.

- *packedancestry* files: When loading *packedancestry* files, missing
  alleles will be converted from 'X' to NA

## Note

Helper functions for accessing `gen_tibble` object attributes and
checking gen_tibble ploidy can be found in gt_helper_functions.R

## Examples

``` r
# Create a gen_tibble from a .bed file
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)
#> # A gen_tibble: 79 loci
#> # A tibble:     176 × 3
#>    id    population  genotypes
#>    <chr> <chr>      <vctr_SNP>
#>  1 Ale04 Ale         [0,.,...]
#>  2 Ale05 Ale         [1,0,...]
#>  3 Ale06 Ale         [.,0,...]
#>  4 Ale08 Ale         [.,2,...]
#>  5 Ale13 Ale         [0,.,...]
#>  6 Ale15 Ale         [1,0,...]
#>  7 Ale16 Ale         [0,.,...]
#>  8 Ale17 Ale         [1,.,...]
#>  9 Ale18 Ale         [0,0,...]
#> 10 Ale19 Ale         [2,.,...]
#> # ℹ 166 more rows

# Create a gen_tibble from a .vcf file
vcf_path <-
  system.file("extdata", "anolis",
    "punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
#> # A gen_tibble: 3249 loci
#> # A tibble:     46 × 2
#>    id                 genotypes
#>    <chr>             <vctr_SNP>
#>  1 punc_BM288         [0,0,...]
#>  2 punc_GN71          [2,0,...]
#>  3 punc_H1907         [0,2,...]
#>  4 punc_H1911         [0,2,...]
#>  5 punc_H2546         [0,1,...]
#>  6 punc_IBSPCRIB0361  [0,0,...]
#>  7 punc_ICST764       [0,0,...]
#>  8 punc_JFT459        [0,0,...]
#>  9 punc_JFT773        [0,0,...]
#> 10 punc_LG1299        [0,0,...]
#> # ℹ 36 more rows

# Create a gen_tibble from a matrix of genotypes:
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  valid_alleles = c("A", "T", "C", "G"),
  quiet = TRUE
)
#> # A gen_tibble: 6 loci
#> # A tibble:     3 × 3
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 b     pop1        [2,1,...]
#> 3 c     pop2        [2,2,...]
```
