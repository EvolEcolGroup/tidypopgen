# Encoding of phased genotypes in a `gen_tibble`

Phased diploid genotypes (read with `phased = TRUE` in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md))
are stored in the same file backed matrix (a
[bigstatsr::FBM.code256](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html))
as unphased ones, using one byte per genotype:

## Usage

``` r
code_phased

phased_bytes
```

## Format

`code_phased` is a numeric vector of length 256. `phased_bytes` is an
integer matrix with 256 rows (bytes 0 to 255) and two columns (`hap1`
and `hap2`).

## Details

|      |          |                    |
|------|----------|--------------------|
| byte | genotype | dosage (`code256`) |
| 0    | 0\|0     | 0                  |
| 1    | 0\|1     | 1                  |
| 2    | 1\|1     | 2                  |
| 3    | missing  | `NA`               |
| 4    | 1\|0     | 1                  |

Bytes 0 to 3 have the same meaning as in unphased diploid data (where
byte 3 is the missing value), so that code that reads or writes the raw
bytes keeps working; only the second orientation of a heterozygote
(`1|0`) uses a new byte. The `code256` of the matrix, `code_phased`,
translates each byte into a dosage, so that all functions that decode
genotypes through `code256` see phased data as ordinary unphased
dosages. `code_phased` is also how phased data are recognised
([`is_phased()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/is_phased.md)):
it differs from the codes used for unphased
([bigsnpr::CODE_012](https://privefl.github.io/bigsnpr/reference/snp_codes.html)),
imputed
([bigsnpr::CODE_IMPUTE_PRED](https://privefl.github.io/bigsnpr/reference/snp_codes.html))
and dosage
([bigsnpr::CODE_DOSAGE](https://privefl.github.io/bigsnpr/reference/snp_codes.html))
data, as it is the only one that decodes byte 4 as 1.

`phased_bytes` gives, for each byte, the alleles carried by the first
and second haplotype of an individual (`NA` for missing, and for bytes
that are not used).

## See also

[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
[`is_phased()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/is_phased.md),
[`show_haplotypes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_haplotypes.md)

## Examples

``` r
code_phased[1:6]
#> [1]  0  1  2 NA  1 NA
phased_bytes[1:6, ]
#>      hap1 hap2
#> [1,]    0    0
#> [2,]    0    1
#> [3,]    1    1
#> [4,]   NA   NA
#> [5,]    1    0
#> [6,]   NA   NA
```
