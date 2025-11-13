# An arrange method for grouped `gen_tibble` objects

An arrange method for grouped `gen_tibble` objects

## Usage

``` r
# S3 method for class 'grouped_gen_tbl'
arrange(..., deparse.level = 1)
```

## Arguments

- ...:

  a gen_tibble and a data.frame or tibble

- deparse.level:

  an integer controlling the construction of column names.

## Value

a grouped `gen_tibble`

## Examples

``` r
test_gt <- load_example_gt("grouped_gen_tbl")
test_gt %>% arrange(id)
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 3
#> # Groups:       population [3]
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 b     pop1        [2,1,...]
#> 3 c     pop2        [2,.,...]
#> 4 d     pop2        [1,0,...]
#> 5 e     pop1        [1,2,...]
#> 6 f     pop3        [0,0,...]
#> 7 g     pop3        [0,1,...]
test_gt <- load_example_gt("grouped_gen_tbl_sf")
test_gt %>% arrange(id)
#> Simple feature collection with 7 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 41 xmax: 2 ymax: 51
#> Geodetic CRS:  WGS 84
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 6
#> # Groups:       population [3]
#>   id    population longitude latitude  genotypes    geometry
#>   <chr> <chr>          <dbl>    <dbl> <vctr_SNP> <POINT [°]>
#> 1 a     pop1               0       51  [1,1,...]      (0 51)
#> 2 b     pop1               0       51  [2,1,...]      (0 51)
#> 3 c     pop2               2       49  [2,.,...]      (2 49)
#> 4 d     pop2               2       49  [1,0,...]      (2 49)
#> 5 e     pop1               0       51  [1,2,...]      (0 51)
#> 6 f     pop3               2       41  [0,0,...]      (2 41)
#> 7 g     pop3               2       41  [0,1,...]      (2 41)
```
