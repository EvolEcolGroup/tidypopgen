# A mutate method for grouped `gen_tibble` objects

A mutate method for grouped `gen_tibble` objects

## Usage

``` r
# S3 method for class 'grouped_gen_tbl'
mutate(..., deparse.level = 1)
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
test_gt %>% mutate(region = "East")
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 4
#> # Groups:       population [3]
#>   id    population  genotypes region
#>   <chr> <chr>      <vctr_SNP> <chr> 
#> 1 a     pop1        [1,1,...] East  
#> 2 b     pop1        [2,1,...] East  
#> 3 c     pop2        [2,.,...] East  
#> 4 d     pop2        [1,0,...] East  
#> 5 e     pop1        [1,2,...] East  
#> 6 f     pop3        [0,0,...] East  
#> 7 g     pop3        [0,1,...] East  
test_gt <- load_example_gt("grouped_gen_tbl_sf")
test_gt %>% mutate(region = "East")
#> Simple feature collection with 7 features and 6 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 41 xmax: 2 ymax: 51
#> Geodetic CRS:  WGS 84
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 7
#> # Groups:       population [3]
#>   id    population longitude latitude genotypes    geometry region
#> * <chr> <chr>          <dbl>    <dbl> <vctr_SN> <POINT [°]> <chr> 
#> 1 a     pop1               0       51 [1,1,...]      (0 51) East  
#> 2 b     pop1               0       51 [2,1,...]      (0 51) East  
#> 3 c     pop2               2       49 [2,.,...]      (2 49) East  
#> 4 d     pop2               2       49 [1,0,...]      (2 49) East  
#> 5 e     pop1               0       51 [1,2,...]      (0 51) East  
#> 6 f     pop3               2       41 [0,0,...]      (2 41) East  
#> 7 g     pop3               2       41 [0,1,...]      (2 41) East  
```
