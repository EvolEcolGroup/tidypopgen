# Return a single Q matrix from a `gt_admix` object

This function retrieves a single Q matrix from a `gt_admix` object based
on the specified k value and run number.

## Usage

``` r
get_q_matrix(x, ..., k, run)
```

## Arguments

- x:

  A `gt_admix` object containing multiple Q matrices

- ...:

  Not used

- k:

  The k value of the desired Q matrix

- run:

  The run number of the desired Q matrix

## Value

A single Q matrix from the `gt_admix` object

## Examples

``` r
# Read example gt_admix obejct
admix_obj <-
  readRDS(system.file("extdata", "anolis", "anole_adm_k3.rds",
    package = "tidypopgen"
  ))

# Extract a Q matrix
get_q_matrix(admix_obj, k = 3, run = 1)
#>            .Q1      .Q2      .Q3
#>  [1,] 0.000010 0.000010 0.999980
#>  [2,] 0.000010 0.000010 0.999980
#>  [3,] 0.000010 0.999980 0.000010
#>  [4,] 0.000010 0.999980 0.000010
#>  [5,] 0.000010 0.999980 0.000010
#>  [6,] 0.999980 0.000010 0.000010
#>  [7,] 0.754153 0.000010 0.245837
#>  [8,] 0.999980 0.000010 0.000010
#>  [9,] 0.999980 0.000010 0.000010
#> [10,] 0.999980 0.000010 0.000010
#> [11,] 0.000010 0.999980 0.000010
#> [12,] 0.000010 0.967888 0.032102
#> [13,] 0.000010 0.224456 0.775534
#> [14,] 0.000010 0.999980 0.000010
#> [15,] 0.000010 0.000010 0.999980
#> [16,] 0.000010 0.999980 0.000010
#> [17,] 0.000247 0.000010 0.999743
#> [18,] 0.000010 0.999980 0.000010
#> [19,] 0.076782 0.000010 0.923208
#> [20,] 0.000010 0.000010 0.999980
#> [21,] 0.332415 0.000010 0.667575
#> [22,] 0.000010 0.000010 0.999980
#> [23,] 0.000010 0.000010 0.999980
#> [24,] 0.000010 0.000010 0.999980
#> [25,] 0.999980 0.000010 0.000010
#> [26,] 0.999980 0.000010 0.000010
#> [27,] 0.999980 0.000010 0.000010
#> [28,] 0.999980 0.000010 0.000010
#> [29,] 0.999980 0.000010 0.000010
#> [30,] 0.000016 0.999974 0.000010
#> [31,] 0.000010 0.000010 0.999980
#> [32,] 0.039869 0.000010 0.960121
#> [33,] 0.999980 0.000010 0.000010
#> [34,] 0.000010 0.999980 0.000010
#> [35,] 0.000010 0.000010 0.999980
#> [36,] 0.000010 0.221311 0.778679
#> [37,] 0.003767 0.265606 0.730627
#> [38,] 0.999980 0.000010 0.000010
#> [39,] 0.999980 0.000010 0.000010
#> [40,] 0.000010 0.000010 0.999980
#> [41,] 0.271323 0.000010 0.728667
#> [42,] 0.999980 0.000010 0.000010
#> [43,] 0.999980 0.000010 0.000010
#> [44,] 0.743485 0.000010 0.256505
#> [45,] 0.000010 0.993175 0.006815
#> [46,] 0.000010 0.999980 0.000010
#> attr(,"class")
#> [1] "q_matrix" "q_matrix" "matrix"   "array"   
#> attr(,"id")
#>  [1] "punc_BM288"        "punc_GN71"         "punc_H1907"       
#>  [4] "punc_H1911"        "punc_H2546"        "punc_IBSPCRIB0361"
#>  [7] "punc_ICST764"      "punc_JFT459"       "punc_JFT773"      
#> [10] "punc_LG1299"       "punc_LSUMZH12577"  "punc_LSUMZH12751" 
#> [13] "punc_LSUMZH13910"  "punc_LSUMZH14100"  "punc_LSUMZH14336" 
#> [16] "punc_LSUMZH15476"  "punc_MPEG20846"    "punc_MPEG21348"   
#> [19] "punc_MPEG22415"    "punc_MPEG24758"    "punc_MPEG26102"   
#> [22] "punc_MPEG28489"    "punc_MPEG29314"    "punc_MPEG29943"   
#> [25] "punc_MTR05978"     "punc_MTR12338"     "punc_MTR12511"    
#> [28] "punc_MTR15267"     "punc_MTR17744"     "punc_MTR18550"    
#> [31] "punc_MTR20798"     "punc_MTR21474"     "punc_MTR21545"    
#> [34] "punc_MTR25584"     "punc_MTR28048"     "punc_MTR28401"    
#> [37] "punc_MTR28593"     "punc_MTR34227"     "punc_MTR34414"    
#> [40] "punc_MTR976723"    "punc_MTR978312"    "punc_MTRX1468"    
#> [43] "punc_MTRX1478"     "punc_MUFAL9635"    "punc_PJD409"      
#> [46] "punc_UNIBAN1670"  
#> attr(,"group")
#>  [1] "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#>  [5] "Amazonian_Forest" "Atlantic_Forest"  "Atlantic_Forest"  "Atlantic_Forest" 
#>  [9] "Atlantic_Forest"  "Atlantic_Forest"  "Amazonian_Forest" "Amazonian_Forest"
#> [13] "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#> [17] "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#> [21] "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#> [25] "Atlantic_Forest"  "Atlantic_Forest"  "Atlantic_Forest"  "Atlantic_Forest" 
#> [29] "Atlantic_Forest"  "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#> [33] "Atlantic_Forest"  "Amazonian_Forest" "Amazonian_Forest" "Amazonian_Forest"
#> [37] "Amazonian_Forest" "Atlantic_Forest"  "Atlantic_Forest"  "Amazonian_Forest"
#> [41] "Amazonian_Forest" "Atlantic_Forest"  "Atlantic_Forest"  "Atlantic_Forest" 
#> [45] "Amazonian_Forest" "Amazonian_Forest"
```
