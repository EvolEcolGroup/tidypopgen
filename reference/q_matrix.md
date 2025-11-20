# Convert a standard matrix to a `q_matrix` object

Takes a single Q matrix that exists as either a matrix or a data frame
and returns a `q_matrix` object.

## Usage

``` r
q_matrix(x)
```

## Arguments

- x:

  A matrix or a data frame

## Value

A `q_matrix` object

## Examples

``` r
# Read in a single .Q file
q_mat <- read.table(system.file("extdata", "anolis", "anolis_ld_run1.3.Q",
  package = "tidypopgen"
))
class(q_mat)
#> [1] "data.frame"

# Convert to a Q matrix object
q_mat <- q_matrix(q_mat)
class(q_mat)
#> [1] "q_matrix" "matrix"   "array"   
```
