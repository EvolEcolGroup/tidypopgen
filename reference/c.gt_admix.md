# Combine method for gt_admix objects

Combine method for gt_admix objects

## Usage

``` r
# S3 method for class 'gt_admix'
c(..., match_attributes = TRUE)
```

## Arguments

- ...:

  A list of `gt_admix` objects

- match_attributes:

  boolean, determining whether all attributes (id, group and algorithm)
  of the `gt_admix` objects to be combined must be an exact match (TRUE,
  the default), or whether non-matching attributes should be ignored
  (FALSE)

## Value

A `gt_admix` object with the combined data

## Examples

``` r
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")

  # Create a second gt_admix object
  admix_obj2 <- example_gt %>% gt_snmf(k = 2:4, project = "force")

  # Combine the two gt_admix objects
  new_admix_obj <- c(admix_obj, admix_obj2)
  summary(new_admix_obj)
}
#> Admixture results for multiple runs:         
#> k 1 2 3 4
#> n 1 2 2 1
#> with slots:
#> $Q for Q matrices
#> $P for  matrices
#> $log for logs from the algorithm
```
