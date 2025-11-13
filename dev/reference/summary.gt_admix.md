# Summary method for gt_admix objects

Summary method for gt_admix objects

## Usage

``` r
# S3 method for class 'gt_admix'
summary(object, ...)
```

## Arguments

- object:

  a `gt_admix` object

- ...:

  unused (necessary for compatibility with generic function)

## Value

A summary of the `gt_admix` object

## Examples

``` r
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")

  # Print a summary
  summary(admix_obj)
}
#> Admixture results for multiple runs:       
#> k 1 2 3
#> n 1 1 1
#> with slots:
#> $Q for Q matrices
#> $P for  matrices
#> $log for logs from the algorithm
```
