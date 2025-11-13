# Augment data with information from a q_matrix object

Augment for `q_matrix` accepts a model object and a dataset and adds Q
values to each observation in the dataset. Q values are stored in
separate columns, which is given name with the pattern ".Q1",".Q2", etc.
For consistency with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html),
a column ".rownames" is also returned; it is a copy of 'id', but it
ensures that any scripts written for data augmented with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html)
will work out of the box (this is especially helpful when adapting
plotting scripts).

## Usage

``` r
# S3 method for class 'q_matrix'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  A `q_matrix` object

- data:

  the `gen_tibble` used to run the clustering algorithm

- ...:

  Not used. Needed to match generic signature only.

## Value

A
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
containing the original data along with additional columns containing
each observation's Q values.

## Examples

``` r
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")

  # Extract a Q matrix
  q_mat_k3 <- get_q_matrix(admix_obj, k = 3, run = 1)

  # Augment the gen_tibble with Q values
  augment(q_mat_k3, data = example_gt)
}
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 7
#>   id    population  genotypes .rownames .Q1        .Q2        .Q3       
#>   <chr> <chr>      <vctr_SNP> <chr>     <q_matrix> <q_matrix> <q_matrix>
#> 1 a     pop1        [1,1,...] a         9.998e-01  9.999e-05  9.999e-05 
#> 2 b     pop1        [2,1,...] b         9.999e-05  9.998e-01  9.999e-05 
#> 3 c     pop2        [2,.,...] c         9.999e-05  9.999e-05  9.998e-01 
#> 4 d     pop2        [1,0,...] d         9.998e-01  9.999e-05  9.999e-05 
#> 5 e     pop1        [1,2,...] e         9.998e-01  9.999e-05  9.999e-05 
#> 6 f     pop3        [0,0,...] f         9.999e-05  9.999e-05  9.998e-01 
#> 7 g     pop3        [0,1,...] g         9.999e-05  9.999e-05  9.998e-01 
```
