# Autoplots for `q_matrix` objects

This autoplot will automatically rearrange individuals according to
their id and any grouping variables if an associated 'data' gen_tibble
is provided. To avoid any automatic re-sorting of individuals, set
`arrange_by_group` and `arrange_by_indiv` to FALSE.

## Usage

``` r
# S3 method for class 'q_matrix'
autoplot(
  object,
  data = NULL,
  annotate_group = TRUE,
  arrange_by_group = TRUE,
  arrange_by_indiv = TRUE,
  reorder_within_groups = FALSE,
  ...
)
```

## Arguments

- object:

  A Q matrix object (as returned by
  [`q_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/q_matrix.md)).

- data:

  An associated tibble (e.g. a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)),
  with the individuals in the same order as the data used to generate
  the Q matrix (unless the Q matrix has the ids of the individuals as an
  attribute, as from
  [`get_q_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/get_q_matrix.md)).
  If the Q matrix has no group information and `data` is grouped, the
  groups of `data` are used.

- annotate_group:

  Boolean determining whether to annotate the plot with the group
  information

- arrange_by_group:

  Boolean determining whether to arrange the individuals by group.
  Groups are ordered as by the grouped functions (see
  [`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html)):
  if the grouping variable in the `gen_tibble` or the metadata of the
  `gt_admix` object is a factor, by the levels of the factor; numbers in
  numeric order; otherwise alphabetically. With several grouping
  variables, groups are ordered by the first variable, then by the
  second, and so on.

- arrange_by_indiv:

  Boolean determining whether to arrange the individuals by their
  individual id within each group (groups are arranged first if
  `arrange_by_group` is TRUE, otherwise they are kept in their current
  order). This only applies if `annotate_group` is TRUE and group
  information is available. If `id` in the `get_tibble` or the metadata
  of the `gt_admix` object is a factor, it will be ordered by the levels
  of the factor; else it will be ordered alphabetically.

- reorder_within_groups:

  Boolean determining whether to reorder the individuals within each
  group based on their ancestry proportion (note that this is not
  advised if you are making multiple plots, as you would get a different
  order for each plot!). Groups keep the order set by
  `arrange_by_group`. If TRUE, `annotate_group` must also be TRUE.

- ...:

  not currently used.

## Value

a barplot of individuals, coloured by ancestry proportion

## Examples

``` r
# Read example gt_admix obejct
admix_obj <-
  readRDS(system.file("extdata", "anolis", "anole_adm_k3.rds",
    package = "tidypopgen"
  ))

# Extract a Q matrix
q_mat_k3 <- get_q_matrix(admix_obj, k = 3, run = 1)

# Basic autoplot
autoplot(q_mat_k3, annotate_group = FALSE, arrange_by_group = FALSE)


# To arrange individuals by group and by Q proportion
autoplot(q_mat_k3,
  annotate_group = TRUE, arrange_by_group = TRUE,
  arrange_by_indiv = TRUE, reorder_within_groups = TRUE
)

```
