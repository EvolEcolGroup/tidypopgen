# Changelog

## tidypopgen 0.4.3

CRAN release: 2026-01-23

- fix occasional test failure whentemporary file names contained certain
  patterns

## tidypopgen 0.4.2

CRAN release: 2026-01-22

- allow least squares projection of multiple components at once in
  `gt_project_pca()`

## tidypopgen 0.4.1

CRAN release: 2025-12-20

- duplicate-distance checks removed in
  [`is_loci_table_ordered()`](https://evolecolgroup.github.io/tidypopgen/reference/is_loci_table_ordered.md)
- make `qc_report_*` functions work with pseudohaploids

## tidypopgen 0.4.0

CRAN release: 2025-10-24

- store chromosomes as factors, and remove redundant column storing
  chromosome as integer

## tidypopgen 0.3.3

- only use the FBM backend (rather than a complete bigsnpr object)
- fix a few minor bugs

## tidypopgen 0.3.2

CRAN release: 2025-08-27

- tidy up for CRAN submission

## tidypopgen 0.3.0

- add
  [`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pseudohaploid.md)
  and adapt appropriate functions to work with pseudohaploids

## tidypopgen 0.2.0

- Lots of optimisations to speed up processing of large datasets.
- Added an `sf` backend to `tidypopgen` for spatial data.
- Windowed statistics.

## tidypopgen 0.1.0

- Initial public version of `tidypopgen` released on r-universe.
