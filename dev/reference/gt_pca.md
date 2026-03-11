# Principal Component Analysis for `gen_tibble` objects

There are a number of PCA methods available for `gen_tibble` objects.
They are mostly designed to work on very large datasets, so they only
compute a limited number of components. For smaller datasets,
[`gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_partialSVD.md)
allows the use of partial (truncated) SVD to fit the PCA; this method is
suitable when the number of individuals is much smaller than the number
of loci. For larger dataset,
[`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_randomSVD.md)
is more appropriate. Finally, there is a method specifically designed
for dealing with LD in large datasets,
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md).
Whilst this is arguably the best option, it is somewhat data hungry, and
so only suitable for very large datasets (hundreds of individuals with
several hundred thousands markers, or larger).

## Note

Monomorphic markers must be removed before PCA is computed. The error
message 'Error: some variables have zero scaling; remove them before
attempting to scale.' indicates that monomorphic markers are present.

Using gt_pca_autoSVD with a small dataset will likely cause an error,
see man page for details.

Rather than accessing elements of a `gt_pca` object directly, it is
better to use `tidy` and `augment`. See
[`gt_pca_tidiers`](https://evolecolgroup.github.io/tidypopgen/dev/reference/tidy_gt_pca.md).

## See also

Other gt_pca_functions:
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md),
[`gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_partialSVD.md),
[`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_randomSVD.md)
