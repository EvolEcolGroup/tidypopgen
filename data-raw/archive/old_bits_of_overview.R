

When creating
summaries, the raw data are used by default:
  ```{r}
missing_gt %>% loci_missingness() %>% hist()
```

We can now run the analysis:
  ```{r}
missing_pca <- missing_gt %>% gt_pca_partialSVD()
missing_pca
```






(but it fails because we have too few SNPs!!!!):
  ```{r, error=TRUE}
missing_pca <- missing_gt %>% gt_pca_autoSVD()
```




We can use a standard PCA:
  ```{r}
missing_pca <- missing_gt %>% gt_pca_partialSVD()
missing_pca
```
Note that the observed genotypes have been kept, and we can go back to using
them (and thus have missing values) with:
  ```{r}
missing_gt <- missing_gt  %>% gt_set_imputed(FALSE)
```


# This should be changed to just a section on imputation, and move the details of analysis into its own vignette.

# Population genetics analyses
Besides using various metrics to describe populations, we often want to perform
more complex analyses, such as Principal Component Analysis or clustering.

We will us the example dataset from `bigsnpr`:
  ```{r}
library(tidypopgen)
library(bigsnpr)
bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
example_gt <- gen_tibble(bedfile, backingfile = tempfile("example_bed"))
```

Let us explore our table:
  ```{r}
example_gt
```

Let's get an overview of the populations:
```{r}
example_gt %>% group_by(population) %>% tally()
```
The loci had already been QC, but not cleaned for LD. We will use a PCA method, ''autoSVD',
that is designed to deal with LD directly via clumping:
  ```{r}
example_pca <- example_gt %>% gt_pca_autoSVD()
```

We can get a quick summary of the resulting object:
  ```{r}
example_pca
```
Objects outputted by different packages and analyes can be difficult to manipulate.
The package `broom` from the `tidyverse` offers methods
to easily extract information in form of tibbles (`tidy`) and merge some of those outputs with the
original data (`augment`).
`tidypopgen` implements such methods for objects created from its functions. So,
we can tidy the results as we would do for a standard PCA object from `prcomp`.
For example, we can get the Eigen values with:
  ```{r}
tidy(example_pca, matrix="eigenvalues")
```

And the scores for each individual:
  ```{r}
tidy(example_pca)
```

We can create quick plots with `autoplot`:
  ```{r}
autoplot(example_pca)
```

And look at the scores with:
  ```{r}
autoplot(example_pca, type = "scores")

```
`autoplots` are delibeartely kept simple, they are just a way to quickly inspect the results.
They are `ggplot2` objects, and so they can be further improved with the usual
`ggplot2` grammar:
  ```{r}
library(ggplot2)
autoplot(example_pca, type = "scores") +
  aes(color = example_gt$population) +
  labs(color = "population")

```

For more advanced/customised plots, we generally want to bring the information into the
`gen_tibble`, so that we can create the plots from scratch. `augment` allows us
to do just that, without any complex wrangling:

  ```{r}
example_gt <- augment(example_pca , data= example_gt, k=2)
example_gt
```

We can then easily plot our data with:
  ```{r}
library(ggplot2)
example_gt %>% ggplot(aes(.fittedPC1, .fittedPC2, color = population)) +
  geom_point(size = 1.5)
```
We can see that the population do separate nicely into 3 main groups on the PCA.

Note that a number of analysis, such as PCA, do not allow for missing data. Let
us consider such a dataset:
  ```{r}
bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))
missing_gt
```

If we attempt to run a PCA on this dataset, we get:
  ```{r, error = TRUE}
missing_pca <- missing_gt %>% gt_pca_autoSVD()
```

It is possible to obviate to this problem by filtering loci with missing data,
but that might lose a lot of loci. The alternative is to inpute the missing
the data. `tidypopgen` provides wrappers for two fast imputation approaches
available in `bigsnpr`, a simple imputation (`gt_impute_simple`) based on the frequency of the alleles
at each locus (by random sampling, or used the mean or mode), and more sophisticated
approach (`gt_impute_xgboost`) that uses boosted trees to try and predict the most likely genotype.
These methods are fine to impute a few missing genotypes, but they should not be used for any
sophisticated imputation (e.g. of low coverage genomes).

We use the simple approach to fix our dataset:
  ```{r}
missing_gt <- gt_impute_simple(missing_gt)
```
We can now run the analysis (but it fails because we have too few SNPs!!!!):
  ```{r, error=TRUE}
missing_pca <- missing_gt %>% gt_pca_autoSVD()
```

We can use a standard PCA:
  ```{r}
missing_pca <- missing_gt %>% gt_pca_partialSVD()
missing_pca
```
Note that the observed genotypes have been kept, and we can go back to using
them (and thus have missing values) with:
  ```{r}
missing_gt <- missing_gt  %>% gt_set_imputed(FALSE)
```

