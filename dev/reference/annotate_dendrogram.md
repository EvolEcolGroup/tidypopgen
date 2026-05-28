# Attach a dendrogram panel to a heatmap

Designed for use with \|\>. Accepts a ggplot or \<heatmap_dendro\> and
returns a \<heatmap_dendro\> with the dendrogram panel inserted directly
adjacent to the heatmap panel — pixel-perfect alignment regardless of
axis label width. Can be chained to add dendrograms on both sides.

## Usage

``` r
annotate_dendrogram(
  plot,
  hc,
  side = c("left", "top"),
  rel_size = 0.15,
  line_color = "grey30",
  line_size = 0.5
)
```

## Arguments

- plot:

  A ggplot (from heatmap_pairwise()) or \<heatmap_dendro\> from a prior
  annotate_dendrogram() call.

- hc:

  hclust object (from order functions that return hclust).

- side:

  "left" (row dendrogram) or "top" (column dendrogram).

- rel_size:

  Size of the dendrogram relative to the panel (null unit).

- line_color:

  Segment colour.

- line_size:

  Segment linewidth.

## Value

A \<heatmap_dendro\> object.
