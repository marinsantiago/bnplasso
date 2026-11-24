# Posterior co-clustering probabilities plot

This function plots the posterior co-clustering probabilities of the
regression coefficients

## Usage

``` r
coclust.probs(
  S.mcmc,
  return.mat = FALSE,
  viridis.pal = "C",
  main = NULL,
  axis.idx = NULL,
  xy.labs = NULL
)
```

## Arguments

- S.mcmc:

  A matrix of size `n.draws`-by-`n.preds`, where each row indicates to
  which cluster the regression coefficients belong to.

- return.mat:

  logical. Whether the matrix of posterior co-clustering probabilities
  should be returned. Default is `FALSE`

- viridis.pal:

  Character describing which "viridis" palette should be employed.
  Options are "A", "B", "C", "D", "E", "F", "G", and "H".

- main:

  Title of the plot.

- axis.idx:

  Sequence of number to appear in the axis of the plot.

- xy.labs:

  A title for the x and y axes.

## Author

Santiago Marin
