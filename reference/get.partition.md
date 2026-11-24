# Recover a partition

This function recovers the partition of the regression coefficients
induced by the nonparametric Bayesian Lasso. This function is based on
the implementation from the `"BNPmix"` package.

## Usage

``` r
get.partition(S.mcmc, loss = "VI")
```

## Arguments

- S.mcmc:

  A matrix of size `n.draws`-by-`n.preds`, where each row indicates to
  which cluster the regression coefficients belong to.

- loss:

  A loss function defined on the space of partitions. It can be either
  the variation of information loss function (`"VI"`) or the Binder loss
  function (`"Binder"`). Default is `"VI"`. See Wade and
  Ghahramani (2018) for additional details.

## References

S. Wade and Z. Ghahramani (2018), Bayesian cluster analysis: Point
estimation and credible balls (with discussion). *Bayesian Analysis*,
13(2):559-626.

## Author

Santiago Marin
