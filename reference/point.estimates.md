# Point estimates of the regression coefficients

Computes point estimates of the regression coefficients for an object of
class `lmBayes` or `'spmBayes'`.

## Usage

``` r
point.estimates(object, type = "sn", retain = "mean")
```

## Arguments

- object:

  An object of class `'lmBayes'` or `'spmBayes'`.

- type:

  A character string denoting which algorithm should be used to recover
  the point estimates. The options are: (i) `"sn"` for the scaled
  neighborhood criterion from Li and Lin (2010), (ii) `"cred.int"` which
  excludes a predictor if the 50% credible interval contains zero, (iii)
  `"post.mode"` for the posterior mode, or (iv) `"post.mean"` for the
  posterior mean. Default is `"sn"`.

- retain:

  A character string denoting which posterior summary should be retained
  if `type = "sn"` or `type = "cred.int"`. The options are: (i) `"mode"`
  to retain the posterior mode or (ii) `"mean"` to retain the posterior
  mean.

## Value

A numeric vector containing the point estimates of the regression
coefficients.

## References

Q. Li, and N. Lin (2010), The Bayesian elastic net. *Bayesian Analysis*,
5(1):151-170.

## Author

Santiago Marin
