# WAIC

This function computes the Watanabe–Akaike information criterion (WAIC)
for an object of class `lmBayes` or `spmBayes`. The WAIC is converted
into deviance scale so that it is comparable with other information
criteria like AIC and BIC.

## Usage

``` r
widely.aic(object)
```

## Arguments

- object:

  An object of class `'lmBayes'` or `spmBayes`.

## Value

The estimated WAIC.

## Author

Santiago Marin
