# Posterior predictive distribution for new data

Compute the posterior predictive distribution for an object of class
`lmBayes`.

## Usage

``` r
# S3 method for class 'lmBayes'
predict(object, X.new, ...)
```

## Arguments

- object:

  An object of class `'lmBayes'`.

- X.new:

  A new matrix of predictors, where each row is a new observation.

- ...:

  Further arguments passed to.

## Value

A matrix where each row corresponds to an MCMC draw and each column to
an observation in the new data.

## Author

Santiago Marin
