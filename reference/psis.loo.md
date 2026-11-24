# PSIS-LOO

This function computes the Pareto-smoothed importance sampling
leave-one-out information criterion (PSIS-LOO) for an object of class
`lmBayes` or `spmBayes`. The PSIS-LOO is converted into deviance scale
so that it is comparable with other information criteria like AIC and
BIC.

## Usage

``` r
psis.loo(object, comp.r_eff = TRUE)
```

## Arguments

- object:

  An object of class `'lmBayes'` or `'spmBayes'`.

- comp.r_eff:

  Logical. Whether the function should estimate the relative effective
  sample size for the likelihood of each observation. If `FALSE`, all
  the relative effective sample sizes are set to one. Default is `TRUE`.

## Value

The estimated PSIS-LOO.

## Author

Santiago Marin
