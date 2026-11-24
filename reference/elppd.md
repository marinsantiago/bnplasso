# Expected log pointwise predictive density (elppd)

This function computes the expected log pointwise predictive density
(elppd) for an object of class `lmBayes`. Given a test (held-out) set,
\\\\y\_{\text{test},i},\\
\mathbf{x}\_{\text{test},i}\\\_{i=1}^{n\_{\text{test}}}\\, the elppd is
defined as \$\$\text{elppd} =
\frac{1}{n\_{\text{test}}}\sum\_{i=1}^{n\_{\text{test}}}
\log\left(\frac{1}{S}\sum\_{s=1}^{S}p\left({y}\_{\text{test}, i} \|
\mathbf{x}\_{\text{test}, i}^{\_{'}}\boldsymbol{\beta}\_{(s)},
\sigma^{2}\_{(s)}\right)\right),\$\$ where \\p()\\ is the likelihood of
the observed data and \\\boldsymbol{\beta}\_{(s)}\\ and
\\\sigma^{2}\_{(s)}\\ are the \\s\\-th draws of the coefficient vector
and the sampling variance in the MCMC algorithm, respectively.

## Usage

``` r
elppd(object, X_new, y_new)
```

## Arguments

- object:

  An object of class `lmBayes`.

- X_new:

  A matrix of predictors from the test data, where each row is an
  observation.

- y_new:

  A vector of responses from the test data of size \\n\_{\text{test}}\\.

## Value

The elppd as a numeric scalar.

## References

A. Gelman, J. Carlin, H. Stern, D. Dunson, and A. Vehtari (2013),
*Bayesian Data Analysis*, 3. Chapman & Hall/CRC.

## Author

Santiago Marin
