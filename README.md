# bnplasso <img src="man/figures/bnplasso.png" alt="bnplasso" width="140" height="150" align="right"> 

<!-- badges: start -->

![R-CMD-check](https://github.com/marinsantiago/BOBgmms/workflows/R-CMD-check/badge.svg)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!-- badges: end -->

</br>

## Overview

This repository contains the R package `bnplasso` (developer's version), which implements the [Nonparametric Bayesian Lasso](https://arxiv.org/abs/2411.08262) (Marin et al., 2025+).

The main routine of the package, ```bnplasso.lm()```, returns an object of S3 class, ```"lmBayes"```, which is supported by various methods 
like ```print()```, ```summary()```, ```plot()```, ```fitted()```, ```residuals()```, ```coef()```, 
and ```predict()```, allowing users to quickly visualize, evaluate, and analyze the output in a familiar fashion. 
The package also implements the Bayesian Lasso (Park and Casella, 2008) and the Bayesian 
adaptive Lasso (Leng et al., 2014), through the functions ```blasso.lm()``` 
and ```balasso.lm()```, respectively.

## Installation

You can install the developer's version via `devtools` as:

``` r
# install.packages("devtools")
devtools::install_github("marinsantiago/bnplasso")
```

On the other hand, if you wish to install the package from the `bnplasso.zip` file in the supplementary materials to Marin et al. (2025+):

  1. Decompress the zip file `bnplasso.zip`. The folder `bnplasso` should result.
  
  2. In R, set your working directory to the folder `bnplasso`.
  
  3. Run the following R code:
  
``` r
# install.packages("devtools")
devtools::build()
devtools::install()
```

## Examples

Detailed guidelines for using the package functions are referred to their help pages in R. Additional examples are available at [https://github.com/marinsantiago/bnplasso-examples](https://github.com/marinsantiago/bnplasso-examples).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our paper:

```
@misc{marin_bnplasso,
  title         = {Adaptive Shrinkage with a Nonparametric Bayesian Lasso}, 
  author        = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  year          = {2024},
  eprint        = {2411.08262},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME}
}
```

## Disclaimer

```
# Disclaimer:
#
# The software is provided "as is", without warranty of any kind, express or implied,
# including but not limited to the warranties of merchantability, fitness for a particular
# purpose and noninfringement. In no event shall the authors or copyright holders be liable
# for any claim, damages, or other liability, whether in an action of contract, 
# tort or otherwise, arising from, out of, or in connection with the software or the use
# or other dealings in the software.
```


## <a name="refs"></a> References

Leng, C., Tran, M.-N., and Nott, D. (2014). "Bayesian Adaptive Lasso." *Annals of the Institute of Statistical Mathematics*, 66, 221–244

Marin, S., Loong, B., and Westveld, A. H. (2025+), "Adaptive Shrinkage with a Nonparametric Bayesian Lasso." *Journal of Computational and Graphical Statistics* (to appear).

Park, T. and Casella, G. (2008). "The Bayesian Lasso." *Journal of the American Statistical Association*, 103 (482), 681–686.