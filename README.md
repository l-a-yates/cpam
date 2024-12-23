cpam: Shape-constrained changepoint additive models for time series
omics data
================

- [Analysis of time series omics data with
  cpam](#analysis-of-time-series-omics-data-with-cpam)
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)

## Analysis of time series omics data with cpam

**cpam** (**c**hange**p**oint **a**dditive **m**odels) fits smooth
shape-constrained additive models to omics count data.

## Overview

Brief description (2-3 sentences) of what cpam does and its key
features.

## Installation

``` r
# Installation code here
remotes::install_github("l-a-yates/cpam")
```

## Usage

Load package and create a tibble for the experimental design.

``` r
library(cpam)

exp_design
#> # A tibble: 50 × 4
#>    sample  time path                                condition
#>    <chr>  <dbl> <chr>                               <chr>    
#>  1 JHSS01     0 output/kallisto/JHSS01/abundance.h5 treatment
#>  2 JHSS02     0 output/kallisto/JHSS02/abundance.h5 treatment
#>  3 JHSS03     0 output/kallisto/JHSS03/abundance.h5 treatment
#>  4 JHSS04     0 output/kallisto/JHSS04/abundance.h5 treatment
#>  5 JHSS05     0 output/kallisto/JHSS05/abundance.h5 treatment
#>  6 JHSS06     5 output/kallisto/JHSS06/abundance.h5 treatment
#>  7 JHSS07     5 output/kallisto/JHSS07/abundance.h5 treatment
#>  8 JHSS08     5 output/kallisto/JHSS08/abundance.h5 treatment
#>  9 JHSS09     5 output/kallisto/JHSS09/abundance.h5 treatment
#> 10 JHSS10     5 output/kallisto/JHSS10/abundance.h5 treatment
#> # ℹ 40 more rows
```

You will need a tibble with transcript-to-gene mapping. E.g.,

``` r
t2g
#> # A tibble: 54,013 × 2
#>    target_id   gene_id  
#>    <chr>       <chr>    
#>  1 AT1G01010.1 AT1G01010
#>  2 AT1G01020.2 AT1G01020
#>  3 AT1G01020.6 AT1G01020
#>  4 AT1G01020.1 AT1G01020
#>  5 AT1G01020.4 AT1G01020
#>  6 AT1G01020.5 AT1G01020
#>  7 AT1G01020.3 AT1G01020
#>  8 AT1G03987.1 AT1G03987
#>  9 AT1G01030.2 AT1G01030
#> 10 AT1G01030.1 AT1G01030
#> # ℹ 54,003 more rows
```
