
## <span style="font-weight: normal">**cpam** (**c**hange**p**oint **a**dditive **m**odels)</span>

<img src="assets/example_trends.png" width="800" height="600" />

**cpam** is an R package for omics time series analysis.

## Key features

- **Changepoint detection**: Identify sharp transitions in gene
  expression.
- **smooth trends**: Model gene expression as a smooth function of time.
- **Shape-constrained trends**: Cluster genes into biologically
  meaningful temporal shape classes.
- **Quantification uncertainty**: Account for uncertainty in expression
  estimates.
- **Transcript-level analysis**: Perform gene- or transcript-level
  inferences. Aggregate $p$-values at the gene level for improved power.
- **Case-only or case-control time series**: Analyse time series data
  with or without controls.
- **User-friendly**: Sensible defaults and an interactive *shiny*
  interface.

Our new package **cpam** provides a comprehensive framework for
analysing time series omics data that combines statistical rigour with
practical utility. The method leverages modern statistical approaches
while remaining user-friendly, through sensible defaults and an
interactive interface. Researchers can directly address key questions in
time series analysis—when changes occur, what patterns they follow, and
how responses are related. While we have focused on transcriptomics, the
framework is applicable to other high-dimensional time series
measurements, making it a valuable addition to the omics analysis
toolkit.

## Installation

``` r
# Installation code here
remotes::install_github("l-a-yates/cpam")
```

## Usage

Load package and create a tibble for the experimental design. Here we
are using kallisto-generated counts, so we provide the file path for
each sample.

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

You will also need a tibble with transcript-to-gene mapping,,unless your
counts are aggregated at the gene level (transcript-level analysis is
recommened). E.g., for *Arabidopsis thaliana*:

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

Now you can run **cpam**

``` r
  cpo <- prepare_cpam(exp_design = exp_design,
                      t2g = t2g,
                      import_type = "kallisto",
                      num_cores = 5) # about 1 min
  cpo <- compute_p_values(cpo) # 30 secs to 5 min
  cpo <- estimate_changepoint(cpo) # 1-10 min
  cpo <- select_shape(cpo) # 1-20 min
```
