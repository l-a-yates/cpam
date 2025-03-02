
## cpam <span style="font-weight: normal">(**c**hange**p**oint **a**dditive **m**odels)</span>

<!-- badges: start -->

[![R-CMD-check](https://github.com/l-a-yates/cpam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/l-a-yates/cpam/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

### An R package for omics time series analysis.

Read the methods paper [here](https://doi.org/10.1101/2024.12.22.630003)

<img src="man/figures/example_trends.png" width="800" height="600" />

Application of **cpam** to RNA-seq time series of *Arabidopsis* plants
treated with excess-light.

## Key features

- **Changepoint detection**: Identify sharp transitions in expression.
- **smooth trends**: Model expression as a smooth function of time.
- **Shape-constrained trends**: Cluster targets into biologically
  meaningful temporal shape classes.
- **Quantification uncertainty**: Account for uncertainty in expression
  estimates.
- **Transcript-level analysis**
  - Perform gene- or transcript-level inferences.
  - Aggregate $p$-values at the gene level for improved power.
- **Case-only or case-control time series**: Analyse time series data
  with or without controls.
- **User-friendly**: Sensible defaults and an interactive *shiny*
  interface.

Our new package **cpam** provides a comprehensive framework for
analysing time series omics data. The method uses modern statistical
approaches while remaining user-friendly, through sensible defaults and
an interactive interface. Researchers can directly address key questions
in time series analysis—when changes occur, what patterns they follow,
and how responses are related. While we have focused on transcriptomics,
the framework is applicable to other high-dimensional time series
measurements.

If you encounter issues or have suggestions for improvements, please
open an [issue](https://github.com/l-a-yates/cpam/issues). We welcome
questions and discussion about using **cpam** for your research through
[Discussions](https://github.com/l-a-yates/cpam/discussions/). Our goal
is to work with users to make **cpam** a robust and valuable tool for
time series omics analysis. We can also be contacted via the email
addresses listed in our paper
[here](https://doi.org/10.1101/2024.12.22.630003).

## Installation

The package is available on CRAN and can be installed using the
following command:

``` r
install.packages("cpam")
```

## Usage

### Step 1: Load the package

``` r
library(cpam)
```

### Step 2: Create a tibble for the experimental design.

In this *Arabidopsis thaliana* time series example, we used the software
[kallisto](https://doi.org/10.1038/nbt.3519) to generate counts from
RNA-seq data. To load the counts, we provide the file path for each
kallisto output file (alternatively you can provide the counts directly
as count matrix, or use other quantification software)

``` r
exp_design_path
#>    sample time                                path condition
#> 1  JHSS01    0 output/kallisto/JHSS01/abundance.h5 treatment
#> 2  JHSS02    0 output/kallisto/JHSS02/abundance.h5 treatment
#> 3  JHSS03    0 output/kallisto/JHSS03/abundance.h5 treatment
#> 4  JHSS04    0 output/kallisto/JHSS04/abundance.h5 treatment
#> 5  JHSS05    0 output/kallisto/JHSS05/abundance.h5 treatment
#> 6  JHSS06    5 output/kallisto/JHSS06/abundance.h5 treatment
#> 7  JHSS07    5 output/kallisto/JHSS07/abundance.h5 treatment
#> 8  JHSS08    5 output/kallisto/JHSS08/abundance.h5 treatment
#> 9  JHSS09    5 output/kallisto/JHSS09/abundance.h5 treatment
#> 10 JHSS10    5 output/kallisto/JHSS10/abundance.h5 treatment
#> 11 JHSS16   10 output/kallisto/JHSS16/abundance.h5 treatment
#> 12 JHSS17   10 output/kallisto/JHSS17/abundance.h5 treatment
#> 13 JHSS18   10 output/kallisto/JHSS18/abundance.h5 treatment
#> 14 JHSS19   10 output/kallisto/JHSS19/abundance.h5 treatment
#> 15 JHSS20   10 output/kallisto/JHSS20/abundance.h5 treatment
#> 16 JHSS26   20 output/kallisto/JHSS26/abundance.h5 treatment
#> 17 JHSS27   20 output/kallisto/JHSS27/abundance.h5 treatment
#> 18 JHSS28   20 output/kallisto/JHSS28/abundance.h5 treatment
#> 19 JHSS29   20 output/kallisto/JHSS29/abundance.h5 treatment
#> 20 JHSS30   20 output/kallisto/JHSS30/abundance.h5 treatment
#> 21 JHSS36   30 output/kallisto/JHSS36/abundance.h5 treatment
#> 22 JHSS37   30 output/kallisto/JHSS37/abundance.h5 treatment
#> 23 JHSS38   30 output/kallisto/JHSS38/abundance.h5 treatment
#> 24 JHSS39   30 output/kallisto/JHSS39/abundance.h5 treatment
#> 25 JHSS40   30 output/kallisto/JHSS40/abundance.h5 treatment
#> 26 JHSS46   45 output/kallisto/JHSS46/abundance.h5 treatment
#> 27 JHSS47   45 output/kallisto/JHSS47/abundance.h5 treatment
#> 28 JHSS48   45 output/kallisto/JHSS48/abundance.h5 treatment
#> 29 JHSS49   45 output/kallisto/JHSS49/abundance.h5 treatment
#> 30 JHSS50   45 output/kallisto/JHSS50/abundance.h5 treatment
#> 31 JHSS56   60 output/kallisto/JHSS56/abundance.h5 treatment
#> 32 JHSS57   60 output/kallisto/JHSS57/abundance.h5 treatment
#> 33 JHSS58   60 output/kallisto/JHSS58/abundance.h5 treatment
#> 34 JHSS59   60 output/kallisto/JHSS59/abundance.h5 treatment
#> 35 JHSS60   60 output/kallisto/JHSS60/abundance.h5 treatment
#> 36 JHSS66   90 output/kallisto/JHSS66/abundance.h5 treatment
#> 37 JHSS67   90 output/kallisto/JHSS67/abundance.h5 treatment
#> 38 JHSS68   90 output/kallisto/JHSS68/abundance.h5 treatment
#> 39 JHSS69   90 output/kallisto/JHSS69/abundance.h5 treatment
#> 40 JHSS70   90 output/kallisto/JHSS70/abundance.h5 treatment
#> 41 JHSS76  180 output/kallisto/JHSS76/abundance.h5 treatment
#> 42 JHSS77  180 output/kallisto/JHSS77/abundance.h5 treatment
#> 43 JHSS78  180 output/kallisto/JHSS78/abundance.h5 treatment
#> 44 JHSS79  180 output/kallisto/JHSS79/abundance.h5 treatment
#> 45 JHSS80  180 output/kallisto/JHSS80/abundance.h5 treatment
#> 46 JHSS86  240 output/kallisto/JHSS86/abundance.h5 treatment
#> 47 JHSS87  240 output/kallisto/JHSS87/abundance.h5 treatment
#> 48 JHSS88  240 output/kallisto/JHSS88/abundance.h5 treatment
#> 49 JHSS89  240 output/kallisto/JHSS89/abundance.h5 treatment
#> 50 JHSS90  240 output/kallisto/JHSS90/abundance.h5 treatment
```

### Step 3: Obtain a table with the transcript-to-gene mapping

N.B. This is not needed if your counts are aggregated at the gene level,
but transcript-level analysis with aggregation of $p$-values to the gene
level is recommended. E.g., for *Arabidopsis thaliana*:

``` r
t2g_arabidopsis
#>      target_id   gene_id
#> 1  AT1G01010.1 AT1G01010
#> 2  AT1G01020.2 AT1G01020
#> 3  AT1G01020.6 AT1G01020
#> 4  AT1G01020.1 AT1G01020
#> 5  AT1G01020.4 AT1G01020
#> 6  AT1G01020.5 AT1G01020
#> 7  AT1G01020.3 AT1G01020
#> 8  AT1G03987.1 AT1G03987
#> 9  AT1G01030.2 AT1G01030
#> 10 AT1G01030.1 AT1G01030
#> 11 AT1G01040.1 AT1G01040
#> 12 AT1G01040.2 AT1G01040
#> 13 AT1G03993.1 AT1G03993
#> 14   at1g01046 AT1G01046
#> 15 AT1G01050.1 AT1G01050
#> 16 AT1G01050.2 AT1G01050
#> 17 AT1G03997.1 AT1G03997
#> 18 AT1G01060.6 AT1G01060
#> 19 AT1G01060.3 AT1G01060
#> 20 AT1G01060.1 AT1G01060
#> 21 AT1G01060.7 AT1G01060
#> 22 AT1G01060.2 AT1G01060
#> 23 AT1G01060.4 AT1G01060
#> 24 AT1G01060.8 AT1G01060
#> 25 AT1G01060.5 AT1G01060
#> 26 AT1G01070.1 AT1G01070
#> 27 AT1G01070.2 AT1G01070
#> 28 AT1G04003.1 AT1G04003
#> 29 AT1G01080.1 AT1G01080
#> 30 AT1G01080.3 AT1G01080
#> 31 AT1G01080.2 AT1G01080
#> 32 AT1G01090.1 AT1G01090
#> 33 AT1G01100.2 AT1G01100
#> 34 AT1G01100.1 AT1G01100
#> 35 AT1G01100.4 AT1G01100
#> 36 AT1G01100.3 AT1G01100
#> 37 AT1G01110.2 AT1G01110
#> 38 AT1G01110.1 AT1G01110
#> 39 AT1G01120.1 AT1G01120
#> 40 AT1G01130.1 AT1G01130
#> 41 AT1G01140.1 AT1G01140
#> 42 AT1G01140.2 AT1G01140
#> 43 AT1G01140.3 AT1G01140
#> 44 AT1G01150.1 AT1G01150
#> 45 AT1G01150.2 AT1G01150
#> 46 AT1G01150.3 AT1G01150
#> 47 AT1G01160.1 AT1G01160
#> 48 AT1G01160.2 AT1G01160
#> 49 AT1G04007.1 AT1G04007
#> 50 AT1G01170.2 AT1G01170
```

### Step 4: Run **cpam**

``` r
  cpo <- prepare_cpam(exp_design = exp_design_path,
                      count_matrix = NULL,
                      t2g = t2g_arabidopsis,
                      model = "case-only",
                      import_type = "kallisto",
                      num_cores = 5)
  cpo <- compute_p_values(cpo) 
  cpo <- estimate_changepoint(cpo) 
  cpo <- select_shape(cpo) 
```

### Step 5: Visualise the results

Load the shiny app for an interactive visualisation of the results:

``` r
  visualise(cpo) # not shown in vignette
```

Or plot one gene at a time:

``` r
  plot_cpam(cpo, gene_id = "AT3G23280")
```

<img src="man/figures/example_gene_plot.png" width="650" height="550" />

Isoform 1 (AT3G23280.1) has a changepoint at 67.5 min and has a
monotonic increasing concave (micv) shape. Isoform 2 (AT3G23280.2) has
no changepoint and has an unconstrained thin-plate (tp) shape. <br><br>

We can generate a results table which has $p$-values, shapes, log-fold
changes and counts with many optimal filters (see tutorials):

``` r
  results(cpo)
#> # A tibble: 15,279 × 25
#>    target_id   gene_id     p    cp shape lfc.0 lfc.5 lfc.10 lfc.20 lfc.30 lfc.45
#>    <chr>       <chr>   <dbl> <dbl> <chr> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#>  1 AT1G01910.1 AT1G01…     0     0 micv      0 1.01   1.70   2.38   2.60   2.73 
#>  2 AT1G01910.2 AT1G01…     0    10 cv        0 0      0      0.553  0.775  0.790
#>  3 AT1G01910.5 AT1G01…     0    10 cx        0 0      0     -3.20  -4.57  -4.82 
#>  4 AT1G02610.1 AT1G02…     0    45 mdcx      0 0      0      0      0      0    
#>  5 AT1G02610.2 AT1G02…     0    10 cx        0 0      0     -0.645 -1.16  -1.71 
#>  6 AT1G02610.3 AT1G02…     0    10 mdcx      0 0      0     -1.48  -2.11  -2.25 
#>  7 AT1G04080.1 AT1G04…     0    10 cv        0 0      0      2.75   3.85   3.97 
#>  8 AT1G04080.2 AT1G04…     0    45 micv      0 0      0      0      0      0    
#>  9 AT1G04080.3 AT1G04…     0     0 micv      0 0.268  0.445  0.603  0.638  0.656
#> 10 AT1G04080.5 AT1G04…     0    10 cx        0 0      0     -2.17  -3.04  -3.10 
#> # ℹ 15,269 more rows
#> # ℹ 14 more variables: lfc.60 <dbl>, lfc.90 <dbl>, lfc.180 <dbl>,
#> #   lfc.240 <dbl>, counts.0 <dbl>, counts.5 <dbl>, counts.10 <dbl>,
#> #   counts.20 <dbl>, counts.30 <dbl>, counts.45 <dbl>, counts.60 <dbl>,
#> #   counts.90 <dbl>, counts.180 <dbl>, counts.240 <dbl>
```

### Tutorials

For a quick-to-run introductory example, we have provided a small
simulated data set as part of the package.

- [Introductory
  Example](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/example.html)

The following two tutorials use real-world data to demonstrate the
capabilities of the **cpam** package. In addition, they provide code to
reproduce the results for the case studies presented in the
[manuscript](https://doi.org/10.1101/2024.12.22.630003) accompanying the
**cpam** package.

- [Arabidopsis Case
  Study](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/crisp.html)
- [Human Embryo Case
  Study](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/torre.html)

## Acknowledgements

This work was supported by the Australian Research Council, Centre of
Excellence for Plant Success in Nature and Agriculture (CE200100015).

<img src="man/figures/plant_success_logo_horizontal.jpg" width="90%" height="auto" alt="Plant Success Logo">
