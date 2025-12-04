# Removes lowly expressed genes

Removes lowly expressed genes

## Usage

``` r
ts_filter(data, min_reads = 5, min_prop = 3/5)
```

## Arguments

- data:

  A tibble or data.frame containing columns:

  - target_id (character): Transcript identifiers

  - time (numeric): Time point of measurement

  - counts (numeric): Read counts

- min_reads:

  minimum reads per transcript per sample

- min_prop:

  minimum proportion of samples that exceed `min_read` at a given time
  point (default: 3/5)

## Value

a character vector of transcript IDs to keep

## Details

Identifies targets that show strong and consistent expression in at
least one timepoint. For each timepoint, the function calculates the
proportion of samples where a targets exceeds `min_reads`. Targets are
retained if they meet the minimum proportion (`min_prop`) at any
timepoint in the experiment.

## Examples

``` r
data <- dplyr::tibble(
  target_id = rep(paste0("t", 1:3), each = 6),
  time = rep(c(0, 4, 8), 6),
  counts = c(6,6,6, 0,0,0, 6,0,6, 0,6,6, 6,6,6, 6,0,0)
)
ts_filter(data)
#> [1] "t2" "t3"
```
