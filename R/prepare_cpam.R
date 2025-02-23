#' Prepare a cpam object
#'
#' @param exp_design a dataframe or tibble with the experimental design, containing at least a 'time' and a 'sample' column
#' @param count_matrix a matrix of counts. Column names must be in 'sample' column of `exp_design`,
#' @param t2g a transcript to gene dataframe or tibble with columns target_id and gene_id
#' @param import_type software used for quantification, one of "kallisto", "salmon" ,.... Ignored if `count_matrix` is supplied.
#' @param model_type "case-only" (default) or "case-control"
#' @param bootstrap logical; load bootstrap samples, also called inferential replicates, if available, and rescale counts.
#' @param filter_fun filter function to remove lowly expressed genes (default is `filter_fun()`)
#' @param filter_fun_args arguments for filter function
#' @param regularize logical; use empirical Bayes regularization of dispersions (default is TRUE)
#' @param gene_level logical; aggregate counts to gene level before data preparation and modelling (default is FALSE)
#' @param aggregate_to_gene logical; aggregate p values from transcript- to gene-level
#' @param condition_var string; column name in `exp_design` for the condition variable (for `model_type` = "case_control" only)
#' @param case_value value of `condition_var` that indicates the "case". All other values are deemed to be control
#' @param num_cores integer; number of cores to use for parallel computation
#' @param normalize logical; use model offsets based on sampling depth and gene length
#' @param fixed_effects a model formula of the form `~ effect1 + effect2`
#' @param intercept_cc string; intercept for case-control model: "1" (default) for common intercept  or "condition"
#'
#' @details
#' This function prepares a cpam object for analysis.
#'  The function loads count data from files or a matrix,
#'  filters lowly expressed genes, computes normalisation factors,
#'  and estimates dispersions. Many of these steps can be customised or turned off.
#'
#'  When bootstrap samples (inferential replicates) are available, it loads and
#'  summarises these using means, standard errors, and estimated overdispersions.
#'  The latter are a measure of quantification uncertainty and they are used to
#'  rescale the counts which accouts for this uncertainty during the modelling steps.
#'
#'  The data within the cpam object are accessible via the slots.
#'
#' @return an object of class cpam
#' @export
#'
#' @examples
#' \dontrun{
#'
#' library(cpam)
#' library(dplyr)
#'
#' # Example Experimental Design
#' exp_design <- tibble(sample = paste0("s",1:50),
#'                      time = rep(c(0:4,10),
#'                      path = paste0("path/",sample,"/abundance.h5"))
#'
#' # Example Transcript-to-Gene Mapping
#' t2g <- readr::read_csv("path/to/t2g.csv")
#'
#' # Prepare a cpam object
#' cpo <- prepare_cpam(
#'  exp_design = exp_design,
#'  t2g = t2g,
#'  import_type = "kallisto",
#'  num_cores = 5)
#'
#' # Print the object
#' cpo
#' }
#'
#' @references
#'  Pedro L Baldoni, Yunshun Chen, Soroor Hediyeh-zadeh, Yang Liao, Xueyi Dong,
#'  Matthew E Ritchie, Wei Shi, Gordon K Smyth,
#'  Dividing out quantification uncertainty allows efficient assessment of
#'  differential transcript expression with edgeR,
#'  Nucleic Acids Research, Volume 52, Issue 3, 9 February 2024,
#'  Page e13, https://doi.org/10.1093/nar/gkad1167
#'
#'  Yunshun Chen, Lizhong Chen, Aaron T L Lun, Pedro L Baldoni, Gordon K Smyth,
#'  edgeR v4: powerful differential analysis of sequencing data with
#'  expanded functionality and improved support for small counts
#'  and larger datasets, Nucleic Acids Research, Volume 53, Issue 2,
#'  27 January 2025, https://doi.org/10.1093/nar/gkaf018
#'
#'
prepare_cpam <- function(exp_design,
                         count_matrix = NULL,
                         t2g = NULL,
                         import_type = NULL,
                         model_type = c("case-only","case-control"),
                         bootstrap = TRUE, # use bootstraps if available
                         filter_fun = "ts_filter",
                         filter_fun_args = list(min_reads = 5, min_prop = 3/5),
                         regularize = TRUE,
                         gene_level = FALSE,
                         aggregate_to_gene = !gene_level,
                         condition_var = "condition",
                         case_value = "treatment",
                         num_cores = 1,
                         normalize = TRUE,
                         fixed_effects = NULL,
                         intercept_cc = c("1",condition_var)
){

  model_type <- match.arg(model_type)
  intercept_cc <- match.arg(intercept_cc)
  import <- ifelse(is.null(count_matrix),T,F)

  validate_inputs(exp_design, count_matrix, t2g, import_type, model_type,
                  condition_var, case_value, fixed_effects,
                  aggregate_to_gene, gene_level, import)

  if (model_type == "case-control") {
    exp_design <- exp_design %>%
      dplyr::mutate(case = as.numeric(.data[[condition_var]] == case_value))
  }

  fixed_effects <- validate_fixed_effects(fixed_effects,exp_design)

  num_cores <- validate_cores(num_cores)

  # Import data if needed, otherwise use provided count matrix
  if (import) {
    message(paste0("Loading ", length(exp_design$path), " samples"))
    result <- import_count_data(exp_design, t2g, import_type, gene_level, bootstrap)
    counts_raw <- result$counts_raw
    bootstrap <- result$bootstrap
    overdispersion.prior <- result$overdispersion.prior
    nboot <- result$nboot
    catch <- result$catch
    boot <- result$boot
  } else {
    counts_raw <- process_count_matrix(count_matrix, t2g, gene_level)
    overdispersion.prior <- nboot <- NULL
    bootstrap <- F
  }

  # Convert to long format and join with experiment design
  data_long <- convert_to_long_format(counts_raw, exp_design)

  # Handle gene aggregation
  data_long <- handle_gene_aggregation(data_long, t2g, aggregate_to_gene)

  # Handle bootstrap if needed
  if (bootstrap) {
    data_long <- process_bootstrap_data(data_long, catch, boot)
  }

  # Filter low count genes
  message("Filtering low count genes")
  target_to_keep <- do.call(filter_fun, args = c(list(data = data_long), filter_fun_args))
  count_matrix_filtered <- counts_raw[target_to_keep, ]

  # Normalize if requested
  norm_factor <- compute_normalization_factors(count_matrix_filtered, normalize)

  # Filter and update data_long
  data_long <-   data_long %>%
    dplyr::filter(.data$target_id %in% target_to_keep) %>%
    dplyr::mutate(norm_factor = norm_factor[sample]) %>%
    dplyr::relocate("target_id")

  # Estimate dispersions if regularization is requested
  if (regularize) {
    dispersions <- estimate_dispersions_wrapper(count_matrix_filtered, exp_design, bootstrap, catch)
    data_long <- data_long %>% dplyr::mutate(disp = dispersions[.data$target_id] %>% as.numeric)
  }

  # create cpam object
  list(exp_design = exp_design %>% dplyr::select(-tidyr::any_of(c("path"))),
       count_matrix_raw = counts_raw,
       count_matrix_filtered = count_matrix_filtered,
       target_to_keep = target_to_keep,
       data_long = data_long,
       t2g = t2g,
       regularize = regularize,
       overdispersion.prior = overdispersion.prior,
       model_type = model_type,
       condition_var = condition_var,
       case_value = case_value,
       bootstrap = bootstrap,
       nboot = nboot,
       filter = list(filter_fun = filter_fun, filter_fun_args = filter_fun_args),
       gene_level = gene_level,
       aggregate_to_gene = aggregate_to_gene,
       times = sort(unique(exp_design$time)),
       num_cores = num_cores,
       fixed_effects = fixed_effects,
       intercept_cc = intercept_cc,
       bss = c("micv","mdcx","cv","cx","micx","mdcv","tp")
  ) %>%
    `class<-`("cpam")
}


#' Removes lowly expressed genes
#'
#' @param data A tibble or data.frame containing columns:
#'   \itemize{
#'     \item target_id (character): Transcript identifiers
#'     \item time (numeric): Time point of measurement
#'     \item counts (numeric): Read counts
#'   }
#' @param min_reads minimum reads per transcript per sample
#' @param min_prop minimum proportion of samples that exceed `min_read` at a
#' given time point (default: 3/5)
#'
#' @details
#' Identifies targets that show strong and consistent expression in at least one timepoint.
#' For each timepoint, the function calculates the proportion of samples
#' where a targets exceeds `min_reads`. Targets are retained if they meet
#' the minimum proportion (`min_prop`) at any timepoint in the experiment.
#'
#' @return a character vector of transcript IDs to keep
#' @export
#'
#' @examples
#' data <- dplyr::tibble(
#'   target_id = rep(paste0("t", 1:3), each = 6),
#'   time = rep(c(0, 4, 8), 6),
#'   counts = c(6,6,6, 0,0,0, 6,0,6, 0,6,0, 6,6,6, 0,0,0)
#' )
#' ts_filter(data)
#'
ts_filter <- function(data, min_reads = 5, min_prop = 3/5) {
  data %>%
    dplyr::summarise(k = mean(.data$counts >= min_reads) >= min_prop, .by = c("target_id","time")) %>%
    dplyr::summarise(keep = any(.data$k), .by = "target_id") %>%
    dplyr::filter(.data$keep) %>%
    dplyr::pull("target_id")
}

#' Summarize bootstrap samples
#'
#' Internal function to process inferential replicates from tximport objects.
#'
#' @param txi A tximport object containing inferential replicates
#' @return A tibble with bootstrap summary statistics
summarise_bootstraps <- function(txi){
  message("Summarising bootstrap samples")
  txi$infReps %>%
    purrr::set_names(colnames(txi$counts)) %>%
    purrr::map(`rownames<-`, rownames(txi$counts)) %>%
    purrr::imap(~ matrixStats::rowQuantiles(.x,probs = c(pnorm(-1),pnorm(1))) %>%
                  `colnames<-`(c("q_lo","q_hi")) %>%
                  dplyr::as_tibble(rownames = "target_id")) %>%
    dplyr::bind_rows(.id = "sample")
}

#' Estimate dispersions for count data
#'
#' Internal function to calculate tagwise dispersions using edgeR.
#'
#' @param counts Count matrix (genes/transcripts × samples)
#' @param exp_design Experimental design with time column
#' @return Named vector of tagwise dispersion estimates
#' @noRd
estimate_dispersions <- function(counts, exp_design){
  # update for case-control series (i.e., use condition in the design matrix)
  message("Estimating dispersions using edgeR")
  if(length(unique(exp_design$time)) > 3){
    design <- stats::model.matrix(~poly(time,2), data = exp_design)
  } else design <- stats::model.matrix(~1, data = exp_design)

  disp <- edgeR::DGEList(counts,
                         sample = data.frame(time = exp_design$time)) %>%
    {edgeR::estimateDisp.default (.,design = design)} %>%
    {.$tagwise.dispersion}

  names(disp) <- rownames(counts)

  disp
}


#' Validate input parameters for prepare_cpam
#'
#' @param exp_design Experimental design data frame
#' @param count_matrix Count matrix
#' @param t2g Transcript to gene mapping
#' @param import_type Import type
#' @param model_type Model type
#' @param condition_var Condition variable
#' @param case_value Case value
#' @param fixed_effects Fixed effects formula
#' @param aggregate_to_gene Whether to aggregate to gene level
#' @param gene_level Whether data is at gene level
#' @param import Whether to import data
#'
#' @return NULL, throws error if validation fails
#' @keywords internal
validate_inputs <- function(exp_design, count_matrix, t2g, import_type, model_type,
                            condition_var, case_value, fixed_effects,
                            aggregate_to_gene, gene_level, import) {
  # Check time column exists with sufficient distinct values
  if (is.null(exp_design[["time"]])) {
    stop("'exp_design' must include a 'time' column")
  }
  if (length(unique(exp_design$time)) < 3) {
    stop("The experimental design must have at least three distinct time points")
  }

  # Check t2g requirement
  if (is.null(t2g) & aggregate_to_gene) {
    stop("'t2g' must be supplied if 'aggregate_to_gene' is true")
  }

  # Import-specific validations
  if (import) {
    if (is.null(import_type)) {
      stop("'import_type' must be given if count matrix is not supplied")
    }
    if (is.null(exp_design["path"])) {
      stop("'exp_design' must contain a 'path' column if count matrix is not supplied")
    }
    if (is.null(t2g) & !gene_level) {
      stop("'t2g' must be supplied if 'gene_level' is false")
    }
  } else {
    if (!is.null(import_type)) {
      message("'import_type' is being ignored since a count matrix has been supplied")
    }
    if (is.null(colnames(count_matrix)) | sum(!colnames(count_matrix) %in% exp_design$sample) > 0) {
      stop("Column names of 'count_matrix' must match samples in 'exp_design'")
    }
  }

  # Case-control specific validations
  if (model_type == "case-control") {
    if (is.null(exp_design[[condition_var]])) {
      stop(paste("type = 'case-control', but the column", condition_var, "does not exist in 'exp_design'"))
    }
    if (!is.null(exp_design[["case"]])) {
      stop("The column name 'case' is reserved and cannot be used in 'exp_design'")
    }
    if(all(exp_design[[condition_var]]!=case_value)){
      stop("type = 'case-control', but there are no case samples")
    }
    if(all(exp_design[[condition_var]]==case_value)){
      stop("type = 'case-control', but there are no control samples")
    }
    df <- exp_design %>% dplyr::select(tidyr::all_of(c("time",condition_var))) %>% dplyr::distinct()
    if(sum(exp_design[[condition_var]]==case_value) != sum(exp_design[[condition_var]]==case_value)){
      stop("The observed time points for case and control samples must be equal")
    }
  }
}


#' Validate and adjust number of cores
#'
#' @param num_cores Number of cores requested
#' @return Validated number of cores
#' @keywords internal
validate_cores <- function(num_cores) {
  if (!num_cores %% 1 == 0) {
    stop("num_cores must be integer values")
  }

  if(.Platform$OS.type == "windows"){
    warning("Parallel processing is not supported on Windows. Setting num_cores = 1")
    return(1)
  }

  available_cores <- parallel::detectCores()
  if (num_cores > available_cores) {
    suggested_cores <- floor(available_cores / 4)
    warning(paste0("num_cores is greater than ", available_cores,
                   ", the available number of cores. Setting num_cores = ",
                   suggested_cores, " (1/4 x number cores)"),
            call. = FALSE)
    num_cores <- suggested_cores
  }

  return(num_cores)
}


#' Import count data from files
#'
#' @param exp_design Experimental design
#' @param t2g Transcript to gene mapping
#' @param import_type Import type
#' @param gene_level Whether to aggregate to gene level
#' @param bootstrap Whether to use bootstrap samples
#' @return List containing imported data and related variables
#' @keywords internal
import_count_data <- function(exp_design, t2g, import_type, gene_level, bootstrap) {
  if (gene_level) {
    message("Summarizing to gene level")
    bootstrap <- FALSE
  }

  txi <- tximport::tximport(
    exp_design$path,
    txIn = TRUE,
    type = import_type,
    txOut = !gene_level,
    dropInfReps = !bootstrap,
    varReduce = FALSE,
    countsFromAbundance = "no",
    tx2gene = t2g
  )

  colnames(txi$counts) <- exp_design$sample
  counts_raw <- txi$counts

  # Handle bootstrap samples
  catch <- boot <- nboot <- overdispersion.prior <- NULL
  if (bootstrap && is.null(txi$infReps)) {
    message("No inferential replicates found, setting bootstrap = FALSE")
    bootstrap <- FALSE
  }

  if (bootstrap) {
    catch <- calculate_overdispersions(txi)
    overdispersion.prior <- catch$overdispersion.prior
    boot <- summarise_bootstraps(txi)
    nboot <- txi$infReps[[1]] %>% ncol
  }

  list(
    counts_raw = counts_raw,
    bootstrap = bootstrap,
    overdispersion.prior = overdispersion.prior,
    nboot = nboot,
    catch = catch,
    boot = boot
  )
}


#' Process count matrix when importing is not needed
#'
#' @param count_matrix Count matrix
#' @param t2g Transcript to gene mapping
#' @param gene_level Whether to aggregate to gene level
#' @return Processed count matrix
#' @keywords internal
process_count_matrix <- function(count_matrix, t2g, gene_level) {
  if (gene_level & !is.null(t2g)) {
    count_matrix %>%
      tidyr::as_tibble(rownames = "target_id") %>%
      dplyr::left_join(t2g, by = "target_id") %>%
      dplyr::group_by(.data[["gene_id"]]) %>%
      dplyr::summarise(dplyr::across(-.data[["target_id"]], sum)) %>%
      (\(x) {
        as.matrix(x[, -1]) %>% `rownames<-`(dplyr::pull(x, "gene_id"))
      })
  } else {
    count_matrix
  }
}

#' Convert count matrix to long format and join with experiment design
#'
#' @param counts_raw Raw count matrix
#' @param exp_design Experimental design
#' @return Long format data frame
#' @keywords internal
convert_to_long_format <- function(counts_raw, exp_design) {
  counts_raw %>%
    tidyr::as_tibble(rownames = "target_id") %>%
    tidyr::pivot_longer(-"target_id", names_to = "sample", values_to = "counts_raw") %>%
    dplyr::mutate(counts = counts_raw) %>%
    dplyr::left_join(exp_design %>% dplyr::select(-tidyr::any_of(c("path"))), by = "sample") %>%
    dplyr::arrange(.data$target_id)
}


#' Handle gene aggregation in long format data
#'
#' @param data_long Long format data
#' @param t2g Transcript to gene mapping
#' @param aggregate_to_gene Whether to aggregate to gene level
#' @return Updated long format data with gene ID information
#' @keywords internal
handle_gene_aggregation <- function(data_long, t2g, aggregate_to_gene) {
  if (aggregate_to_gene) {
    data_long %>% dplyr::left_join(t2g, by = c("target_id"))
  } else {
    data_long %>% dplyr::mutate(gene_id = .data$target_id)
  }
}

#' Process bootstrap data
#'
#' @param data_long Long format data
#' @param catch Overdispersion calculations
#' @param boot Bootstrap summary
#' @return Updated data with bootstrap information
#' @keywords internal
process_bootstrap_data <- function(data_long, catch, boot) {
  data_long %>%
    dplyr::mutate(
      overdispersions = catch$overdispersion[.data$target_id],
      counts_scaled = .data$counts_raw / .data$overdispersions,
      counts = .data$counts_scaled
    ) %>%
    dplyr::left_join(boot, by = c("sample", "target_id"))
}

#' Compute normalization factors
#'
#' @param count_matrix_filtered Filtered count matrix
#' @param normalize Whether to perform normalization
#' @return Vector of normalization factors
#' @keywords internal
compute_normalization_factors <- function(count_matrix_filtered, normalize) {
  if (normalize) {
    DESeq2::estimateSizeFactorsForMatrix(count_matrix_filtered)
  } else {
    rep(1, ncol(count_matrix_filtered)) %>%
      `names<-`(colnames(count_matrix_filtered))
  }
}

#' Wrapper for estimating dispersions
#'
#' @param count_matrix_filtered Filtered count matrix
#' @param exp_design Experimental design
#' @param bootstrap Whether bootstrap was used
#' @param catch Overdispersion calculations from bootstrap
#' @return Estimated dispersions
#' @keywords internal
estimate_dispersions_wrapper <- function(count_matrix_filtered, exp_design, bootstrap, catch) {
  if (bootstrap) {
    estimate_dispersions(
      count_matrix_filtered / catch$overdispersion[rownames(count_matrix_filtered)],
      exp_design
    )
  } else {
    estimate_dispersions(count_matrix_filtered, exp_design)
  }
}
