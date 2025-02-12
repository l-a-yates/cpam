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
#'  # Print the object
#'  cpo
#'  }
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

  # check that design matrix is non-singular
  if(is.null(exp_design[["time"]])) stop("'exp_design' must include a 'time' column")
  if(length(unique(exp_design$time))<3) stop("The experimental design must have a least three distinct time points")

  # make this a helper function check_fixed_effects()
  if(!is.null(fixed_effects)){
    e <- try(stats::model.frame(fixed_effects, exp_design), silent = T)
    if(inherits(e, "try-error")){
      stop((paste0("Fixed effects do not match the supplied design: ",
                   stringr::str_split_i(as.character(attr(e,"condition")),":",2))))
    }
    fixed_effects <- deparse(fixed_effects[[2]])
  }
  if(is.null(t2g) & aggregate_to_gene) stop("'t2g' must be supplied if 'aggregate_to_gene' is true")

  if(import){
    if(is.null(import_type)) stop("'import_type' must be given if count matrix is not supplied")
    if(is.null(exp_design["path"])) stop("'exp_design' must contain a 'path' column if count matrix is not supplied")
    if(is.null(t2g) & !gene_level) stop("'t2g' must be supplied if 'gene_level' is false")
  } else{
    if(!is.null(import_type)) message("'import_type' is being ignored since a count matrix has been supplied")
    #if(is.null(t2g) & gene_level) message("'t2g' has not been supplied. This is fine if the counts in the supplied count matrix are already aggregated to the gene level")
    if(is.null(colnames(count_matrix)) | sum(!colnames(count_matrix) %in% exp_design$sample)>0){
      stop("Column names of 'count_matrix' must match samples in 'exp_design'")}
    bootstrap <- F
  }


  if(model_type == "case-control"){
    if(is.null(exp_design[[condition_var]])){
      stop(paste("type = 'case-control', but the column",condition_var,"does not exist in 'exp_design"))
    }

    if(!is.null(exp_design[["case"]])) stop("The column name 'case' is reserved and cannot be used in 'exp_design'")
    exp_design <- exp_design %>% dplyr::mutate(case = as.numeric(exp_design[[condition_var]] == case_value))
    # check the case values by timepoint
  }

  if(!num_cores%%1==0) stop("num_cores must be integer values")
  if(num_cores > parallel::detectCores()){
    warning(paste0("num_cores is greater than ",parallel::detectCores(),
                   ", the available number of cores. Setting num_cores = ",
                   floor(parallel::detectCores()/4),
                   " (1/4 x number cores)"),
            call. = F)
    num_cores <- floor(parallel::detectCores()/4)}


  if(import){
    message(paste0("Loading ",length(exp_design$path), " samples"))
    if(gene_level){
      message("Summarizing to gene level")
      aggregate_to_gene <- F
      bootstrap <- F
    }
    txi <- tximport::tximport(
      exp_design$path,
      txIn = T,
      type = import_type,
      txOut = !gene_level,
      dropInfReps = !bootstrap,
      varReduce = F,
      countsFromAbundance = "no",
      tx2gene = t2g
    )

    colnames(txi$counts) <- exp_design$sample
    counts_raw <- txi$counts

    if(bootstrap & is.null(txi$infReps)){
      message("No inferential replicates found, setting bootstrap = FALSE")
      bootstrap <- F
    }

    if(bootstrap){
      catch <- calculate_overdispersions(txi)
      overdispersion.prior = catch$overdispersion.prior
      boot <- summarise_bootstraps(txi)
      nboot <- txi$infReps[[1]] %>% ncol
    } else {
      overdispersion.prior <- nboot <- NULL
    }
  } else {
    if(gene_level & !is.null(t2g)){
      counts_raw <-
        count_matrix %>%
        tidyr::as_tibble(rownames = "target_id") %>%
        dplyr::left_join(t2g, by = "target_id") %>%
        dplyr::group_by(.data[["gene_id"]]) %>%
        dplyr::summarise(dplyr::across(-.data[["target_id"]], sum)) %>%
        (\(x){as.matrix(x[,-1]) %>% `rownames<-`(dplyr::pull(x,"gene_id"))})
    } else {
      counts_raw <- count_matrix
    }
    overdispersion.prior <- nboot <- NULL
  }

  data_long <-
    counts_raw %>%
    tidyr::as_tibble(rownames = "target_id") %>%
    tidyr::pivot_longer(-.data$target_id, names_to = "sample", values_to = "counts_raw") %>%
    dplyr::mutate(counts = counts_raw) %>%
    dplyr::left_join(exp_design %>% dplyr::select(-tidyr::any_of(c("path"))), by = "sample") %>%
    dplyr::arrange(.data$target_id)

  if(aggregate_to_gene){
    data_long <- data_long %>% dplyr::left_join(t2g, by = c("target_id"))
  } else {
    data_long <- data_long %>% dplyr::mutate(gene_id = .data$target_id)
  }

  if(bootstrap){
    data_long <-
      data_long %>%
      dplyr::mutate(overdispersions = catch$overdispersion[.data$target_id],
                    counts_scaled = counts_raw/.data$overdispersions,
                    counts = .data$counts_scaled) %>%
      dplyr::left_join(boot, by = c("sample","target_id"))
  }

  message("Filtering low count genes")
  target_to_keep <- do.call(filter_fun, args = c(list(data = data_long),filter_fun_args))

  #if(aggregate_to_gene){
  #  genes_to_keep <-  dplyr::filter(t2g, .data$target_id %in% target_to_keep) %>% dplyr::pull(.data$gene_id) %>% unique
  #  target_to_keep <-
  #    data_long %>%
  #    dplyr::filter(.data$gene_id %in% genes_to_keep) %>%
  #    dplyr::filter(sum(round(.data$counts))>0, .by = "target_id") %>%
  #    dplyr::pull(.data$target_id) %>% unique
  #}

  count_matrix_filtered <- counts_raw[target_to_keep,]
  if(normalize){
    norm_factor <- DESeq2::estimateSizeFactorsForMatrix(count_matrix_filtered)
  } else norm_factor <- rep(1, nrow(exp_design)) %>% `names<-`(exp_design$sample)

  data_long <-
    data_long %>%
    dplyr::filter(.data$target_id %in% target_to_keep) %>%
    dplyr::mutate(norm_factor = norm_factor[sample]) %>%
    dplyr::relocate(.data$target_id)

  if(regularize){
    if(bootstrap){
      dispersions <- estimate_dispersions(count_matrix_filtered/catch$overdispersion[rownames(count_matrix_filtered)],exp_design)
    } else {
      dispersions <- estimate_dispersions(count_matrix_filtered,exp_design)
    }
    data_long <-
      data_long %>%
      dplyr::mutate(disp = dispersions[.data$target_id] %>% as.numeric)
  }

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
    dplyr::pull(.data$target_id)
}

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
