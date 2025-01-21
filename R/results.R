# scaled: should the counts for `min_count` be scaled or not (only applicable to bootstrapped data) - default is T

#' Title
#'
#' @param cpo a cpam object
#' @param p_threshold numerical; threshold
#' @param p_type character; choose the type of p-value. Options are "p_gam" (default) or "p_mvn".
#' @param min_lfc numerical; maximum absolute log (base 2) fold change must exceed this minimum value; default is 0
#' @param min_count numerical; maximum of the modelled counts evaluated at the set of observed time points must exceed this minimum value for
# at least one isoform - default is 0
#' @param aggregate_to_gene  logical; filter by gene-aggregated p-values
#' @param add_lfc logical; add log (base 2) fold changes for each time point
#' @param add_counts logical; add modelled counts for each time point
#' @param cp_type character; model-selection rule used to select the changepoint
#' @param shape_type character; "shape1" to include unconstrained or otherwise "shape2"
#' @param summarise_to_gene logical; return gene-level results only
#'
#' @return a tibble
#' @export
#'
#' @examples 1+1
results <- function(cpo,
                    p_threshold = 0.05,
                    p_type = c("p_gam","p_mvn"),
                    min_lfc = 0,
                    min_count = 0,
                    aggregate_to_gene = cpo$aggregate_to_gene,
                    add_lfc = T,
                    add_counts = T,
                    cp_type = c("cp_1se","cp_min"),
                    shape_type = c("shape1","shape2"),
                    summarise_to_gene = F){

  p_type <- match.arg(p_type)
  cp_type <- match.arg(cp_type)
  shape_type <- match.arg(shape_type)

  if(aggregate_to_gene){
    if(!cpo$aggregate_to_gene) stop("The cpam object must be prepared with `aggregate_to_gene = T' to use this option")
    level <- "gene"
  } else level <- "target"

  if(summarise_to_gene & !cpo$gene_level){
    if(add_lfc | add_counts) message("add_lfc and add_counts are not available for summarise_to_gene = T")
    add_lfc <- F
    add_counts <- F
  }

  if(p_type == "p_gam"){
    p <- paste0("q_val_",level)
    p_table <-
      cpo$p_table %>%
      dplyr::filter(.data[[p]] < p_threshold) %>%
      dplyr::select(dplyr::any_of(c("target_id","gene_id")), p = {{p}})
  } else{
    p <- paste0("q_mvn_",level)
    p_table <-
      cpo$p_mvn %>%
      dplyr::filter(.data[[p]] < p_threshold) %>%
      dplyr::select(dplyr::any_of(c("target_id","gene_id")), p = {{p}})
  }

  keep_p_filtered <- p_table %>% dplyr::pull(.data$target_id)

  if(min_count >0 | add_counts){
    if(is.null(cpo$pred)) stop("Shapes must be selected before results can be filtered by min count or counts added")
    counts_by_time <-
      cpo$pred %>%
      {`rownames<-`(as.matrix(dplyr::select(.,-"target_id")),.$target_id)}
  }

  if(min_count > 0){
    keep_min_count <-
      matrixStats::rowMaxs(counts_by_time) %>%
      {.[. > min_count]} %>%
      names

    keep_p_filtered <- intersect(keep_p_filtered,keep_min_count)
  }

  if(min_lfc > 0){
    if(is.null(cpo$lfc)) stop("Shapes must be selected before results can be filtered by log-fold change (lfc)")
    keep_min_lfc <-
      cpo$lfc %>%
      {`rownames<-`(as.matrix(dplyr::select(.,-.data$target_id)),.$target_id)} %>%
      matrixStats::rowMaxs() %>%
      {.[abs(.)> min_lfc]} %>%
      names

    keep_p_filtered <- intersect(keep_p_filtered,keep_min_lfc)
  }

  p_table <- p_table %>% dplyr::filter(.data$target_id %in% keep_p_filtered)

  if(!is.null(cpo$changepoints)){
    p_table <-
      p_table %>%
      dplyr::left_join(cpo$changepoints %>% dplyr::select(.data$target_id, cp = {{cp_type}}),
                       by = "target_id")
  }

  if(!is.null(cpo$shapes)){
    p_table <-
      p_table %>%
      dplyr::left_join(cpo$shapes %>% dplyr::select(.data$target_id, shape = {{shape_type}}),
                       by = "target_id")

    if(!is.null(p_table$cp)){
      p_table <-
        p_table %>% dplyr::mutate(cp = dplyr::if_else(.data$shape == "null",max(cpo$times),.data$cp))
    }

  }

  if(add_lfc){
    if(is.null(cpo$lfc)) stop("Shapes must be selected before results can be filtered by log-fold change (lfc)")
    p_table <- p_table %>%
      dplyr::left_join(cpo$lfc %>%
                  dplyr::rename_with(~ paste0("lfc.", .x, recycle0 = TRUE),
                              .cols = -.data$target_id),
                  by = "target_id")
  }

  if(add_counts){
    if(is.null(cpo$lfc)) stop("Shapes must be selected before results can be filtered by log-fold change (lfc)")

    p_table <- p_table %>%
      dplyr::left_join(counts_by_time %>%
                         dplyr::as_tibble(rownames = "target_id") %>%
                         dplyr::rename_with(~ paste0("counts.", .x, recycle0 = TRUE),
                              .cols = -.data$target_id),
                       by = "target_id")
  }


  if(summarise_to_gene & !cpo$gene_level){
    p_table <-  p_table %>% dplyr::select(dplyr::any_of(c("gene_id","p"))) %>% dplyr::distinct()
  }

  p_table %>% dplyr::arrange(.data$p)
}
