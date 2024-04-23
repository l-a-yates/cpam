
#' Launches a Shiny app to visualise the data and fitted models of a cpam object
#'
#' @param cpo a cpam object
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#' to load into Shiny app
#' @param degs_only logical; display DEGs only
#' @param deg_threshold numerical; threshold for DEGs
#'
#' @return None
#' @export
#'
#' @examples 1+1
visualize <- function(cpo,
                      subset = NULL,
                      degs_only = T,
                      deg_threshold = 0.05){

  if(is.null(subset)){
    if(degs_only & !is.null(cpo$p_table)){
      if(cpo$gene_level){
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_target <= deg_threshold) %>%
          dplyr::pull(.data$target_id)
        cpo$data_long <- dplyr::filter(cpo$data_long,.data$target_id %in% subset)
      } else {
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_gene <= deg_threshold) %>%
          dplyr::pull(.data$gene_id)
        cpo$data_long <- dplyr::filter(cpo$data_long,.data$gene_id %in% subset)
      }
      if(length(subset)==0) return(message("There are no differentially genes to plot for the given threshold.
                                         Try a different threshold or set 'degs_only = F'"))
    }
  }

  app(cpo)
} # end visualize


app <- function(cpo){

  genes_all <- cpo$data_long %>% dplyr::pull(.data$gene_id) %>% unique
  times <- cpo$times
  times <- c("auto",times)

  bss <- c(auto = "auto",
           null = "null",
           `convex (cx)` = "cx",
           `concave (cv)` = "cv",
           #`increasing convex (micx)` = "micx",
           `increasing concave (micv)` = "micv",
           `decreasing convex (mdcx)` = "mdcx",
           #`decreasing concave (mdcv)` = "mdcv",
           `unconstrained (tp)` = "tp"
           )


# Define UI
ui <- shiny::fluidPage(shiny::withMathJax(),
                       shinyjs::useShinyjs(),
                       shiny::sidebarLayout(
                         shiny::sidebarPanel(
                           shiny::checkboxGroupInput("fits", "Fitted trend", c(`mean trend` = "mean", `CI` = "CI"),
                                       inline = T, selected = c("mean","CI")),
                           shiny::checkboxGroupInput("datas", "Data", c(`mean` = "mean", `bootstrapped quantile` = "CI"),
                                       inline = T, selected = c("mean")),
                           shiny::selectizeInput("gene_id", "Gene ID", choices = genes_all[1], selected = genes_all[1]),
                           shiny::selectInput("tx_id", "Transcript ID", choices = "all"),
                           shiny::checkboxInput('gene_level', 'Aggregate counts to gene level', value = F),
                           shiny::checkboxInput('remove_null', 'Plot DETs only', value = F),
                           shiny::numericInput("remove_null_threshold", "Adjusted P value threshold for DETs",
                                 value = 0.05, min = 0.01, max = 0.99, step = 0.01),
                           #shiny::checkboxInput('logged', 'Log transformed', value = T),
                           #shiny::checkboxInput('diff_var', 'Differential variance', value = F),
                           #shiny::numericInput("diff_var_alpha", "p-value threshold for differential variance",
                           #      value = 0.05, min = 0.01, max = 0.99, step = 0.01),
                           shiny::checkboxInput('facet', 'Transcripts in separate plots', value = T),
                           shiny::checkboxInput('shape_type', 'Include \'unconstrained\' (\'tp\') as a candidate shape', value = T),
                           shiny::hr(style="border-color: black;"),
                           shiny::strong("Manual settings for fitted trends"),
                           shiny::br(),shiny::br(),
                           shiny::selectInput("cp_fix", "Changepoint", choices = times,selected = "AUTO"),
                           shiny::selectInput("bs", "Shape", choices = bss),
                           shiny::checkboxInput('sp_select', 'Set smoothness', value = F),
                           shiny::sliderInput("sp", "Smoothness",
                                min = 0, max = 0.15, value = 0.02,
                                step = NULL,
                                animate = shiny::animationOptions(interval = 100)),
                    width = 4
                  ),

                  shiny::mainPanel(
                    shiny::plotOutput("myPlot"),
                    width = 15
                  )
                )
)


# Define server
server <-  function(input, output, session) {

  session$onSessionEnded(function() {
    shiny::stopApp()
  })


  shiny::observe(shiny::updateSelectizeInput(session, 'gene_id',
                       choices = genes_all,
                       server = TRUE))
  shiny::observe(shiny::updateSelectInput(session, 'tx_id',
                            choices = c("all",(cpo$t2g %>%
                                                 dplyr::filter(.data$gene_id == input$gene_id) %>%
                                                 dplyr::arrange(.data$target_id) %>%
                                                 dplyr::pull(.data$target_id))),
                            selected = "all"))

  shiny::observe(if(cpo$gene_level){
    shinyjs::disable("tx_id")
    shinyjs::disable("facet")
    shinyjs::disable("remove_null")
    shinyjs::disable("remove_null_threshold")
    shinyjs::disable("gene_level")
  })
  #shiny::observe(if(!cpo$bootstrap) shinyjs::disable("tx_id"))
  #shiny::observe(shinyjs::disable("diff_var"))
  #shiny::observe(shinyjs::disable("diff_var_alpha"))

  # observe(if(tp3){
  #   disable("sp_select")
  #   disable("sp")
  #   disable("bs")
  #   disable("diff_var")
  #   disable("diff_var_alpha")
  # })

  #observe(if(!cpo$bootstrap) disable(selector = '#datas .checkbox:nth-child(SE) label'))

  output$myPlot <- shiny::renderPlot({
    plot_gene_co(cpo,
              gene_id = input$gene_id,
              target_id = if(input$tx_id=="all") {NULL} else input$tx_id,
              gene_level_plot = input$gene_level,
              show_fit = "mean" %in% input$fits,
              show_fit_ci = "CI" %in% input$fits,
              show_data = "mean" %in% input$datas,
              show_data_ci = "CI" %in% input$datas,
              shape_type = if(input$shape_type) "shape1" else "shape2",
              #diff_var = input$diff_var,
              #diff_var_alpha = input$diff_var_alpha,
              #logged = input$logged,
              remove_null = input$remove_null,
              null_threshold = input$remove_null_threshold,
              sp = if(input$sp_select){(input$sp)} else {NULL},
              cp_fix = if(input$cp_fix!="auto"){as.numeric(input$cp_fix)} else {-999},
              bs = input$bs,
              facet = input$facet)

  }, height = 900, width = 700) #, outputArgs = list(height = "100%")
}#end server

# Run the application
shiny::shinyApp(ui = ui, server = server,  options = list(launch.browser = T, height = 1200, width = 1800))


}

