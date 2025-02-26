
#' Launches a Shiny app to visualise the data and fitted models of a cpam object
#'
#' @param cpo A cpam object containing count data, model fits, and optional changepoint/shape estimates
#' @param subset Character vector; names of targets or genes (if `cpo$gene_level = TRUE`)
#'        to load into the Shiny app. If NULL, all genes/targets are included based on `degs_only`.
#' @param degs_only Logical; if TRUE, display only differentially expressed genes/targets
#'        with adjusted p-value below `deg_threshold`. Default is TRUE.
#' @param deg_threshold Numeric; significance threshold for differentially expressed genes/targets.
#'        Only used when `degs_only = TRUE`. Default is 0.05.
#' @param p_type character; choose the type of p-value. Options are "p_gam" (default)
#'  or "p_mvn" (see [`compute_p_values()`] for details).
#' @param shape_type character; "shape1" to include unconstrained or otherwise "shape2".
#' Default is "shape1". In some instances, all of the transcripts for a gene may be "null" shaped,
#' but the p-value for the gene may still be significant. This is due to the different
#' methods of determining significance for the changepoints and the gene-level p-values.
#' Here, conservatively, we remove these null-shaped genes from the DEG list.
#'
#' @return None (launches Shiny app in browser)
#' @export
#' @aliases visualize
#'
#' @examples
#' \dontrun{
#'
#' # Launch visualization with all genes
#' visualise(cpo, degs_only = FALSE)
#'
#' # Launch with only significant genes
#' visualise(cpo, deg_threshold = 0.05)
#'
#' # Launch with specific genes
#' visualise(cpo, subset = c("transcript_1", "transcript_2"))
#' }
visualise <- function(cpo,
                      subset = NULL,
                      degs_only = T,
                      deg_threshold = 0.05,
                      p_type = c("p_gam","p_mvn"),
                      shape_type = c("shape1","shape2")){

  p_type <- match.arg(p_type)
  shape_type <- match.arg(shape_type)

  if(is.null(subset)){
    if(degs_only & !is.null(cpo$p_table)){
      if(cpo$gene_level){
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_target <= deg_threshold) %>%
          dplyr::pull(.data$target_id)
        cpo$data_long <- dplyr::filter(cpo$data_long,.data$target_id %in% subset)
      } else {
        subset <-
          cpo$p_table %>%
          {
            if(!is.null(cpo$shapes)){
              dplyr::left_join(.,cpo$shapes %>% dplyr::select(.data$target_id, shape = {{shape_type}}),
                               by = "target_id") %>%
                dplyr::filter(!all(.data$shape=="null"),.by = "gene_id")
            } else .
          } %>%
          dplyr::filter(.data$q_val_gene <= deg_threshold) %>%
          dplyr::pull(.data$gene_id)
        cpo$data_long <- dplyr::filter(cpo$data_long,.data$gene_id %in% subset)
      }
      if(length(subset)==0) return(message("There are no differentially genes to plot for the given threshold.
                                         Try a different threshold or set 'degs_only = F'"))
    }
  } else{
    if(!all(subset %in% cpo$data_long$gene_id)) stop("one or more of the supplied gene_ids are invalid")
    cpo$data_long <- dplyr::filter(cpo$data_long,.data$gene_id %in% subset)
  }

  app(cpo)
} # end visualize

visualize <- visualise


#' Internal Shiny App for CPAM Visualization
#'
#' @param cpo A CPAM object containing fitted models and expression data
#' @return A Shiny app object
#' @keywords internal
app_old <- function(cpo){

  genes_all <- cpo$data_long %>% dplyr::pull(.data$gene_id) %>% unique
  times <- cpo$times
  times <- c("auto",times)

  bss <- c(auto = "auto",
           null = "null",
           `log-linear (lin)` = "lin",
           `convex (cx)` = "cx",
           `concave (cv)` = "cv",
           #`increasing convex (micx)` = "micx",
           `increasing concave (icv)` = "micv",
           `decreasing convex (dcx)` = "mdcx",
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
                                       inline = T, selected = c("mean","CI")),
                           shiny::selectizeInput("gene_id", "Gene ID", choices = genes_all[1], selected = genes_all[1]),
                           shiny::selectInput("tx_id", "Transcript ID", choices = "all"),
                           #shiny::checkboxInput('gene_level', 'Aggregate counts to gene level', value = F),
                           shiny::checkboxInput('remove_null', 'Plot DETs only', value = F),
                           shiny::numericInput("remove_null_threshold", "Adjusted P value threshold for DETs",
                                 value = 0.05, min = 0.01, max = 0.99, step = 0.01),
                           #shiny::checkboxInput('logged', 'Log transformed', value = T),
                           #shiny::checkboxInput('diff_var', 'Differential variance', value = F),
                           #shiny::numericInput("diff_var_alpha", "p-value threshold for differential variance",
                           #      value = 0.05, min = 0.01, max = 0.99, step = 0.01),
                           shiny::checkboxInput('facet', 'Transcripts in separate plots', value = T),
                           shiny::checkboxInput('common_y_scale', 'Common scale for y-axis', value = T),
                           shiny::checkboxInput('shape_type', 'Include \'unconstrained\' (\'tp\') as a candidate shape', value = T),
                           shiny::hr(style="border-color: black;"),
                           shiny::strong("Manual settings for fitted trends"),
                           shiny::br(),shiny::br(),
                           shiny::selectInput("cp_fix", "Changepoint", choices = times,selected = "auto"),
                           shiny::selectInput("bs", "Shape", choices = bss, selected = "auto"),
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
                            choices = c("all",(if(cpo$gene_level){input$gene_id } else { cpo$t2g %>%
                                                 dplyr::filter(.data$gene_id == input$gene_id) %>%
                                                 dplyr::arrange(.data$target_id) %>%
                                                 dplyr::pull(.data$target_id)})),
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
    plot_cpam(cpo,
              gene_id = input$gene_id,
              target_id = if(input$tx_id=="all") {NULL} else input$tx_id,
              #gene_level_plot = input$gene_level,
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
              facet = input$facet,
              common_y_scale = input$common_y_scale
              )

  }, height = 900, width = 1000) #, outputArgs = list(height = "100%")
}#end server

# Run the application
shiny::shinyApp(ui = ui, server = server,  options = list(launch.browser = T, height = "100%", width = 1800))

}



#------------------------------

app <- function(cpo) {
  # Prepare data
  genes_all <- cpo$data_long %>% dplyr::pull(.data$gene_id) %>% unique()
  times <- cpo$times
  times <- c("auto", times)

  bss <- c(
    auto = "auto",
    null = "null",
    `log-linear (lin)` = "lin",
    `convex (cx)` = "cx",
    `concave (cv)` = "cv",
    `increasing concave (micv)` = "micv",
    `decreasing convex (mdcx)` = "mdcx",
    `unconstrained (tp)` = "tp"
  )

  # Helper function for tooltips
  tooltip <- function(id, title) {
    bslib::tooltip(
      id,
      title,
      placement = "right"
    )
  }

  # UI Definition
  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(
      version = 5,
      bootswatch = "litera",
      bg = "#fff",
      fg = "#000",
      info = "#4582EC",
      #primary = "#0d6efd",
      "enable-shadows" = TRUE,
      "spacer" = "1.5rem"
    ),
    shiny::withMathJax(),
    shinyjs::useShinyjs(),

    # Add custom CSS
    shiny::tags$head(
      shiny::tags$style(
        shiny::HTML("
  .card { margin-bottom: 1rem; }
  .form-group { margin-bottom: 1rem; }
  .tooltip { font-size: 14px; }
  button.collapse-button[aria-expanded='false'] .collapsed-icon { display: inline-block; }
  button.collapse-button[aria-expanded='false'] .expanded-icon { display: none; }
  button.collapse-button[aria-expanded='true'] .collapsed-icon { display: none; }
  button.collapse-button[aria-expanded='true'] .expanded-icon { display: inline-block; }
  .shiny-input-container .checkbox-inline {
    margin-right: 30px !important;
    padding-right: 10px;
  }
  /* Make dropdowns appear above other content */
  .selectize-dropdown,
  .select-dropdown {
    z-index: 1000 !important;
    position: absolute !important;
  }

  /* Ensure the card doesn't clip the dropdowns */
  .card-body {
    overflow: visible !important;
  }
  .card {
    overflow: visible !important;
  }
")

      )
    ),

    shiny::navbarPage(
      title = "CPAM Visualization",
      id = "nav",

      shiny::tabPanel(
        "Plot",
        shiny::sidebarLayout(
          shiny::sidebarPanel(

            bslib::card(
              bslib::card_header("Gene Selection"),
              bslib::card_body(
                shiny::selectizeInput("gene_id",
                                      label = shiny::tags$span("Gene ID"),
                                      choices = genes_all[1],
                                      selected = genes_all[1],
                                      width = "100%"
                ),
                shiny::selectInput("tx_id",
                                   label = shiny::tags$span("Transcript ID"),
                                   choices = "all",
                                   width = "100%"
                )
              )
            ),

            # Display Options Card
            bslib::card(
              bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                "Display Options",
                shiny::tags$button(
                  class = "btn btn-link collapse-button",
                  type = "button",
                  `data-bs-toggle` = "collapse",
                  `data-bs-target` = "#displayOptionsCard",
                  `aria-expanded` = "false",
                  style = "box-shadow: none;",
                  shiny::icon("plus", class = "collapsed-icon"),
                  shiny::icon("minus", class = "expanded-icon")
                )
              ),
              bslib::card_body(
                class = "collapse",
                id = "displayOptionsCard",
                shiny::div(style = "margin-bottom: -1rem;",
                           shiny::div(class = "d-flex align-items-baseline",
                shiny::tags$span("Trend:", class = "me-3"),
                shiny::checkboxGroupInput("fits",
                                          label = NULL,
                                          c(`mean` = "mean", `standard error` = "CI"),
                                          inline = TRUE,
                                          selected = c("mean", "CI")
                    )
                  )
                ),
                shiny::div(class = "d-flex align-items-baseline",
                    shiny::tags$span("Data:", class = "me-3"),
                    shiny::checkboxGroupInput("datas",
                                          label = NULL,
                                          c(`mean` = "mean", `standard error` = "CI"),
                                          inline = TRUE,
                                          selected = c("mean", "CI")
                    )
                ),
                shiny::checkboxInput('facet',
                                     'Transcripts in separate plots',
                                     value = TRUE
                ),
                shiny::checkboxInput('common_y_scale',
                                     'Common y scale',
                                     value = FALSE
                )
              )
            ),

            # Filtering Card
            bslib::card(
              bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                "Filtering",
                shiny::tags$button(
                  class = "btn btn-link collapse-button",
                  type = "button",
                  `data-bs-toggle` = "collapse",
                  `data-bs-target` = "#filteringCard",
                  `aria-expanded` = "false",
                  style = "box-shadow: none;",
                  shiny::icon("plus", class = "collapsed-icon"),
                  shiny::icon("minus", class = "expanded-icon")
                )
              ),
              bslib::card_body(
                class = "collapse",
                id = "filteringCard",
                shiny::checkboxInput('remove_null', 'Remove non-significant transcripts',
                                     value = FALSE
                ),
                shiny::numericInput("remove_null_threshold",
                                    "Adjusted p-value threshold",
                                    value = 0.05, min = 0.01, max = 0.99,
                                    step = 0.01
                )
              )
            ),

            # Manual Fitting Options Card
            bslib::card(
              bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                "Advanced Options",
                shiny::tags$button(
                  class = "btn btn-link collapse-button",
                  type = "button",
                  `data-bs-toggle` = "collapse",
                  `data-bs-target` = "#manualFittingCard",
                  `aria-expanded` = "false",
                  style = "box-shadow: none;",
                  shiny::icon("plus", class = "collapsed-icon"),
                  shiny::icon("minus", class = "expanded-icon")
                )
              ),
              bslib::card_body(
                class = "collapse",
                id = "manualFittingCard",
                shiny::checkboxInput('shape_type',
                                     "Include unconstrained (tp) as candidate shape",
                                     value = TRUE
                ),
                shiny::selectInput("cp_type", "Changepoint selection method ('auto')",
                                   choices = c(`One standard-error rule`="cp_1se",`Minimum score` = "cp_min"),
                                   selected = "cp_1se"
                ),
                shiny::selectInput("cp_fix", "Changepoint (set manually)",
                                   choices = times,
                                   selected = "auto"
                ),
                shiny::selectInput("bs", "Shape (set manually)",
                                   choices = bss,
                                   selected = "auto"
                ),
                shiny::checkboxInput('sp_select',
                                     'Manually set smoothness',
                                     value = FALSE
                ),
                shiny::conditionalPanel(
                  condition = "input.sp_select == true",
                  shiny::sliderInput("sp", "Smoothness",
                                     min = 0, max = 0.25,
                                     value = 0.02,
                                     step = NULL
                  )
                )
              )
            ),

            # Download buttons
            bslib::card(
              bslib::card_header("Export"),
              shiny::downloadButton("downloadPlot", "Download Plot",
                                    class = "btn-primary mb-2 w-100"),
              shiny::downloadButton("downloadData", "Download Data",
                                    class = "btn-secondary w-100")
            ),
            width = 4
          ),

          shiny::mainPanel(
            bslib::card(
              full_screen = TRUE,
              bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                "Results",
                shiny::uiOutput("plotTitle"),
                shiny::div(class = "d-flex align-items-center gap-2",  # Right side with font controls
                    shiny::tags$span("Size:"),
                    shiny::actionButton("decrease_font", "", icon = shiny::icon("minus"),
                                 class = "btn btn-outline-secondary btn-sm"),
                    shiny::tags$span(
                      shiny::textOutput("current_font_size", inline = TRUE)
                    ),
                    shiny::actionButton("increase_font", "", icon = shiny::icon("plus"),
                                 class = "btn btn-outline-secondary btn-sm")
                )
              ),
              shiny::plotOutput("myPlot", height = "800px"),
              shiny::htmlOutput("plotInfo")
            ),
            width = 8
          )
        )
      )

      # shiny::tabPanel(
      #   "Help",
      #   bslib::card(
      #     bslib::card_header("How to use this app"),
      #     bslib::card_body(
      #       shiny::markdown("
      #         ## CPAM Visualization Tool
      #
      #         This tool allows you to visualize and explore your CPAM (Change Point Analysis for Multiple conditions) results.
      #
      #         ### Key Features
      #         - Visualize gene expression patterns over time
      #         - Compare different fitted models
      #         - Explore confidence intervals
      #         - Filter by significance
      #         - Export results
      #
      #         ### Getting Started
      #         1. Select a gene from the dropdown menu
      #         2. Choose display options for fitted trends and data
      #         3. Adjust filtering and plot settings as needed
      #
      #         ### Need Help?
      #         If you encounter any issues or have questions, please refer to the package documentation
      #         or contact the package maintainer.
      #       ")
      #     )
      #   )
      # )
    )
  )

  # Server Definition
  server <- function(input, output, session) {
    # Clean up on session end
    session$onSessionEnded(function() {
      shiny::stopApp()
    })

    loading <- shiny::reactiveVal(TRUE)

    # Add reactive value for font size
    font_size <- shiny::reactiveVal(16)  # Default size

    # Update font size display
    output$current_font_size <- shiny::renderText({
      paste0(font_size(), "px")
    })

    # Increase font size button
    shiny::observeEvent(input$increase_font, {
      current <- font_size()
      if(current < 24) {  # Maximum size
        font_size(current + 2)
      }
    })

    # Decrease font size button
    shiny::observeEvent(input$decrease_font, {
      current <- font_size()
      if(current > 8) {  # Minimum size
        font_size(current - 2)
      }
    })

    # Reactive values for storing state
    rv <- shiny::reactiveValues(
      plot_data = NULL,
      current_gene = NULL,
      loading_genes = TRUE,
      loading_transcripts = TRUE
    )

    # Update gene selection dropdown with loading state
    shiny::observe({
      shiny::req(genes_all)  # Wait for genes_all to be available

      if(rv$loading_genes) {
        shiny::updateSelectizeInput(session, 'gene_id',
                                    choices = c("Loading genes..." = ""),
                                    server = TRUE
        )
        shiny::invalidateLater(100)

      } else {
        shiny::updateSelectizeInput(session, 'gene_id',
                                    choices = genes_all,
                                    selected = genes_all[1],
                                    server = TRUE
        )
      }
      rv$loading_genes <- FALSE
    })

    # Update transcript selection based on gene
    shiny::observe({
      shiny::req(input$gene_id != "", !rv$loading_genes)

      if(rv$loading_transcripts) {
        shiny::updateSelectInput(session, 'tx_id',
                                 choices = c("Loading transcripts..." = ""),
                                 selected = ""
        )
        shiny::invalidateLater(100)
      } else {
        tx_choices <- if(cpo$gene_level) {
          input$gene_id
        } else {
          cpo$t2g %>%
            dplyr::filter(.data$target_id %in% cpo$target_to_keep) %>%
            dplyr::filter(.data$gene_id == input$gene_id) %>%
            dplyr::arrange(.data$target_id) %>%
            dplyr::pull(.data$target_id)
        }

        shiny::updateSelectInput(session, 'tx_id',
                                 choices = c("all", tx_choices),
                                 selected = "all"
        )
      }
      rv$loading_transcripts <- FALSE
    })


    # # Update gene selection dropdown
    # shiny::observe({
    #   shiny::updateSelectizeInput(session, 'gene_id',
    #                               choices = genes_all,
    #                               server = TRUE
    #   )
    # })

    # Update transcript selection based on gene
    # shiny::observe({
    #   tx_choices <- if(cpo$gene_level) {
    #     input$gene_id
    #   } else {
    #     cpo$t2g %>%
    #       dplyr::filter(.data$gene_id == input$gene_id) %>%
    #       dplyr::arrange(.data$target_id) %>%
    #       dplyr::pull(.data$target_id)
    #   }
    #
    #   shiny::updateSelectInput(session, 'tx_id',
    #                            choices = c("all", tx_choices),
    #                            selected = "all"
    #   )
    # })

    # Disable controls based on conditions
    shiny::observe({
      if(cpo$gene_level) {
        shinyjs::disable("tx_id")
        shinyjs::disable("facet")
        shinyjs::disable("remove_null")
        shinyjs::disable("remove_null_threshold")
      }
    })

    # Generate plot title
    output$plotTitle <- shiny::renderUI({
      shiny::HTML(sprintf(
        "Visualizing %s %s",
        if(input$tx_id == "all") "all transcripts for" else paste(input$tx_id,"for gene "),
        input$gene_id
      ))
    })

    # Main plot
    output$myPlot <- shiny::renderPlot({

      shiny::req(
        input$gene_id != "",
        input$tx_id != "",
        !rv$loading_genes,
        !rv$loading_transcripts
      )

      shiny::withProgress(
        message = 'Generating plot',
        detail = 'This may take a while...',
        value = 0,
        {
          plot_cpam(cpo,
                    gene_id = input$gene_id,
                    target_id = if(input$tx_id=="all") {NULL} else input$tx_id,
                    show_fit = "mean" %in% input$fits,
                    show_fit_ci = "CI" %in% input$fits,
                    show_data = "mean" %in% input$datas,
                    show_data_ci = "CI" %in% input$datas,
                    shape_type = if(input$shape_type) "shape1" else "shape2",
                    cp_type = input$cp_type,
                    remove_null = input$remove_null,
                    null_threshold = input$remove_null_threshold,
                    sp = if(input$sp_select){(input$sp)} else {NULL},
                    cp_fix = if(input$cp_fix!="auto"){as.numeric(input$cp_fix)} else {-999},
                    bs = input$bs,
                    facet = input$facet,
                    common_y_scale = input$common_y_scale,
                    base_size = font_size()
          )
        }
      )
    }, height = 800, width = 1000)

    # Plot information
    output$plotInfo <- shiny::renderText({
      if(!is.null(input$gene_id)) {
        sprintf(
          "<div class='alert alert-info'>
            Currently displaying: %s<br>
            Selected shape: %s<br>
            Changepoint: %s
          </div>",
          input$gene_id,
          input$bs,
          input$cp_fix
        )
      }
    })

    # Download handlers
    output$downloadPlot <- shiny::downloadHandler(
      filename = function() {
        paste0("cpam_plot_", input$gene_id, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        grDevices::pdf(file, width = 12, height = 8)
        print(plot_cpam(cpo,
                        gene_id = input$gene_id,
                        target_id = if(input$tx_id=="all") {NULL} else input$tx_id,
                        show_fit = "mean" %in% input$fits,
                        show_fit_ci = "CI" %in% input$fits,
                        show_data = "mean" %in% input$datas,
                        show_data_ci = "CI" %in% input$datas,
                        shape_type = if(input$shape_type) "shape1" else "shape2",
                        remove_null = input$remove_null,
                        null_threshold = input$remove_null_threshold,
                        sp = if(input$sp_select){(input$sp)} else {NULL},
                        cp_fix = if(input$cp_fix!="auto"){as.numeric(input$cp_fix)} else {-999},
                        bs = input$bs,
                        facet = input$facet,
                        common_y_scale = input$common_y_scale
        ))
        grDevices::dev.off()
      }
    )

    output$downloadData <- shiny::downloadHandler(
      filename = function() {
        paste0("cpam_data_", input$gene_id, "_", format(Sys.time(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        data <- cpo$data_long %>%
          dplyr::filter(.data$gene_id == input$gene_id)
        utils::write.csv(data, file, row.names = FALSE)
      }
    )
  }

  # Run the application
  shiny::shinyApp(
    ui = ui,
    server = server,
    options = list(
      launch.browser = TRUE,
      height = "100%",
      width = "100%"
    )
  )
}

