
visualize <- function(tso, 
                      degs_only = T, 
                      deg_threshold = 0.05, 
                      deg_type = "qval_gene"){
  
  if ( !require('shiny') ) {
    stop("'visualize()' requires 'shiny'. Please install it using install.packages('shiny')")
  }
  
  if(degs_only){
    degs <- 
      tso$results %>% 
      filter(get(deg_type) < deg_threshold) %>% 
      pull(gene_id) %>% 
      unique()
    
    tso$data_filtered <- 
      tso$data_filtered %>% 
      filter(gene_id %in% degs)
    
    if(length(degs)==0) return(message("There are no differentially genes to plot for the given threshold. Try a different threshold or set 'degs_only = F'"))
      
  }
  app(tso)
} # end visualize


app <- function(tso){
  

  genes_all <- tso$data_filtered %>% pull(gene_id) %>% unique
  times <- tso$times
  tp3 <- length(times)==3
  if(tp3) times <- c(`change-change` = times[1],
                     `change-flat` = -times[2],
                     `flat-change` = times[2],
                     `null` = times[3])
  times <- c("auto",times)
  
  bss <- c(auto = "auto",
           null = "null",
           `convex (cx)` = "cx",
           `concave (cv)` = "cv",
           `increasing convex (micx)` = "micx",
           `increasing concave (micv)` = "micv",
           `decreasing convex (mdcx)` = "mdcx",
           `decreasing concave (mdcv)` = "mdcv",
           `unconstrained (tp)` = "tp")


# Define UI
ui <- fluidPage(withMathJax(),
                useShinyjs(),
                sidebarLayout(
                  sidebarPanel(
                    checkboxGroupInput("fits", "Fitted trend", c(`mean trend` = "mean", `2SE` = "SE"), 
                                       inline = T, selected = c("mean","SE")),
                    checkboxGroupInput("datas", "Data", c(`mean` = "mean", `2SE (inferential replicates)` = "SE"), 
                                       inline = T, selected = c("mean")),
                    selectizeInput("gene_id", "Gene ID", choices = genes_all[1], selected = genes_all[1]),
                    selectInput("tx_id", "Transcript ID", choices = "all"),
                    checkboxInput('gene_level', 'Aggregate counts to gene level', value = F),
                    checkboxInput('logged', 'Log transformed', value = T),
                    checkboxInput('diff_var', 'Differential variance', value = F),
                    numericInput("diff_var_alpha", "p-value threshold for differential variance",
                                 value = 0.05, min = 0.01, max = 0.99, step = 0.01),
                    checkboxInput('facet', 'Transcripts in separate plots', value = T),
                    hr(style="border-color: black;"),
                    strong("Manual settings for fitted trends"),
                    br(),br(),
                    selectInput("cp_fix", "Changepoint", choices = times,selected = "AUTO"),
                    selectInput("bs", "Shape", choices = bss),
                    checkboxInput('sp_select', 'Set smoothness', value = F),
                    sliderInput("sp", "Smoothness", 
                                min = 0, max = 0.15, value = 0.02, 
                                step = NULL,
                                animate = animationOptions(interval = 100)),
                    width = 4
                  ),
                  
                  mainPanel(
                    plotOutput("myPlot"),
                    width = 15
                  )
                )
)


# Define server
server <-  function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })

  
  observe(updateSelectizeInput(session, 'gene_id', 
                       choices = genes_all, 
                       server = TRUE))
  observe(updateSelectInput(session, 'tx_id',
                            choices = c("all",(tso$t2g %>% 
                                                 filter(gene_id == input$gene_id) %>% 
                                                 arrange(tx_id) %>% 
                                                 pull(tx_id))),
                            selected = "all"))
  
  observe(if(tp3){
    disable("sp_select")
    disable("sp")
    disable("bs")
    disable("diff_var")
    disable("diff_var_alpha")
  })
  
  #observe(if(!tso$bootstrap) disable(selector = '#datas .checkbox:nth-child(SE) label'))
  
  output$myPlot <- renderPlot({
    plot_gene_co(tso,
              gene_id_ = input$gene_id,
              tx_id_ = input$tx_id,
              gene_level = input$gene_level,
              show.fit = "mean" %in% input$fits,
              show.fit.se = "SE" %in% input$fits,
              show.data = "mean" %in% input$datas,
              show.data.se = "SE" %in% input$datas,
              diff_var = input$diff_var,
              diff_var_alpha = input$diff_var_alpha,
              logged = input$logged,
              sp = if(input$sp_select){(input$sp)} else {NULL},
              cp_fix = if(input$cp_fix!="auto"){as.numeric(input$cp_fix)} else {-999},
              bs = input$bs,
              facet = input$facet)
    
  }, height = 900, width = 700) #, outputArgs = list(height = "100%")
}#end server

# Run the application 
shinyApp(ui = ui, server = server,  options = list(launch.browser = T, height = 1200, width = 1800))


}

