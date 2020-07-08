#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(dplyr)
library(feather)
library(ggplot2)
library(pheatmap)
library(DT)
library(rcytoscapejs2)
library(glue)
load("data/joint_cortex/cortex_prep.Rda") # a list, data_cortex
load("data/joint_pons/pons_prep.Rda")     # a list, data_pons
load("data/joint_cortex/common_prep.Rda") # metadata and colour_palettes
source("functions.R")

# datasets

#---------------------------------------------------------------------------
# TF has ext and weight suffix in these datas
# a vector that contains all TFs(ext or regular) in the activity data(by cluster/cell)
#TF_active <- read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds")

# import feather files later based on input$TF
#activity_cluster <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cluster.feather")
#activity_cell <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather")

# metadata, load in data_prep.R
#metadata <- read_tsv("data/joint_cortex/metadata_20190716.tsv")


#------------------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("joint cortex and pons data app"),
    
    sidebarLayout(
        sidebarPanel(
            # choose which datasets to analyze for the whole app
            radioButtons("region", "Brain region",
                         # use the names in the vector to display
                         # use the character "joint_cortex" to match the path to import data
                         choices = c("Forebrain" = "cortex",
                                     "Pons" = "pons"),
                         selected = "cortex"),
            
            #actionButton("update_tf", label = "Update transcription factors to see the plots"),
            selectInput(inputId = "TF",
                        label = "transcription factor",
                        choices = data_cortex$unique_active_TFs_bare,
                        multiple = TRUE,
                        selected = c("Arx","Lef1")),
            
            # 1. table and network graph of related TF and genes
            conditionalPanel(condition = "input.tabs == 'table and network'",
                             radioButtons("show", "Node display option",
                                          # use the names in the vector to display
                                          # use the character "joint_cortex" to match the path to import data
                                          choices = c("show all nodes" = "all",
                                                      #"color by user's input tf list" = "pathway",
                                                      "neglect graynodes " = "neglect",
                                                      
                                                      "stop showing" = "stop"),
                                          selected = "stop"),#,
                             checkboxInput("show_pathway","color by user's input genes list",
                                           TRUE),
                             selectInput("input_pathway", "input your interested genes pathway",
                                         choices = data_cortex$unique_active_TFs_bare,
                                         multiple = TRUE,
                                         selected = c("Arx","Lef1"))
                             #actionButton("update_graph", label = "See the network graph")
                             ),
            # 2. heatmap and clustering
            conditionalPanel(condition = "input.tabs == 'heatmap and clustering'",
                             numericInput(inputId = "num_cell_plot", label = "number of cells to visualize",
                                          value = 300),
                             # numericInput(inputId = "num_cluster_plot", label = "number of clusters to visualize",
                             #              value = 50),
                             checkboxGroupInput("method", "Plot by cluster or cells",
                                   choices = c("cluster" = "Cluster",
                                               "cell" = "Cell")             
    
                             )),
            # 3. time series plot
            conditionalPanel(condition = "input.tabs == 'time series'"),
            # Update everything
            actionButton("update", label = "Update"),
        ),
        mainPanel(tabsetPanel(
            
            tabPanel(title = "table and network",
                     dataTableOutput("table"),
                     textOutput("desc"),
                     rcytoscapejsOutput("network", width = "1200px",height = "600px"),
                     
                     value = "table and network"
            ),
            tabPanel("heatmap and clustering",
                     plotOutput("heatmap_cell"),
                     plotOutput("heatmap_cluster"),
                     plotOutput("cluster1"),
                     plotOutput("cluster2"),
                     
                     value = "heatmap and clustering"
            ),
            tabPanel("time series",
                     plotOutput("timeseries1"),
                     plotOutput("timeseries2"),
                     plotOutput("timeseries3"),
                     value = "time series"),
            id = "tabs"
        ))
    ),
    
    
)

server <- function(input, output, session) {
  # Dynamic UI, change the selectInput tf lists on display
  observeEvent(input$region,{
    if(input$region == "cortex"){
      updateSelectInput(session, inputId = "TF", choices = data_cortex$unique_active_TFs_bare, 
                        selected = c("Arx","Lef1"))
      updateSelectInput(session, inputId = "input_pathway", choices = unique(data_cortex$TF_target_gene_info$gene), 
                          selected = c("Dlx6","Sox6") )
      
    }
    else{
      updateSelectInput(session, inputId = "TF", choices = data_pons$unique_active_TFs_bare, 
                        selected = c("Lhx5","Pax7"))
      updateSelectInput(session, inputId = "input_pathway", choices = unique(data_pons$TF_target_gene_info$gene), 
                        selected = c("Gad2"))
      
    }
    updateRadioButtons(session, "show", selected = "stop")
  })
  
  
  
  input_new <- eventReactive(input$update,{
    
    l <- list()
    if(input$region == "cortex"){
      l <- data_cortex
    }
    else if(input$region == "pons"){
      l <- data_pons
      }
    
    l$tf <- input$TF
    l$region <- input$region
    l$input_pathway <- input$input_pathway
    # l has following elements with same names for both options above:
    # l contains ...
    
    # We will use the same name attributes to retrieve data
    return (l)
    })
  
  #input_tf <- reactive(input_new()$tf)
  
  # -----------------------------Tab1:table and network------------------------------------------
    output$table <- renderDataTable({
        # process data, filter the lines with our interested TF
      datatable(dplyr::filter(input_new()$TF_target_gene_info, TF %in% input_new()$tf))
    })
    
    
    output$desc <- renderText({
      text <- "Orange nodes are genes that express its own transcription factors; " %>%
        paste("Purple nodes in the center are your input transcription factors; ") %>%
        paste("grey nodes are other genes.")
    })
    nodeData <- reactive(
      #input$show,
      {

      if(input$show_pathway){input_pathway <- input_new()$input_pathway}
      else{input_pathway <- c()}
      
      if(input$show == "all"){
        create_network(input_new()$tf, input_new()$TF_target_gene_info,
                       input_new()$unique_active_TFs_bare,
                       input_pathway = input_pathway)$nodes
      }
      else if(input$show == "neglect"){
        create_network(input_new()$tf, input_new()$TF_target_gene_info,
                       input_new()$unique_active_TFs_bare,
                       input_pathway = input_pathway)$nodes %>%
          filter(color!="lightgrey")
      }
        
    })
    output$network <- renderRcytoscapejs({
      req(nodeData())
      nodeData <- nodeData()
      edgeData <- create_network(input_new()$tf, input_new()$TF_target_gene_info,
                                 input_new()$unique_active_TFs_bare)$edges
      network <- createCytoscapeJsNetwork(nodeData, edgeData)
      rcytoscapejs2(network$nodes, network$edges)
 
    })
    
    
    # -----------------------------Tab2-------------------------------------------
   
    output$heatmap_cell <- renderPlot({
      req("Cell" %in% input$method)
      plot_heatmap(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext,input_new()$cell_metadata,
                   cell_plot_num = input$num_cell_plot)
    })
    
    output$heatmap_cluster <- renderPlot({
      req("Cluster" %in% input$method)
      plot_heatmap(input_new()$tf, "Cluster",input_new()$region, input_new()$TF_and_ext,input_new()$cell_metadata)
                   
    })
  
    # The cluster scatterplot is always plot by cells, so we use an independent reactive
    # value for this plot
    activity_data_cluster <- reactive({
      # use the feature of feather data to read certain col to optimize speed
      create_activity_data(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext)
    })
    
    output$cluster1 <- renderPlot({
      req(length(input_new()$tf)>0)
      plot_UMAP(tf_number = 1,input_new()$cell_metadata, activity_data_cluster())
      
    })
    
    output$cluster2 <- renderPlot({
      req(length(input_new()$tf)>1)
      plot_UMAP(tf_number = 2,input_new()$cell_metadata, activity_data_cluster())
      
    })
    
    
    
    # --------------------------------------Tab3: timeseries-------------------------------------------
    output$timeseries1 <- renderPlot({
      req(length(input_new()$tf)>0)
      # binary_active_TFs is loaded at beginning by data_prep.R
      TF <- translate_tf(input_new()$tf[1],input_new()$binary_active_TFs)
      req(TF)
      plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity)
 
    })
    output$timeseries2 <- renderPlot({
      req(length(input_new()$tf)>1)
      TF <- translate_tf(input_new()$tf[2],input_new()$binary_active_TFs)
      req(TF)
      plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity)
      
    })
    output$timeseries3 <- renderPlot({
      req(length(input_new()$tf)>2)
      TF <- translate_tf(input_new()$tf[3],input_new()$binary_active_TFs)
      req(TF)
      plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity)
      
    })
    
    
    
    
      
      
    
}
# Run the application 
shinyApp(ui = ui, server = server)
