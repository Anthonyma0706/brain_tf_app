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
load("data/joint_cortex/cortex_prep.Rda")
source("functions.R")

# datasets

#---------------------------------------------------------------------------


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
    titlePanel("joint cortex app"),
    
    sidebarLayout(
        sidebarPanel(
            # choose which datasets to analyze for the whole app
            radioButtons("region", "Brain region",
                         # use the names in the vector to display
                         # use the character "joint_cortex" to match the path to import data
                         choices = c("Forebrain" = "joint_cortex",
                                     "Pons" = "joint_pons"),
                         selected = "Forebrain"),
            
            actionButton("update_tf", label = "Update transcription factors to see the plots"),
            selectInput(inputId = "TF",
                        label = "transcription factor",
                        choices = unique_TF,
                        multiple = TRUE,
                        selected = c("Arx","Lef1")),
            
            # 1. table and network graph of related TF and genes
            conditionalPanel(condition = "input.tabs == 'table and network'",
                             radioButtons("show", "Node display option",
                                          # use the names in the vector to display
                                          # use the character "joint_cortex" to match the path to import data
                                          choices = c("show all nodes" = "all",
                                                      "neglect graynodes " = "neglect"))#,
                             #actionButton("update_graph", label = "See the network graph")
                             ),
            # 2. heatmap and clustering
            conditionalPanel(condition = "input.tabs == 'heatmap and clustering'",
                             numericInput(inputId = "num_cell", label = "number of cells to visualize",
                                          value = 50),
                             checkboxGroupInput("method", "Plot by cluster or cells",
                                   choices = c("cluster" = "Cluster",
                                               "cell" = "Cell")             
    
                             )),
            # 3. time series plot
            conditionalPanel(condition = "input.tabs == 'time series'")
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

server <- function(input, output) {
  
  input_tf <- eventReactive(input$update_tf,{input$TF})
  
  # -----------------------------Tab1:table and network------------------------------------------
    output$table <- renderDataTable({
        # process data, filter the lines with our interested TF
      datatable(dplyr::filter(TF_target_gene, TF %in% input_tf()))
    })
    
    
    output$desc <- renderText({
      text <- "Orange nodes are genes that express its own transcription factors; " %>%
        paste("Purple nodes in the center are your input transcription factors; ") %>%
        paste("grey nodes are other genes.")
    })
    nodeData <- eventReactive(input$show,{
      req(input$show)
      if(input$show == "all"){
        create_network(input_tf())$nodes
      }
      else{
        create_network(input_tf())$nodes %>%
          filter(color!="lightgrey")
      }
    })
    output$network <- renderRcytoscapejs({
      
      nodeData <- nodeData()
      edgeData <- create_network(input_tf())$edges
      network <- createCytoscapeJsNetwork(nodeData, edgeData)
      rcytoscapejs2(network$nodes, network$edges)
 
    })
    
    
    # -----------------------------Tab2-------------------------------------------
    # activity_data_heatmap <- reactive({
    #   # use the feature of feather data to read certain col to optimize speed
    #   create_activity_data(input_tf(), input$method)
    # })
    # eventReactive(input$update_heatmap, {
    #   
    # })
    
    output$heatmap_cell <- renderPlot({
      req("Cell" %in% input$method)
        act_cell <- create_activity_data(input_tf(), "Cell") %>%
          mutate(Cluster = gsub("_"," ",forebrain_data[["Sample_cluster"]])) %>%
          filter(!grepl("BLACKLIST", Cluster)) %>% # filter out bad samples
          sample_n(300) %>%
          tibble::column_to_rownames(var = "Cell") # make that column name as row name ...
        anno_row_cell <- select(act_cell, Cluster)
        #anno_row_cell$Cluster <- as.factor(anno_row_cell$Cluster)
        act_cell <- select(act_cell, -Cluster)# %>% # must remove Cluster data before plotting
       
        pheatmap::pheatmap(t(act_cell),
                           show_colnames = FALSE,
                           scale = "none",
                           border_color = NA,
                           color = colorRampPalette(c("blue", "white", "red"))(100),
                           main = "Plot by Cells",
                           annotation_col = anno_row_cell,
                           # change the default color annotation
                           annotation_colors = hm_anno$side_colors, 
                           annotation_legend = TRUE,
                           cellwidth = 2,
                           cellheight = 10)
    })
    
    output$heatmap_cluster <- renderPlot({
      req("Cluster" %in% input$method)
        act <- create_activity_data(input_tf(), "Cluster") %>%
          sample_n(30) %>%
          tibble::column_to_rownames(var = "Cluster") # make that column name as row name ...
        
        pheatmap::pheatmap(t(act),
                           #show_colnames = FALSE,
                           scale = "none",
                           border_color = NA,
                           color = colorRampPalette(c("blue", "white", "red"))(100),
                           main = "Plot by Clusters",
                           annotation_col = hm_anno$anno_row,
                           # change the default color annotation
                           annotation_colors = hm_anno$side_colors, 
                           annotation_legend = TRUE,
                           cellwidth = 10,
                           cellheight = 10)
      
      
      
      
    })
  
    # The cluster scatterplot is always plot by cells, so we use an independent reactive
    # value for this plot
    activity_data_cluster <- reactive({
      # use the feature of feather data to read certain col to optimize speed
      create_activity_data(input_tf(), "Cell")
    })
    
    output$cluster1 <- renderPlot({
      req(length(input_tf())>0)
      activity_test1 <- activity_data_cluster()[,2][[1]] # now we only plot the first tf input,
      # more plots could be generated later
      #activity_test1 <- create_activity_data("Arx")[,2][[1]]
      # add a activity column
      forebrain_with_activity <- mutate(forebrain_data, activity_1 = activity_test1) 
      
      ggplot(data = forebrain_with_activity, mapping = aes(x=UMAP1,y=UMAP2))+
        geom_point(aes(color = activity_1))+
        scale_color_gradientn(colors = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100)))+
        theme_bw()
      
    })
    
    output$cluster2 <- renderPlot({
      req(length(input_tf())>1)
      activity_test1 <- activity_data_cluster()[,3][[1]] # now we only plot the first tf input,
      # more plots could be generated later
      #activity_test1 <- create_activity_data("Arx")[,2][[1]]
      # add a activity column
      forebrain_with_activity <- mutate(forebrain_data, activity_1 = activity_test1) 
      
      ggplot(data = forebrain_with_activity, mapping = aes(x=UMAP1,y=UMAP2))+
        geom_point(aes(color = activity_1))+
        scale_color_gradientn(colors = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100)))+
        theme_bw()
      
    })
    
    
    
    # --------------------------------------Tab3: timeseries-------------------------------------------
    output$timeseries1 <- renderPlot({
      req(length(input_tf())>0)
      #tf_df <- as_tibble(rownames(activity))
      # tf_df is loaded at beginning using data_prep.R
      TF <- translate_tf(input_tf()[1],tf_df)
      req(TF)
      plot_timeseries(TF,cell_metadata_cortex, binary_activity)
 
    })
    output$timeseries2 <- renderPlot({
      req(length(input_tf())>1)
      TF <- translate_tf(input_tf()[2],tf_df)
      req(TF)
      plot_timeseries(TF,cell_metadata_cortex, binary_activity)
      
    })
    output$timeseries3 <- renderPlot({
      req(length(input_tf())>2)
      TF <- translate_tf(input_tf()[3],tf_df)
      req(TF)
      plot_timeseries(TF,cell_metadata_cortex, binary_activity)
      
    })
    
    
    
    
      
      
    
}
# Run the application 
shinyApp(ui = ui, server = server)
