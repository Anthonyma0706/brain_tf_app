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

source("functions.R")

# datasets
forebrain_data <- read_tsv("data/joint_cortex/Forebrain_join.2D.tsv")
#---------------------------------------------------------------------------
# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene <- as_tibble(read_rds("data/joint_cortex/joint_cortex.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF <- unique(TF_target_gene[["TF"]])

#---------------------------------------------------------------------------
# TF has ext and weight suffix in these datas
# a vector that contains all TFs(ext or regular) in the activity data(by cluster/cell)
TF_active <- read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds")

# import feather files later based on input$TF
#activity_cluster <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cluster.feather")
#activity_cell <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather")

# metadata
metadata <- read_tsv("data/joint_cortex/metadata_20190716.tsv")


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
            
            selectInput(inputId = "TF",
                        label = "transcription factor",
                        choices = unique_TF,
                        multiple = TRUE,
                        selected = c("Arx","Lef1")),
            
            # 1. table and network graph of related TF and genes
            conditionalPanel(condition = "input.tabs == 'table and network'",
              
                             actionButton("update_network", label = "Update")),
            # 2. heatmap and clustering
            conditionalPanel(condition = "input.tabs == 'heatmap and clustering'",
                             numericInput(inputId = "num_cell", label = "number of cells to visualize",
                                          value = 20),
                             
                             actionButton("update_heatmap", label = "Update")),
            # 3. time series plot
            conditionalPanel(condition = "input.tabs == 'time series'",
                            
                             actionButton("update_timeseries", label = "Update"))
        ),
        mainPanel(tabsetPanel(
            
            tabPanel(title = "table and network",
                     dataTableOutput("table"),
                     rcytoscapejsOutput("network", width = "1200px",height = "600px"),
                     value = "table and network"
            ),
            tabPanel("heatmap and clustering",
                     plotOutput("heatmap"),
                     plotOutput("cluster"),
                     textOutput("tf_not_exist"),
                     
                     value = "heatmap and clustering"
            ),
            tabPanel("time series",
                     value = "time series"),
            id = "tabs"
        ))
    ),
    
    
)

server <- function(input, output) {

    output$table <- renderDataTable({
        # process data, filter the lines with our interested TF
      datatable(dplyr::filter(TF_target_gene, TF %in% input$TF))
      #datatable(activity_data())
    })
    output$network <- renderRcytoscapejs({
      
      nodeData <- create_network(input$TF)$nodes
      edgeData <- create_network(input$TF)$edges
      network <- createCytoscapeJsNetwork(nodeData, edgeData)
      rcytoscapejs2(network$nodes, network$edges)
      
    })
    
    
    
    activity_data <- reactive({
      # use the feature of feather data to read certain col to optimize speed
      create_activity_data(input$TF)
    })
    
    
    output$heatmap <- renderPlot({
      act <- activity_data() %>%
        sample_n(input$num_cell) %>%
        tibble::column_to_rownames(var = "Cell")
      pheatmap::pheatmap(t(act),
                         scale = "none",
                         border_color = NA,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         main = "heatmap",
                         #annotation_col = hm_anno$anno_row,
                         # change the default color annotation
                         #annotation_colors = hm_anno$side_colors, 
                         annotation_legend = TRUE,
                         cellwidth = 10,
                         cellheight = 10)
      
    })
    
    output$cluster <- renderPlot({
      activity_test1 <- create_activity_data("Arx")[,2][[1]]
      # add a activity column
      forebrain_with_activity <- mutate(forebrain_data, activity_1 = activity_test1) 
      
      ggplot(data = forebrain_with_activity, mapping = aes(x=UMAP1,y=UMAP2))+
        geom_point(aes(color = activity_1))+
        theme_bw()
      
    })
      
      
    
}
# Run the application 
shinyApp(ui = ui, server = server)
