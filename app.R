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


#------------------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("joint cortex app"),
    
    sidebarLayout(
        sidebarPanel(
            # choose which datasets to analyze for the whole app
            radioButtons("region", "Brain region",
                         # use the names in the vector
                         choices = c("Forebrain" = "joint_cortex",
                                     "Pons" = "joint_pons"),
                         selected = "Forebrain"),
            
            # 1. table and network graph of related TF and genes
            conditionalPanel(condition = "input.tabs == 'table and network'",
                             
                             selectInput(inputId = "TF",
                                         label = "transcription factor",
                                         choices = unique_TF),
                             actionButton("update_network", label = "Update")),
            # 2. heatmap and clustering
            conditionalPanel(condition = "input.tabs == 'heatmap and clustering'",
                             
                             selectInput(inputId = "transcription factor",
                                         label = "TF",
                                         choices = TF_active),
                             actionButton("update_heatmap", label = "Update")),
            # 3. table and network graph of related TF and genes
            conditionalPanel(condition = "input.tabs == 'time series'",
                             
                             selectInput(inputId = "transcription factor",
                                         label = "TF",
                                         choices = TF_active),
                             actionButton("update_timeseries", label = "Update"))
        ),
        mainPanel(tabsetPanel(
            
            tabPanel(title = "table and network",
                     tableOutput("table"),
                     value = "table and network"
            ),
            
            tabPanel("heatmap and clustering",
                     plotOutput("heatmap"),
                     plotOutput("clustering"),
                     
                     value = "heatmap and clustering"
            ),
            
            tabPanel("time series",
                     value = "time series"),
            id = "tabs"
        ))
    ),
    
    
)

server <- function(input, output, session) {
    
    eventReactive(input$region,{
        if(input$region == "joint_cortex"){
            # datasets
            forebrain_data <- read_tsv("data/joint_cortex/Forebrain_join.2D.tsv")
            #---------------------------------------------------------------------------
            # These datasets describe TF and genes that are target of TFs doesn't have ext suffix
            TF_target_gene <- as_tibble(read_rds("data/joint_cortex/joint_cortex.regulon_target_info.Rds"))
            TF_target_gene <- select(TF_target_gene,-logo)
            unique_TF <- unique(TF_target_gene[["TF"]])
            
            #---------------------------------------------------------------------------
            # TF has ext and weight suffix in these datas
            # a vector that contains all TFs(ext or regular) in the activity data(by cluster/cell)
            TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
            activity_cluster <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cluster.feather")
            activity_cell <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather")
            
            TF_and_ext <- TF_active %>% 
                rename(name = value) %>%
                mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the gram
                mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
                mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
                mutate(type = str_replace(TF_type, "_ext", "")) %>%
                select(name, type, ext)
            # metadata
            cluster_info <- read_tsv("data/joint_cortex/metadata_20190716.tsv")
        }
        
    })
    
    output$table <- renderTable({ # display a table
        # process data
        filter(TF_target_gene, TF == input$TF)
    })
    
    
    
    
    
    
}
# Run the application 
shinyApp(ui = ui, server = server)
