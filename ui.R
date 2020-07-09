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
               plotlyOutput("cluster1"),
               plotlyOutput("cluster2"),
               
               value = "heatmap and clustering"
      ),
      tabPanel("time series",
               textOutput("tf_timeseries_desc"),
               textOutput("timeseries_desc"),
               plotlyOutput("timeseries1"),
               plotlyOutput("timeseries2"),
               plotlyOutput("timeseries3"),
               plotlyOutput("timeseries4"),
               value = "time series"),
      id = "tabs"
    ))
  ),
  
  
)