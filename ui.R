ui <- fluidPage(
  introjsUI(),
  # Application title
  introBox(
    titlePanel("joint cortex and pons data app"),
    data.step = 1,
    data.intro = "This app displays transcription factor activity inference data from 
    a developmental timecourse of the mouse forebrain in pons in three main tabs including:"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # choose which datasets to analyze for the whole app
      actionButton("help", label = "See instructions"),
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
      # fileInput("file_tf", "Choose CSV File containing your tf list",
      #           accept = c(
      #             "text/csv",
      #             "text/comma-separated-values,text/plain",
      #             ".csv")
      # ),
      
      # 1. table and network graph of related TF and genes
      
      conditionalPanel(condition = "input.tabs == 'table and network'",
                      radioButtons("show", "Node display option",
                                    # use the names in the vector to display
                                    # use the character "joint_cortex" to match the path to import data
                                    choices = c("show all nodes" = "all",
                                                #"color by user's input tf list" = "pathway",
                                                "shrink graynodes" = "shrink",
                                                "neglect graynodes" = "neglect",
                                                
                                                "stop showing" = "stop"),
                                    selected = "stop"),#,
                       checkboxInput("show_pathway","color by user's input genes list",
                                     TRUE),
                       selectInput("input_pathway", "input your interested genes pathway",
                                   choices = data_cortex$unique_active_TFs_bare,
                                   multiple = TRUE,
                                   selected = c("Arx","Lef1"))
                       # fileInput("file_gene", "Choose CSV File containing your genes list",
                       #           accept = c(
                       #             "text/csv",
                       #             "text/comma-separated-values,text/plain",
                       #             ".csv")
                       # )
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
      
      introBox(
      # Update everything
      actionButton("update", label = "Update"),
      data.hint = "click me to update everything!",
      data.step = 2,
      data.intro = "click it to update everything! Do this after you changed your 
      tf input and options. Feel free to QUIT the intro first and update it to see the 
      table and plots",
      data.position = "right"
      ),
    ),
    mainPanel(
      tabsetPanel(
      
      
        tabPanel(
          title = "table and network",
          textOutput("general_desc"),
          introBox(
          
          dataTableOutput("table"),
          data.step = 3,
          data.intro = "Table and network tab: 
        A table of tf and its target gene with motifs and other information"
          ),
          introBox(
            data.intro = "Feel free to quit the intro now, click the 'show all nodes' button
            in the sidebar to see the cytoscape network graph, then we continue",
            data.step = 4
          ),
          
          
          textOutput("desc"),
          tags$style(type="text/css", "#desc {white-space: pre-wrap;}"),
          introBox(
          rcytoscapejsOutput("network", width = "1200px",height = "600px"),
          data.step = 5,
          data.intro = "a network graph visualization displaying detailed information with node color 
          and size: 
          Orange nodes are active transcription factors (tf genes that express their own tf);
          Purple nodes in the center are your input transcription factors;
          Green nodes are your input genes related to input tfs(purple nodes)
          ;  grey nodes are other genes."
          ),
          value = "table and network"
        ),
        
        
        
      
      tabPanel("heatmap and clustering",
               
               plotOutput("heatmap_cell"),
               downloadButton("download_hm_cell", "Heatmap by cell (Png)"),
               plotOutput("heatmap_cluster"),
               downloadButton("download_hm_cluster", "Heatmap by cluster (Png)"),
               imageOutput("color_hm_palette", width = "6in", height = "4in"),
               
               fluidRow(
                textOutput("cluster_UMAP_desc"),
                column(width = 6, plotOutput("cluster1",width = "5in", height = "5in"),
                       downloadButton("download_UMAP_1", "UMAP scatterplot 1 (Png)")),
               
                column(width = 6, plotOutput("cluster2", width = "5in",height = "5in"),
                       downloadButton("download_UMAP_2", "UMAP scatterplot 2 (Png)")),
              
               ),
               
               value = "heatmap and clustering"
      ),
      tabPanel("time series",
               textOutput("tf_timeseries_desc"),
               textOutput("timeseries_desc"),
               
               fluidRow(
               plotlyOutput("timeseries1"),
               downloadButton("download_ribbon_1", "Timeseries ribbon plot (Png)"),
               plotOutput("timeseries2"),
               imageOutput("timeseries_color"),
               #plotOutput("timeseries3"),
               #plotlyOutput("timeseries4")
               ),
               value = "time series"),
      id = "tabs"
    ))
  ),
  
  
)