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
    l$method <- input$method
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
      text <- "Orange nodes are active transcription factors(tf genes that express their own tf); " %>%
        paste("Purple nodes in the center are your input transcription factors; ") %>%
        paste("Green nodes are your input genes related to input tfs(purple nodes); ") %>%
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
      req("Cell" %in% input_new()$method)
      plot_heatmap(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext,input_new()$cell_metadata,
                   cell_plot_num = input$num_cell_plot)
    })
    
    output$heatmap_cluster <- renderPlot({
      req("Cluster" %in% input_new()$method)
      plot_heatmap(input_new()$tf, "Cluster",input_new()$region, input_new()$TF_and_ext,input_new()$cell_metadata)
                   
    })
  
    # The cluster scatterplot is always plot by cells, so we use an independent reactive
    # value for this plot
    activity_data_cluster <- reactive({
      # use the feature of feather data to read certain col to optimize speed
      create_activity_data(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext)
    })
    
    output$cluster1 <- renderPlotly({
      req(length(input_new()$tf)>0)
      p <- plot_UMAP(tf_number = 1,input_new()$cell_metadata, activity_data_cluster())
      ggplotly(p)
    })
    
    output$cluster2 <- renderPlotly({
      req(length(input_new()$tf)>1)
      p <- plot_UMAP(tf_number = 2,input_new()$cell_metadata, activity_data_cluster())
      ggplotly(p)
    })
    
    
    
    # --------------------------------------Tab3: timeseries-------------------------------------------
    # tf_nexist_data <- reactive({
    #   tf_nexist <- ""
    #   for(tf in input_new()$tf){
    #     if (tf %in% input_new()$tfs_not_exist_timeseries){
    #       tf_nexist <- paste(tf_nexist,tf,sep = " ")
    #     }
    #   }
    # })
    tf_desc_timeseries <- reactive({
      tf_nexist_string <- ""
      for(tf_n in input_new()$tfs_not_exist_timeseries){
        tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
      }
      text <- glue("We do not have these followning tfs in this tab: {tf_nexist_string}")
      
      # tf_nexist <- ""
      # for(tf in input_new()$tf){
      #   if (tf %in% input_new()$tfs_not_exist_timeseries){
      #     tf_nexist <- paste(tf_nexist,tf,sep = " ")
      #   }
      # }
      # 
      # if(tf_nexist == ""){
      #   text <- "Good! All of your input tfs exist in our timeseries activity datasets!"
      # }
      # else{
      #   # tf_nexist_string <- ""
      #   # for(tf_n in input_new()$tfs_not_exist_timeseries){
      #   #   tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
      #   # }
      #   text <- glue('Those tfs in your input list does not not exist in our 
      #              timeseries datasets: {tf_nexist}.
      #              We do not have these followning tfs in this tab: {tf_nexist_string}')
      # }
    })
    
    output$tf_timeseries_desc <- renderText({
      tf_desc_timeseries()
      
    })
    
    
    output$timeseries_desc <- renderText({
      text <- "Click option: You may double click the color palatte of cell types at the right side to 
      display that cell type ONLY; you could also click on one cell type to eliminate that in the
      plot at left.
      Mouse over the plot to see the cell types"
    
    })
    
    
    
    output$timeseries1 <- renderPlotly({
      req(length(input_new()$tf)>0)
      # binary_active_TFs is loaded at beginning by data_prep.R
      TF <- translate_tf(input_new()$tf[1],input_new()$binary_active_TFs)
      req(TF)
      ggplotly(plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity))
      
    })
    output$timeseries2 <- renderPlotly({
      req(length(input_new()$tf)>1)
      TF <- translate_tf(input_new()$tf[2],input_new()$binary_active_TFs)
      req(TF)
      ggplotly(plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity))
      
    })
    output$timeseries3 <- renderPlotly({
      req(length(input_new()$tf)>2)
      TF <- translate_tf(input_new()$tf[3],input_new()$binary_active_TFs)
      req(TF)
      ggplotly(plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity))
    })
    
    output$timeseries4 <- renderPlotly({
      req(length(input_new()$tf)>3)
      TF <- translate_tf(input_new()$tf[4],input_new()$binary_active_TFs)
      req(TF)
      ggplotly(plot_timeseries(TF,input_new()$timeseries_input_meta, input_new()$binary_activity))
    })
    
    
    
}

