server <- function(input, output, session) {
  # Dynamic UI, change the selectInput tf lists on display
  observeEvent(input$help,
               introjs(session, options = list("nextLabel"="Next",
                                               "prevLabel"="Did you forget something?",
                                               "skipLabel"="Don't be a quitter"),
                       events = list("oncomplete"=I('alert("Glad that is over, congrats!")')))
  )
  
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
    l$num_cell_plot <- input$num_cell_plot
    # l has following elements with same names for both options above:
    # l contains ...
    
    # We will use the same name attributes to retrieve data
    return (l)
    })
  
  #input_tf <- reactive(input_new()$tf)
  
  # -----------------------------Tab1:table and network------------------------------------------
    output$general_desc <- renderText({
      "This app designs for displaying transcription factor and gene data from mice brain (cortex & pons part) in various fancy ways by three main tabs;
                                             
       PROBLEM: There are some transcription factors from your input that may not have the corresponding data
       in the following tabs. (Sometimes you may not see the information of that transcription factor or the plot
       is not updated, etc. That is unfortunately because of the lack of data in the cell activity data in tab2,
       or the binary cell activity data in tab3.  
   "
      
    })
  
    output$table <- renderDataTable({
        # process data, filter the lines with our interested TF
      datatable(dplyr::filter(input_new()$TF_target_gene_info, TF %in% input_new()$tf))
    })
    
    
    output$desc <- renderText({
      text <- "\n  Orange nodes are active transcription factors (tf genes that express their own tf);
Purple nodes in the center are your input transcription factors;
Green nodes are your input genes related to input tfs(purple nodes);
grey nodes are other genes."
    })
    nodeData <- reactive(
      #input$show,
      {

      if(input$show_pathway){input_pathway <- input_new()$input_pathway}
      else{input_pathway <- c()}
      
      if(input$show == "all"){
        create_network(input_new()$tf, input_new()$TF_target_gene_info,
                       input_new()$unique_active_TFs_bare,
                       pathway_genes = input_pathway)$nodes
      }
      else if(input$show == "neglect"){
        create_network(input_new()$tf, input_new()$TF_target_gene_info,
                       input_new()$unique_active_TFs_bare,
                       pathway_genes = input_pathway)$nodes %>%
          filter(color!="lightgrey")
      }
      else if(input$show == "shrink"){
        create_network(input_new()$tf, input_new()$TF_target_gene_info,
                       input_new()$unique_active_TFs_bare,
                       pathway_genes = input_pathway,
                       shrink_gray = TRUE)$nodes
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
    output$color_hm_palette <- renderImage({
      
      expr = list(src = "www/timeseries_color.png",
           alt = "This is alternate text")
      
      
    },
    deleteFile = FALSE)
    
    hm_cell_plot <- reactive({
      req("Cell" %in% input_new()$method)
      plot_heatmap(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext,input_new()$cell_metadata,
                   cell_plot_num = input_new()$num_cell_plot)
      
    })
    
    hm_cluster_plot <- reactive({
      req("Cluster" %in% input_new()$method)
      plot_heatmap(input_new()$tf, "Cluster",input_new()$region, 
                   input_new()$TF_and_ext,input_new()$cell_metadata)
      
    })
    
    output$heatmap_cell <- renderPlot({
      hm_cell_plot()
    })
    
    output$download_hm_cell <- downloadHandler(filename = "heatmap_cell.png",
                                               contentType = "image/png",
                                               content = function(file){
                                                 ggsave(filename = file, plot = hm_cell_plot(),
                                                        width = 20, height = 25)
                                               })
    
    output$heatmap_cluster <- renderPlot({
      hm_cluster_plot()        
    })
    
    output$download_hm_cluster <- downloadHandler(filename = "heatmap_cluster.png",
                                               contentType = "image/png",
                                               content = function(file){
                                                 ggsave(filename = file, plot = hm_cluster_plot(),
                                                        width = 20, height = 25)
                                               })
  
    # The cluster scatterplot is always plot by cells, so we use an independent reactive
    # value for this plot
    activity_data_cluster <- reactive({
      # use the feature of feather data to read certain col to optimize speed
      create_activity_data(input_new()$tf, "Cell",input_new()$region, input_new()$TF_and_ext)
    })
    
    output$cluster_UMAP_desc <- renderText({
      text <- "Now we can only support two plots of your 
      first two transcription factor inputs."
    })
    
    Umap_plot_1 <- reactive({
      req(length(input_new()$tf)>0)
      plot_UMAP(tf_number = 1,input_new()$cell_metadata, activity_data_cluster())
    })
    Umap_plot_2 <- reactive({
      req(length(input_new()$tf)>1)
      plot_UMAP(tf_number = 2,input_new()$cell_metadata, activity_data_cluster())
    })
    output$cluster1 <- renderPlot({
      Umap_plot_1()
      
    })
    
    output$cluster2 <- renderPlot({
      Umap_plot_2()
    })
    
    output$download_UMAP_1 <- downloadHandler(filename = "UMAP1.png",
                                                  contentType = "image/png",
                                                  content = function(file){
                                                    ggsave(filename = file, plot = Umap_plot_1(),
                                                           width = 20, height = 20)
                                                  })
    output$download_UMAP_2 <- downloadHandler(filename = "UMAP2.png",
                                              contentType = "image/png",
                                              content = function(file){
                                                ggsave(filename = file, plot = Umap_plot_2(),
                                                       width = 20, height = 20)
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
    # tf_desc_timeseries <- reactive({
    #   tf_nexist_string <- ""
    #   for(tf_n in input_new()$tfs_not_exist_timeseries){
    #     tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
    #   }
    #   text <- glue("We do not have these followning tfs in this tab: {tf_nexist_string}")
    #   
    #   tf_nexist <- ""
    #   for(tf in input_new()$tf){
    #     if (tf %in% input_new()$tfs_not_exist_timeseries){
    #       tf_nexist <- paste(tf_nexist,tf,sep = " ")
    #     }
    #   }
    # 
    #   if(tf_nexist == ""){
    #     text <- "Good! All of your input tfs exist in our timeseries activity datasets!"
    #   }
    #   else{
    #     tf_nexist_string <- ""
    #     for(tf_n in input_new()$tfs_not_exist_timeseries){
    #       tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
    #     }
    #     text <- glue('Those tfs in your input list does not not exist in our
    #                timeseries datasets: {tf_nexist}.
    #                We do not have these followning tfs in this tab: {tf_nexist_string}')
    #   }
    # })
    
    output$tf_timeseries_desc <- renderText({
      # tf_desc_timeseries()
      tf_nexist_string <- ""
      for(tf_n in input_new()$tfs_not_exist_timeseries){
        tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
      }
      text <- glue("We do not have these followning tfs in this tab: {tf_nexist_string}")
      
      
    })
    
    
    output$timeseries_desc <- renderText({
      text <- "Click option: You may double click the color palatte of cell types at the right side to 
      display that cell type ONLY; you could also click on one cell type to eliminate that in the
      plot at left.
      Mouse over the white vertical line on the plot to see the cell types. 
      We only support four plots of your first four tfs input for now."
    
    })
    
    # we must transform the TF format, from raw form (Arx) to (Arx_extended (21g)) to fetch
    # information
    TF_transformed <- reactive({
      translate_tf(input_new()$tf,input_new()$binary_active_TFs)
      })
    
    # ggplotly_list_plot <- reactive({
    #   req(TF_transformed())
    #   # binary_active_TFs is loaded at beginning by data_prep.R
    #   plot_list <- lapply(TF_transformed(), plot_timeseries, cell_metadata = data_cortex$timeseries_input_meta, 
    #                       activity = data_cortex$binary_activity, make_plotly = TRUE)
    #   # produce a list of ggplotly plots
    #   subplot(plot_list, nrows = 2, margin = 0.04, heights = c(0.6, 0.4), shareX = TRUE, shareY = FALSE)
    #   
    # })
    
    
    ggplot_list_plot <- reactive({
      req(TF_transformed())
      plot_list <- lapply(TF_transformed(), plot_timeseries, cell_metadata = input_new()$timeseries_input_meta, 
                          activity = input_new()$binary_activity, make_plotly = FALSE, show_legend = FALSE)
      plot_grid(plotlist = plot_list)
    })
    
    output$timeseries1 <- renderPlotly({ # a plotly list
      req(length(input_new()$tf)>0)
      plot_timeseries(TF_transformed()[1][1], input_new()$timeseries_input_meta, input_new()$binary_activity,make_plotly = TRUE)
     })
    
    output$download_ribbon_1 <- downloadHandler(filename = "timeseries_ribbon.png",
                                                contentType = "image/png",
                                                content = function(file){
                                                  ggsave(filename = file, plot = ggplot_list_plot(),
                                                         width = 20, height = 15)
                                                })
    
    output$timeseries2 <- renderPlot({ # a ggplot list
      ggplot_list_plot()
      
    })
    
    output$timeseries_color <- renderImage({
      
      list(src = "www/timeseries_color.png",
           alt = "This is alternate text")
      
    },
    deleteFile = FALSE)
    
    
    
    
}

