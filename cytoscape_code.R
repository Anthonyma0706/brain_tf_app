# Florian Wunnemann

# function to create network
create_cytoscape_network <- function(regulon_selected,
                                     tf_list){
  ## Create list of TF and all targets as nodes
  target_genes <- regulon_selected[[1]] # a vector type, use vector type to subset information
  target_genes_tf <- subset(target_genes, target_genes %in% tf_list)
  target_genes_not_tf<- target_genes[! target_genes %in% tf_list] # dataframe/list type
  
  ## Merge regulon TF and different annotated targets
  id <- c(names(regulon_selected),target_genes_tf,target_genes_not_tf)
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  nodeData$color <- c("#9d4097",
                      replicate(length(target_genes_tf),"#4fafc6"),
                      replicate(length(target_genes_not_tf),"lightgrey"))
  
  nodeData$shape <- c("octagon",
                      replicate(length(target_genes_tf),"octagon"),
                      replicate(length(target_genes_not_tf),"ellipse"))
  
  nodeData$nodeLabelColor <- c("white",
                               replicate(length(target_genes_tf),"white"),
                               replicate(length(target_genes_not_tf),"black"))
  
  ## Format interaction between TF and targets
  source <- replicate(length(regulon_selected[[1]]),names(regulon_selected))
  target <- regulon_selected[[1]]
  edgeData <- data.frame(source, target, stringsAsFactors=FALSE)
  
  return(list(nodes = nodeData,
              edges = edgeData))
}

ui
library(rcytoscapejs)
rcytoscapejsOutput("Rcyotscape_selected_network",
                               width = "1200px",
                               height= "600px"),


server
library(rcytoscapejs)
## Check size of regulon and only make a Cytoscape animation if < 100
      
      ## Cytoscape function to render the network
      output$Rcyotscape_selected_network <-  renderRcytoscapejs({
        ## Read in List of transcriptions in data
        tf_table <- read.delim("/home/florian/data/grns/inputTFs_for_GRNBoost.txt",
                             sep="\t",
                             col.names = FALSE,
                             quote="")
        tf_list <- tf_table[,1]
        
        ## Check length of regulon
        if(length(regulon_selection()[[1]]) < 100){
        
        regulon_selected <- regulon_selection()
        
        cyto_data <- create_cytoscape_network(regulon_selected = regulon_selected,
                                 tf_list = tf_list)
        
        Nodedata <- cyto_data$nodes
        EdgeData <- cyto_data$edges
        
        network <- createCytoscapeJsNetwork(Nodedata, EdgeData)
        
        rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE,
                     highlightConnectedNodes = FALSE)
        } else {
          
          regulon_selected <- regulon_selection()
        
          ## Get the top 100 targets based on a specific score
          ## Only show top 100 genes based on NES score and nMotifs
          top100 <- regulonTargetsInfo_selected() %>%
            arrange(desc(nMotifs)) %>%
            top_n(n=100,
                  wt = NES)
          
          top100 <- top100$gene[1:100]
          
          regulon_selected_subset <- list(regulon_selected[[1]][regulon_selected[[1]] %in% top100])
          names(regulon_selected_subset) <- names(regulon_selected)
           
          assign("top100",top100, envir = globalenv())   
          assign("regulon_selected",regulon_selected, envir = globalenv())  
          assign("regulon_selected_subset",regulon_selected_subset, envir = globalenv())
                             
          cyto_data <- create_cytoscape_network(regulon_selected = regulon_selected_subset,
                                                tf_list = tf_list)
          
          assign("file_markers",cyto_data, envir = globalenv()) 
          
          Nodedata <- cyto_data$nodes
          EdgeData <- cyto_data$edges
          
          network <- createCytoscapeJsNetwork(Nodedata, EdgeData)
          
          rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE,
                       highlightConnectedNodes = FALSE)
          }
        })
      
      output$sel_regulon_table <- renderDT({
        selected_regulon_targets <- regulonTargetsInfo_selected()
        return(selected_regulon_targets)
      }, escape = FALSE
      )
      
      
      
      
      
      
library(rcytoscapejs2)
id <- c("Jerry", "Elaine", "Kramer", "George")
name <- id
nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
source <- c("Jerry", "Jerry", "Jerry", "Elaine", "Elaine", "Kramer", "Kramer", "Kramer", "George")
target <- c("Elaine", "Kramer", "George", "Jerry", "Kramer", "Jerry", "Elaine", "George", "Jerry")
edgeData <- data.frame(source, target, stringsAsFactors=FALSE) 
edgeData <- edgeData[1:5,]


network <- createCytoscapeJsNetwork(nodeData, edgeData)
rcytoscapejs2(network$nodes, network$edges)
library(tidyverse)
TF_target_gene <- as_tibble(read_rds("data/joint_cortex/joint_cortex.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF <- unique(TF_target_gene[["TF"]])

TF_interest1 <- filter(TF_target_gene, TF %in% "Arx")[["TF"]]
gene_target1 <- filter(TF_target_gene, TF %in% "Arx")[["gene"]]

source <- TF_interest
target <- gene_target

id <- c(TF_interest, gene_target)
name <- id
nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
edgeData <- data.frame(source, target, stringsAsFactors = FALSE)


create_network <- function(tf){
  TF_interest <- filter(TF_target_gene, TF %in% tf)[["TF"]]
  gene_target <- filter(TF_target_gene, TF %in% tf)[["gene"]]
  
  source <- TF_interest
  target <- gene_target
  
  id <- c(TF_interest, gene_target)
  name <- id
  nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
  edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
  return(list(nodes = nodeData,
              edges = edgeData))
  }
nodeData <- create_network("Arx")$nodes
nodeData$color <- c("#9d4097",
                    replicate(length(target_genes_tf),"#4fafc6"),
                    replicate(length(target_genes_not_tf),"lightgrey"))
nodeData$color <- "#9d4097" %>%
  mutate()

edgeData <- create_network("Arx")$edges


network <- createCytoscapeJsNetwork(nodeData, edgeData)

rcytoscapejs2(network$nodes, network$edges)

create_network <- function(tf){
  TF_interest <- filter(TF_target_gene, TF %in% tf)[["TF"]]
  gene_target <- filter(TF_target_gene, TF %in% tf)[["gene"]]
  
  source <- TF_interest
  target <- gene_target
  
  id <- c(TF_interest, gene_target)
  name <- id
  
  nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
  edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
  return(list(nodes = nodeData,
              edges = edgeData))
}

nodeData <- create_network(c("Arx","Dlx1"))$nodes
edgeData <- create_network(c("Arx","Dlx1"))$edges



      