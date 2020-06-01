plotByValue <- function(tf, cluster_data){
  
  # because the data name is not good(contain space), 
  # we firstly make a vector of that then call it, pick the first one
  TF_activity <- select(cluster_data, contains(tf))[[1]] # a vector, not a list
  
  #arrange(cluster_data, desc(tf_col)) # sort
  cluster_name <- cluster_data[[1]]
  
  ggplot(data = cluster_data)+
    geom_point(mapping = aes(x = cluster_name, y = TF_activity, color = cluster_name))+
    theme(axis.text.x= element_text(angle = 90, hjust = 1))
  
}
# generate a tibble that has two columns indicating whether the tf has ext type,
# the ext column labels whether that data is ext type
TF_and_ext <- function(TF_name_activity_tibble){
  TF_name_activity_tibble %>% 
    rename(name = value) %>%
    mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the gram
    mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
    mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
    mutate(type = str_replace(TF_type, "_ext", "")) %>%
    select(name, type, ext)
}
# three booleans that determine the quality of TF data in the activity datasets
has_regular <- function(TF, TF_ext_data){
  has_regular <- filter(TF_ext_data, type==TF & ext==TF)
  nrow(has_regular)!=0 #boolean
}
has_ext <- function(TF, TF_ext_data){
  is_ext <- filter(TF_ext_data, type==TF & ext=="ext")
  nrow(is_ext)!=0
}
tf_exist <- function(TF){
  for(tf in TF){
    if(has_regular(tf, TF_ext_data) || has_ext(tf, TF_ext_data)){
    }
    else{
      return (tf)
    }
  }
  return (TRUE)
  #!has_regular(TF, TF_ext_data) && !has_ext(TF, TF_ext_data)
}

tf_regular <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext==TF)[[1,1]]
}
tf_ext <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext=="ext")[[1,1]]
}

create_activity_data <- function(tf){ 
  TF_activity_cell_tibble <- as_tibble(TF_active)
  TF_ext_data <- TF_and_ext(TF_activity_cell_tibble)
  # use the feature of feather data to read certain col to optimize speed
  cell_col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather",
                           "Cell")
  # note that the first
  activity_cell <- cell_col
  for(TF in tf){
    tf_to_read <- TF
    if(has_regular(TF, TF_ext_data)){
      tf_to_read <- tf_regular(TF,TF_ext_data)
    }
    else if(has_ext(TF, TF_ext_data)){
      tf_to_read <- tf_ext(TF,TF_ext_data)
    }
    else{
      next
    }
    col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather",
                        tf_to_read)
    activity_cell <- add_column(activity_cell,col)
  }
  activity_cell %>%
    select(Cell, everything()) # move cell col to start
}


#--------------------------------Tab 2: heatmap----------------------------------
# This function takes a colour palette as input,
# and creates the data formats needed to annotate
# a pheatmap with some colours

makePheatmapAnno <- function(palette, column) {
  
  palette <- palette[unique(names(palette))]
  
  # Make dataframe, retrieve the data frame
  anno_row <- data.frame(cluster = names(palette))
  names(anno_row) <- column
  rownames(anno_row) <- anno_row[[1]]
  
  # Make list containing colours
  side_colors <- list(cluster = palette)
  names(side_colors) <- column
  
  return(list(anno_row = anno_row,
              side_colors = side_colors))
  
}

plot_heatmap <- function(metadata, TF_active, activity_cluster){
  
  colour_palette <- metadata %>% 
    # use gsub to change all contents in Cluster (cluster name format)
    mutate(Cluster = gsub("_", " ", Cluster)) %>% 
    # Get two columns
    select(Cluster, Colour) %>% 
    # Convert to vector of colours, where the first column gives the names
    # and the second column is converted to the values
    deframe() # VECTOR , not data frame
  
  activity_cluster %>%
    select(c("Cluster", TF_active[1])) %>% 
    tibble::column_to_rownames(var = "Cluster") # make that column name as row name
  
  # A helper function to prepare a dataframe to annotate the heatmap with colours
  hm_anno <- makePheatmapAnno(colour_palette, "Cluster")
  
  pheatmap::pheatmap(t(act),
                     scale = "none",
                     border_color = NA,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = "My 1st heatmap",
                     annotation_col = hm_anno$anno_row,
                     # change the default color annotation
                     annotation_colors = hm_anno$side_colors, 
                     annotation_legend = TRUE,
                     cellwidth = 10,
                     cellheight = 10)
} 


# -------------------------------------cytoscape----------------------------------------------
# function to create network


create_network <- function(tf){ 
  # takes a vector input that contains user selected TFs
  # good to visualize correlations among multiple TFs
  TF_interest <- filter(TF_target_gene, TF %in% tf)[["TF"]]
  gene_target <- filter(TF_target_gene, TF %in% tf)[["gene"]]
  
  source <- TF_interest
  target <- gene_target
  
  id <- c(TF_interest, gene_target)
  name <- id
  nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
  edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
  
  mutual_target <- edgeData %>% 
    # a character vector that indicates the nodes that are target of multiple selected TFs
    count(target) %>%
    filter(n > 1 & !target %in% tf ) %>%
    .[[1]]
  
  nodeData <- nodeData %>%
    mutate(color = case_when(id %in% tf ~ "#9d4097",
                             id %in% mutual_target ~ "#4fafc6",
                             TRUE ~ "lightgrey"))
  return(list(nodes = nodeData,
              edges = edgeData
              ))
}
tf <- c("Arx","Lef1")

nodeData <- create_network(tf)$nodes
edgeData <- create_network(tf)$edges
network <- createCytoscapeJsNetwork(nodeData, edgeData)
rcytoscapejs2(network$nodes, network$edges)



network_data <- create_network(c("Arx","Lef1"))

nodeData <- network_data$nodes
edgeData <- network_data$edges

unique_edges <- edgeData %>%
    count(target) %>%
    filter(n > 1 & !target %in% c("Arx","Dlx1") ) %>%
  .[[1]]

nodeData <- nodeData %>%
  mutate(color = case_when(id %in% c("Arx","Lef1") ~ "#9d4097",
                           id %in% unique_edges ~ "#4fafc6",
                           TRUE ~ "lightgrey"))

network <- createCytoscapeJsNetwork(nodeData, edgeData)

rcytoscapejs2(network$nodes, network$edges)


