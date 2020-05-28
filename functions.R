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
has_nothing <- function(TF, TF_ext_data){
  !has_regular(TF, TF_ext_data) & !has_ext(TF, TF_ext_data)
}

tf_regular <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext==TF)[[1,1]]
}
tf_ext <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext=="ext")[[1,1]]
}

activity_dataset <- function(tf){ 
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

plot_heatmap <- function(cluster_info, TF_active, activity_cluster){
  
  colour_palette <- cluster_info %>% 
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


#plot_heatmap(cluster_info,TF_active,activity_cluster)



