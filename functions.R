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
identify_tf <- function(TF_name_activity_tibble){
  TF_name_activity_tibble %>% 
    rename(name = value) %>%
    mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the gram
    mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
    mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
    mutate(type = str_replace(TF_type, "_ext", "")) %>%
    select(name, type, ext)
}
# boolean functions that determine the identity of TF data in the activity datasets
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
# read the corresponding data by tf's identity
tf_regular <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext==TF)[[1,1]]
}
tf_ext <- function(TF, TF_ext_data){
  filter(TF_ext_data, type==TF & ext=="ext")[[1,1]]
}

# NODE: TF_and_ext is a dataframe (loaded already) that created in order to identify whether the TF data
# is a regular TF (with high confidence annotation) or ext type(with lower confidence)
# Dogma: if we have regular TF type, we use that data; if we only have ext data, use ext
# if we have no data related to this tf, we will either give an error or do nothing
create_activity_data <- function(tf){ 
  # use the feature of feather data to read certain col to optimize speed
  cell_col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather",
                           "Cell")
  # add certain tf activity data to the Cell column
  activity_cell <- cell_col
  for(TF in tf){ # tf is input tf list, could contain many tfs
    tf_to_read <- TF
    if(has_regular(TF, TF_ext_data)){
      tf_to_read <- tf_regular(TF,TF_ext_data) # a helper to read the corresponding data
    }
    else if(has_ext(TF, TF_ext_data)){
      tf_to_read <- tf_ext(TF,TF_ext_data)
    }
    else{
      next # means we don't have that data, we jump over it and do thing
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


#' Title
#'
#' @param tf 
#'
#' @return
#' @export
#'
#' @examples
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
                             id %in% unique_TF ~ "#D6604D",
                             id %in% mutual_target ~ "#4fafc6",
                             TRUE ~ "lightgrey"))
  return(list(nodes = nodeData,
              edges = edgeData
              ))
}




# make_cell_metadata
create_cell_metadata <- function(metadata_part){
  metadata_part %>% 
    select(Age = orig.ident, Cell, Cluster = ID_20190730_with_blacklist_and_refined) %>% 
    # In this case, we remove the "prefix" of the Cluster column, so that we are
    # simply left with the abbreviation representing the cell type, so that 
    # we can link the cells of the same cell type across ages
    separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
    mutate(Age = factor(Age, levels = c("Forebrain E12.5",
                                        "Forebrain E15.5",
                                        "Forebrain P0",
                                        "Forebrain P3",
                                        "Forebrain P6"))) %>% 
    arrange(Cell)
  
}
#cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)


# since we have two datasets with different tf name, with ext and weights or not
# and those tf names are essential for matching data, this function use a clean tf
# and return a tf with ext and weight type
translate_tf <- function(tf, tf_dataframe){
  identify_tf(tf_dataframe) %>%
    filter(type==tf) %>%
    .[[1,1]]
}
tf_df <- as_tibble(rownames(activity))
  

TF <- "Arx_extended (21g)" # input_tf

#activity <- readRDS("data/joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")

#plot_timeseries(TF,cell_metadata_cortex,activity)

plot_timeseries <- function(TF,cell_metadata, activity){
  activity <- activity[TF, ] %>% 
    {data.frame("TF" = .)} %>% 
    tibble::rownames_to_column(var = "Cell") %>% # the original activity vector has names
    arrange(Cell)
  
  if(!all(cell_metadata$Cell == activity$Cell)) return (-1)
  # Add the TF activity to the new dataframe
  ribbon_df <- cell_metadata
  ribbon_df$TF <- activity$TF
  
  ribbon_df <- ribbon_df %>% 
    filter(!grepl("BLACKLIST", Cluster))
  ribbon_df_celltype_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    # Total cells at each age
    mutate(total = n()) %>% 
    group_by(Age, Cluster) %>%
    # Proportion of TF+ cells per cluster, per age
    mutate(frac = sum(TF > 0) / total) %>% 
    distinct(Age, Cluster, frac) %>% 
    ungroup()
  
  ribbon_df_cum_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    summarize(cumfrac = sum(TF > 0) / n()) %>% 
    ungroup()
  
  timepoints2 <- ribbon_df$Age
  clusters <- ribbon_df$Cluster
  
  df = data.frame(cluster = rep(unique(clusters), length(unique(timepoints2))),
                  stage = do.call(c, lapply(as.character(unique(timepoints2)), rep, times = length(unique(clusters)))))
  
  df$ranking = match(df$cluster, names(colours))
  df = df[order(df$stage, df$ranking),]
  
  df <- left_join(df, select(ribbon_df_celltype_frac, cluster = Cluster, stage = Age, frac)) %>% 
    mutate(frac = replace_na(frac, 0)) %>% 
    left_join(select(ribbon_df_cum_frac, stage = Age, cumfrac))
  
  df$xpos = match(df$stage, unique(timepoints2))
  
  df %>%
    ggplot(aes(x = xpos, y = frac, fill = cluster)) +
    geom_area(stat = "identity") +
    scale_fill_manual(values = colour_palette, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(unique(df$stage)),
                       labels = c("E12.5", "E15.5", "P0", "P3", "P6"),
                       limits = c(1, length(unique(df$stage)))) +
    labs(x = "age", y = "Binary activity", title = TF) +
    guides(fill = guide_legend(ncol = 5)) +
    theme(legend.position = "bottom")
  
}





TF_tbl1 <- as_tibble(TF_active)
TF_tb2<- as_tibble(tf_list)

missing_tf <- TF_tbl1 %>%
  mutate(exist = if_else(value %in% tf_list, 1, 0)) %>%
  filter(exist==0)

mutual_tf <- TF_tbl1 %>%
  mutate(exist = if_else(value %in% tf_list, 1, 0)) %>%
  filter(exist==1)
  
  
  
  
  
  

