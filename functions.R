load("data/joint_cortex/cortex_prep.Rda")
# --------------------------Tab1-cytoscape network visualization----------------------------------------------
# function to create network


#' Create rcytoscape network data
#' 
#' Takes a vector input that contains user selected TFs and output a list of nodeData and edgeData
#' which will be used for createCytoscapeJsNetwork
#' A good to visualize correlations among multiple TFs
#' 
#' @param tf one single tf name character
#'
#' @return a list of nodeData and edgeData that are required for generating a rcytoscapejs network object
#' 
#' @examples 
#' TF <- c("Arx","Lef1")
#' nodeData <- create_network(TF)$nodes
#' edgeData <- create_network(TF)$edges
#' network <- createCytoscapeJsNetwork(nodeData, edgeData)
#' rcytoscapejs2(network$nodes, network$edges)
#' 
create_network <- function(tf){ 
  TF_interest <- filter(TF_target_gene, TF %in% tf)[["TF"]]
  gene_target <- filter(TF_target_gene, TF %in% tf)[["gene"]]
  
  source <- TF_interest
  target <- gene_target
  
  id <- c(TF_interest, gene_target)
  name <- id
  nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
  edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
  
  #unique_TF <- unique(TF_target_gene[["TF"]])
  
  mutual_target <- edgeData %>% 
    # a character vector that indicates the nodes that are target of multiple selected TFs
    count(target) %>%
    filter(n > 1 & !target %in% tf ) %>%
    .[[1]]
  
  nodeData <- nodeData %>%
    mutate(color = case_when(id %in% tf ~ "#9d4097",
                             id %in% unique_TF ~ "#D6604D",
                             id %in% mutual_target ~ "#4fafc6",
                             TRUE ~ "lightgrey")) %>%
    mutate(height = case_when(id %in% tf ~ "100",
                           TRUE ~ "70")) %>%
    mutate(width = case_when(id %in% tf ~ "100",
                            TRUE ~ "70"))
  return(list(nodes = nodeData,
              edges = edgeData
  ))
}
# ------------------------------------------------------------------------------------
#' Identify transcription factor data type
#' 
#' Generate a tibble that has two columns indicating whether the tf has ext type,
#' the ext column labels whether that data is ext type
#'
#' @param TF_name_activity_tibble a dataframe/tibble that saves all names of transcription 
#' factor with suffix ext and weights
#' E2f1_extended (133g) and E2f1 (122g) are two examples of tf data detected
#' tf with no '_extended' attached on are those with high confidence annotation
#' tf with extended are those with relatively low confidence annotation
#' The point is to use the high annotation data if we have it, if not, we use
#' the ext type tf data
#' 
#' @return a tibble/dataframe with columns 'type' and 'ext'
#' type shows the tf name without any suffix -- ex. "E2f1"
#' ext shows whether that data is ext/regular
#' regular tf has its own name while extended data has 'ext' in that column
#' This feature is important in following functions, has_regular, tf_regular ...
#' to help retrieve the correct TF name to be used in other dataframes(tf activity in cell..)
#' 
#' @examples 
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' 
identify_tf <- function(TF_name_activity_tibble){
  TF_name_activity_tibble %>% 
    rename(name = value) %>%
    mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the weight
    mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
    mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
    mutate(type = str_replace(TF_type, "_ext", "")) %>%
    select(name, type, ext)
}

#' has_regular, has_ext
#' 
#' Boolean functions that determine the identity of TF data in the activity datasets
#'
#' @param TF character vector, containing one or more TF names
#' @param TF_and_ext The dataframe/tibble generated using identify_tf, that has three cols:
#' name, type, ext
#'
#' @return A boolean that checks if the datasets has the regular/ext data of that tf input
#'
#' @examples
#' 
#' has_regular("Arx", TF_and_ext) # False
#' has_ext("Arx", TF_and_ext) # True
#' 
#' 
has_regular <- function(TF, TF_and_ext){
  has_regular <- filter(TF_and_ext, type==TF & ext==TF)
  nrow(has_regular)!=0 #boolean
}
has_ext <- function(TF, TF_and_ext){
  is_ext <- filter(TF_and_ext, type==TF & ext=="ext")
  nrow(is_ext)!=0
}

#' tf_exist
#' 
#' This function also uses TF_and_ext, loaded in data_prep.R, can also generate using identify_tf
#' 
#' @param TF character vector, containing one or more TF names
#'
#' @return if the dataset doesn't have data for a tf input, it returns the tf input name; if it has them
#' all, this returns TRUE
#'
#' @examples 
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' TF <- "Arx"
#' tf_exist(TF)
#' 
#' 
tf_exist <- function(TF){
  for(tf in TF){
    if(has_regular(tf, TF_and_ext) || has_ext(tf, TF_and_ext)){}
    else{return (tf)}
  }
  return (TRUE)
}

# read the corresponding data by tf's identity

#' tf_regular, tf_ext
#'
#' @param TF character vector, containing one or more TF names
#' @param TF_and_ext The dataframe/tibble generated using identify_tf, that has three cols:
#' name, type, ext
#'
#' @return The best represented tf with suffix, still using the same logic: 
#' Use the high annotation data if we have it, if not, we use the ext type tf data
#' 
#'
#' @examples
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' TF <- c("Arx","Lef1")
#' tf_regular(TF, TF_and_ext)
#' 
tf_regular <- function(TF, TF_and_ext){
  filter(TF_and_ext, type==TF & ext==TF)[[1,1]]
}
tf_ext <- function(TF, TF_and_ext){
  filter(TF_and_ext, type==TF & ext=="ext")[[1,1]]
}

# --------------------------------Tab2 data--------------------------------
# NOTE: TF_and_ext is a dataframe (loaded already) that created in order to identify 
# whether the TF data is a regular TF (with high confidence annotation) 
# or ext type(with lower confidence)
# Dogma: if we have regular TF type, we use that data; if we only have ext data, use ext
# if we have no data related to this tf, we will either give an error message or do nothing

#' create Cell/Cluster activity data
#'
#' @param tf character vector, containing one or more TF names
#' @param method either by Cell --> use cell data, or by cluster --> use cluster data, 
#' this should be a string indicating the column name
#' @return a dataframe that has a column containing all the cell names and columns of the input tfs
#' the corresponding activity
#' data value (NES) 
#'
#' @examples
#' create_activity_data("Arx", "Cell")
#' create_activity_data("Arx", "Cluster")
create_activity_data <- function(tf, method){ 
  # use the feature of feather data to read certain col to optimize speed
  if(str_detect(method,"(?i)Cell")){
    cell_col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather",
                             "Cell")
  }
  else if(str_detect(method,"(?i)Cluster")){
    cell_col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cluster.feather",
                             "Cluster")
  }
  else{return("Wrong usage")}
  
  # add certain tf activity data to the Cell column
  activity_cell <- cell_col
  for(TF in tf){ # tf is input tf list, could contain many tfs
    tf_to_read <- TF
    if(has_regular(TF, TF_and_ext)){
      tf_to_read <- tf_regular(TF,TF_and_ext) # a helper to read the corresponding data
    }
    else if(has_ext(TF, TF_and_ext)){
      tf_to_read <- tf_ext(TF,TF_and_ext)
    }
    else{
      next # means we don't have that data, we jump over it and do thing
    }
    if(str_detect(method,"(?i)Cell")){
      col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cell.feather",
                          tf_to_read)
    }
    else{
      col <- read_feather("data/joint_cortex/joint_cortex.regulon_activity_per_cluster.feather",
                          tf_to_read)
    }
    
    activity_cell <- add_column(activity_cell,col)
  }
  activity_cell %>%
    select(method, everything()) # move cell col to start
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

# not used for now

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


#' make cell metadata of certain region, cortex/pon
#'
#' @param metadata_part a dataframe, the cell_metadata corresponding to a region: cortex/pon
#'
#' @return a dataframe/tibble of metadata with columns Age, Cell, Prefix and Cluster, used for plotting
#' timeseries
#'
#' @examples
#' cell_metadata_cortex <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")
#' cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)
#' 
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


#' Translate transcription factor name version
#' 
#' Since we have two datasets with different tf name format, (with ext and weights or not)
#' and those tf names are essential for retrieving data in various datasets,
#' this function use a clean tf and return a tf with ext and weight type
#'
#' @param tf a character vector of a tf without any suffix. EX: "Arx" 
#' @param tf_dataframe a one column dataframe that contains all the TF names with suffix
#' 
#' @return If the dataframe contains the tf input, it return a best represented tf name
#'  with ext and weight suffix. If not, it returns FALSE
#'
#' @examples
#' tf_df <- as_tibble(rownames(activity))
#' translate_tf("Arx",tf_df)  # Arx_extended (21g)
#' translate_tf("Brahl",tf_df) # FALSE
#' 
translate_tf <- function(tf, tf_dataframe){
  tf_info <- identify_tf(tf_dataframe)
  if(has_regular(tf, tf_info)){
    tf_regular(tf,tf_info) # a helper to read the corresponding data
  }
  else if(has_ext(tf, tf_info)){
    tf_ext(tf,tf_info)
  }
  else{
    return (FALSE) # means we don't have that data
  }
}

#' Plot timeseries
#'
#' @param TF a character vector that contains one or multiple TF names
#' @param cell_metadata cell_metadata_cortex, loaded in data_prep.R, a dataframe
#' with Age, Cell, Prefix, Cluster columns, this data is specific to forebrain cells/pon cells
#' @param activity loaded in data_prep.R
#'
#' @return a plot that displays the percentage level of that tf along with several the time points
#'
#' @examples
#' tf_df <- as_tibble(rownames(activity))
#' TF <- translate_tf("Lef1",tf_df)
#' binary_activity <- readRDS("data/joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")
#' cell_metadata_cortex <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")
#' cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)
#' plot_timeseries(TF,cell_metadata_cortex, binary_activity)
#' 
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





#TF_tbl1 <- as_tibble(TF_active)
#TF_tb2<- as_tibble(tf_list)

# missing_tf <- TF_tbl1 %>%
#   mutate(exist = if_else(value %in% tf_list, 1, 0)) %>%
#   filter(exist==0)
#
# mutual_tf <- TF_tbl1 %>%
#   mutate(exist = if_else(value %in% tf_list, 1, 0)) %>%
#   filter(exist==1)
  
  
  
  
  

