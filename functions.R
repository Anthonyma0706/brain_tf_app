# --------------------------Tab1-cytoscape network visualization----------------------------------------------
# function to create network


#' Create rcytoscape network data
#' 
#' Takes a vector input that contains user selected TFs and output a list of nodeData and edgeData
#' which will be used for createCytoscapeJsNetwork
#' A good to visualize correlations among multiple TFs
#' 
#' @param tf one single tf name character
#' @param TF_target_gene TF_target_gene data, specific for cortex/pon
#' @param unique_TF unique_TF data, specific for cortex/pon
#'
#' @return a list of nodeData and edgeData that are required for generating a rcytoscapejs network object
#' 
#' @examples 
#' TF <- c("Arx","Lef1")
#' # Note that TF_target_gene and unique_TF will be saved in data_cortex list, by data_prep.R
#' nodeData <- create_network(TF, TF_target_gene_pon, unique_TF)$nodes
#' edgeData <- create_network(TF, TF_target_gene_pon, unique_TF)$edges
#' network <- createCytoscapeJsNetwork(nodeData, edgeData)
#' rcytoscapejs2(network$nodes, network$edges)
#' 
create_network <- function(tf, TF_target_gene, unique_TF, pathway_genes = c(),
                           shrink_gray = FALSE){ 
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
    # you can customize the color using the case_when structure easily,
    # check the tfs in id column that exist in your vector, then you can control its size,
    # shape and color easily
    mutate(color = case_when(id %in% tf ~ "#9d4097", # orange
                             # orange nodes are tfs that are active in this region
                             id %in% pathway_genes ~ "green",
                             id %in% unique_TF ~ "#D6604D", 
                             id %in% mutual_target ~ "#4fafc6",
                             TRUE ~ "lightgrey")) %>%
    mutate(height = case_when(id %in% tf ~ "100",
                           TRUE ~ "70")) %>%
    mutate(width = case_when(id %in% tf ~ "100",
                            TRUE ~ "70"))
  
  if(shrink_gray){
    nodeData <- nodeData %>%
      mutate(height = case_when(color %in% "lightgrey" ~ "40",
                                TRUE ~ "70")) %>%
      mutate(width = case_when(color %in% "lightgrey" ~ "40",
                               TRUE ~ "70"))
    
  }
  
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
#' @param TF_and_ext TF_and_ext data, specific for cortex/pons
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
tf_exist <- function(TF, TF_and_ext){
  for(tf in TF){
    if(has_regular(tf, TF_and_ext) || has_ext(tf, TF_and_ext)){}
    else{return (tf)} # or regurn FALSE
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
#' Make activity data used in tab2 by either Cell or Cluster, the method would be provided by
#' user's input in Shiny app
#' This function uses feather file that will be read by a certain col to maximize speed,
#' so we switch the paths of the feather file for different brain region
#' 
#' @param tf character vector, containing one or more TF names
#' @param method either by Cell --> use cell data, or by cluster --> use cluster data, 
#' this should be a string indicating the column name
#' @param TF_and_ext TF_and_ext data, specific for cortex/pons
#' @return a dataframe that has a column containing all the cell names and columns of the input tfs
#' the corresponding activity
#' data value (NES) 
#'
#' @examples
#' create_activity_data("Arx", "Cell", "cortex", TF_and_ext)
#' create_activity_data("Pax6", "Cluster", "pons", TF_and_ext_pon)


create_activity_data <- function(tf, method, region, TF_and_ext){ 
  # use the feature of feather data to read certain col to optimize speed
  #if(tf_exist(tf, TF_and_ext) != TRUE){return("TF does not exist")}
  if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons")
  
  method2 <- str_to_lower(method)

  # set up the path of the feather file to read
  path <- glue('data/joint_{region}/joint_{region}.regulon_activity_per_{method2}.feather')
  
  # case-insensitive checking
  if(str_detect(method,"(?i)Cell")){cell_col <- read_feather(path, "Cell")}
  else if(str_detect(method,"(?i)Cluster")){cell_col <- read_feather(path,"Cluster")}
  else{return("Wrong usage, method should be Cell/Cluster")}
  
  # add certain tf activity data to the Cell column
  activity <- cell_col
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
    
    if(str_detect(method,"(?i)Cell")){col <- read_feather(path,tf_to_read)}
    else{col <- read_feather(path,tf_to_read)}
    
    activity <- add_column(activity,col)
  }
  activity %>%
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


#' Plot heatmap by cluster/cells
#'
#' @param tf 
#' @param method 
#' @param region 
#' @param TF_and_ext 
#' @param brain_data either forebrain_data or pon_data, eventually will be saved by data_prep.R
#' and loaded at the beginning of app.R as an element in a list
#'
#' @return
#' @export
#'
#' @examples
#' plot_heatmap(c("Arx","Lef1"), "Cluster","cortex", TF_and_ext,forebrain_data)
#' plot_heatmap(c("Arx","Lef1"), "Cell","cortex", TF_and_ext,forebrain_data)
#' plot_heatmap(c("Pax6","Lef1"), "Cluster","pons", TF_and_ext_pon, pon_data)
#' plot_heatmap(c("Pax6","Lef1"), "Cell","pons", TF_and_ext_pon,pon_data)
#' 
plot_heatmap <- function(tf,method, region, TF_and_ext,brain_data, cell_plot_num = 300){
  # sanity checking
  if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons")
  if(!method %in% c("Cell","Cluster")) return("Wrong usage, method should be Cell/Cluster")
  
  if(method == "Cell"){
    # 1. create the activity data for plotting 
    act_cell <- create_activity_data(tf, "Cell",region, TF_and_ext) %>%
      mutate(Cluster = gsub("_"," ",brain_data[["Sample_cluster"]])) %>%
      filter(!grepl("BLACKLIST", Cluster)) %>% # filter out bad samples
      sample_n(cell_plot_num) %>%  # randomly sample it
      tibble::column_to_rownames(var = "Cell") # make that column name as row name ...
    
    anno_row_cell <- select(act_cell, Cluster)
    # change the anno_row, since we change the color palettes
    new_anno_row_cell <- anno_row_cell %>%
      mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
    rownames(new_anno_row_cell) <- rownames(anno_row_cell) # re-assign the rownames
    
    act <- select(act_cell, -Cluster) # must remove Cluster data before plotting
    
    # customized for plotting by cell
    anno_col <- new_anno_row_cell # assign to the same variable for plotting
    cell_width_plot <- 2
    show_colname_plot <- FALSE
  }
  else if(method == "Cluster"){
    act <- create_activity_data(tf, "Cluster",region, TF_and_ext) %>%
      #sample_n(cluster_plot_num) %>% # randomly sample it
      tibble::column_to_rownames(var = "Cluster") # make that column name as row name ...
    
    # change the anno_row, since we change the color palettes
    new_anno_row <- hm_anno$anno_row %>%
      mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
    rownames(new_anno_row) <- rownames(hm_anno$anno_row) # re-assign the rownames
    # note that the rownames correspond to the col names of the matrix t(act_cluster)
    # customized for plotting by cluster
    anno_col <- new_anno_row # this is loaded by data_prep.R
    cell_width_plot <- 10
    show_colname_plot <- TRUE
  }
  pheatmap::pheatmap(t(act),
                     show_colnames = show_colname_plot,
                     scale = "none",
                     border_color = NA,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = glue('Plot by {method}s'),
                     annotation_col = anno_col,
                     # change the default color annotation
                     annotation_colors = hm_anno_new$side_colors, # loaded by data_prep.R
                     annotation_legend = FALSE,
                     cellwidth = cell_width_plot,
                     cellheight = 10)
} 


#' Make UMAP clustering scatterplot
#'
#' @param tf_number Either 1 or 2. In the tf input vector we get from user in Shiny app, there could be
#' multiple tfs, but we only support plotting two tfs since the scatterplot is big
#' @param overall_brain_data metadata (forebrain_data or pon_data), saved in data_cortex
#' and data_pons
#' @param cell_activity_data made by make_cell_metadata given the tf input
#' @param sample_number we eliminate half of the cell samples to relieve the burden of 
#' the RAM to speed up plotting, since we have over 37000 cells(samples), we randomly sample 
#' 13000 to optimize speed, but one can also specify this value to see fewer or more sample points
#'
#' @return a UMAP scatter plot that shows in which cluster(region) the tf expresses the most
#'
#' @examples
#' tf <- c("Arx","Lef1")
#' activity_test_tf1 <- create_activity_data(tf, "Cell","cortex", data_cortex$TF_and_ext)
#' plot_UMAP(tf_number = 1,data_cortex$overall, activity_test_tf1)
#' 
plot_UMAP <- function(tf_number, cell_metadata, cell_activity_data, sample_number = 13000,
                      sample_reduce = TRUE){
  if(tf_number == 1) tf_plot <- 2 # number of col, the first col is Cell, so start from 2
  else if(tf_number == 2) tf_plot <- 3
  else{return(
    "Wrong usage, now we only support plotting two tfs since the scatterplot is big"
  )}

  if(! sample_reduce) sample_number <- 27000
  
  activity_tf <- cell_activity_data[,tf_plot][[1]]
  cell_meta_with_activity <- mutate(cell_metadata, activity_tf = activity_tf) %>%
    sample_n(sample_number)
  ggplot(data = cell_meta_with_activity, mapping = aes(x=UMAP1,y=UMAP2))+
    geom_point(aes(color = activity_tf))+
    scale_color_gradientn(colors = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100)))+
    theme_bw()
  
}



#' make cell metadata of certain region, cortex/pon
#'
#' @param cell_metadata a dataframe, forebrain_data or pons_data, saved in data_cortex / data_pons
#' @param part  a string, either "cortex" or "pons"
#' @return a dataframe/tibble of metadata with columns Age, Cell, Prefix and Cluster, used for plotting
#' timeseries
#'
#' @examples
#' cell_metadata_cortex <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")
#' cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)
#' 
create_metadata_timeseries <- function(cell_metadata, part){
  if(part == "cortex") level <- c("Forebrain E12.5",
                                  "Forebrain E15.5",
                                  "Forebrain P0",
                                  "Forebrain P3",
                                  "Forebrain P6")
  else if (part == "pons") level <- c("Hindbrain E12.5",
                                      "Pons E15.5",
                                      "Pons P0",
                                      "Pons P3",
                                      "Pons P6")
  else{(return("Wrong usage, input either cortex or pons"))}
  
  cell_metadata %>% 
    select(Age = Sample, Cell, Cluster = Sample_cluster) %>% 
    # In this case, we remove the "prefix" of the Cluster column, so that we are
    # simply left with the abbreviation representing the cell type, so that 
    # we can link the cells of the same cell type across ages
    separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
    mutate(Age = factor(Age, levels = level)) %>% 
    arrange(Cell)
  
}

# create_cell_metadata_pon <- function(metadata_part){
#   metadata_part %>% 
#     select(Age = orig.ident, Cell, Cluster = ID_20190715_with_blacklist_and_refined) %>% 
#     # In this case, we remove the "prefix" of the Cluster column, so that we are
#     # simply left with the abbreviation representing the cell type, so that 
#     # we can link the cells of the same cell type across ages
#     separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
#     mutate(Age = factor(Age, levels = c("Hindbrain E12.5",
#                                         "Pons E15.5",
#                                         "Pons P0",
#                                         "Pons P3",
#                                         "Pons P6"))) %>% 
#     arrange(Cell)
#   
# }


#' Translate transcription factor name version
#' 
#' Since we have two datasets with different tf name format, (with ext and weights or not)
#' and those tf names are essential for retrieving data in various datasets,
#' this function use a clean tf and return a tf with ext and weight type
#'
#' @param tf a character vector of a tf without any suffix. EX: "Arx" 
#' @param tf_dataframe a one column dataframe that contains all the TF names with suffix 
#' , get from the rownames of the cell binary activity data for timeseries tab3
#' 
#' @return If the dataframe contains the tf input, it return a best represented tf name
#'  with ext and weight suffix. If not, it returns FALSE
#'
#' @examples
#' tf_df <- data_cortex$ # as_tibble(rownames(activity))
#' translate_tf("Arx",tf_df)  # Arx_extended (21g)
#' translate_tf("Brahl",tf_df) # NULL
#' translate_tf(c("Lef1","Arx"),tf_df) # "Lef1 (22g)"         "Arx_extended (21g)"
translate_tf <- function(tf_list, tf_dataframe){
  tf_info <- identify_tf(tf_dataframe)
  l <- c()
  for(TF in tf_list){
    if(has_regular(TF, tf_info)){
      l <- c(l, tf_regular(TF,tf_info)) # a helper to read the corresponding data
    }
    else if(has_ext(TF, tf_info)){
      l <- c(l, tf_ext(TF,tf_info))
    }
    else{
      next
    }
  }
  if(is.null(l)) return (FALSE) # means we don't have that data at all
  else{return (l)} # return the list
}



#' Plot timeseries
#' @author Selin Jessa and Anthony Ma, most credit to Selin and Anthony puts codes into the function
#' @param TF a character vector that contains one or multiple TF names, that may need to be 
#' transformed by translate_tf function to change its string form
#' @param cell_metadata cell_metadata_cortex, loaded in data_prep.R, a dataframe
#' with Age, Cell, Prefix, Cluster columns, this data is specific to forebrain cells/pon cells
#' @param activity binary_activity data, loaded in data_prep.R
#'
#' @return a plot that displays the percentage level of that tf along with several the time points
#'
#' @examples
#' tf_df <- as_tibble(rownames(activity))
#' TF <- translate_tf("Lef1",tf_df)
#' binary_activity <- data_cortex$binary_activity
#' cell_metadata_cortex <- create_metadata_timeseries(data_cortex$cell_metadata, "cortex")
#' plot_timeseries(TF,cell_metadata_cortex, binary_activity)
#' 
plot_timeseries <- function(TF,cell_metadata, activity, make_plotly = FALSE, show_legend = TRUE){
  
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
  
  plot <- df %>%
    ggplot(aes(x = xpos, y = frac, fill = cluster)) +
    geom_area(stat = "identity", show.legend = show_legend) +
    scale_fill_manual(values = colour_palette, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(unique(df$stage)),
                       labels = c("E12.5", "E15.5", "P0", "P3", "P6"),
                       limits = c(1, length(unique(df$stage)))) +
    labs(x = "age", y = "Proportion", title = TF) +
    guides(fill = guide_legend(ncol = 5)) +
    theme(legend.position = "bottom")
  
  if(make_plotly) {
    return (ggplotly(plot))
  }
  else{return(plot)}
}


