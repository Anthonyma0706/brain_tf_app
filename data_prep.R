library(tidyr)
library(dplyr)
library(readr)
forebrain_data <- read_tsv("data/joint_cortex/Forebrain_join.2D.tsv") # for UMAP cluster
TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))

# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene <- as_tibble(read_rds("data/joint_cortex/joint_cortex.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF <- unique(TF_target_gene[["TF"]])


TF_and_ext <- identify_tf(TF_active)
# TF_and_ext <- TF_active %>% 
#   rename(name = value) %>%
#   mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the gram
#   mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
#   mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
#   mutate(type = str_replace(TF_type, "_ext", "")) %>%
#   select(name, type, ext)

# make color palette
metadata <- read_tsv("data/joint_cortex/metadata_20190716.tsv")

# color palette for heatmap
colour_palette_cluster <- metadata %>% 
  # use gsub to change all contents in Cluster (cluster name format)
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  # Get two columns
  select(Cluster, Colour) %>% 
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe() # VECTOR , not data frame 

# A helper function to prepare a dataframe to annotate the heatmap with colours
hm_anno <- makePheatmapAnno(colour_palette_cluster, "Cluster")

# color palette for timeseries plot, tab3
colour_palette <- metadata %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  separate(Cluster, into = c("Prefix", "Cluster"), sep = " ") %>% 
  # Get two columns
  select(Cluster, Colour) %>% 
  distinct(Cluster, .keep_all = TRUE) %>% 
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe()



# metadata specific for each cell, corresponding to the activity data
cell_metadata_cortex_prep <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")

cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex_prep)
# activity for cortex timeseries graph data
binary_activity <- readRDS("data/joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")
tf_df <- as_tibble(rownames(binary_activity)) #a dataframe that contains all the tf with 
# best representation of its identity in the binary_activity dataset

save(forebrain_data,TF_and_ext,TF_active,metadata,tf_df,TF_target_gene, unique_TF,
     colour_palette_cluster,hm_anno,colour_palette,
     cell_metadata_cortex,binary_activity, file = "data/joint_cortex/cortex_prep.Rda")







