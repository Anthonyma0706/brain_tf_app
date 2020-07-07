library(tidyr)
library(dplyr)
library(readr)
source(function.R)
# ———————————————————————————————————color palette————————————————————————————————————————
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


# ———————————————————————————————————Cortex data————————————————————————————————————————
forebrain_data <- read_tsv("data/joint_cortex/Forebrain_join.2D.tsv") %>% # for UMAP cluster
  mutate(Sample_cluster = str_replace(Sample_cluster," ","_"))
# clean some samples with space in between ...

TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))

# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene <- as_tibble(read_rds("data/joint_cortex/joint_cortex.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF <- unique(TF_target_gene[["TF"]])


TF_and_ext <- identify_tf(TF_active)

timeseries_input_meta_cortex <- create_metadata_timeseries(forebrain_data, "cortex")


# metadata specific for each cell, corresponding to the activity data
#cell_metadata_cortex_prep <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")

#cell_metadata_cortex_test <- create_cell_metadata_cortex(forebrain_data)

#cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex_prep)
# activity for cortex timeseries graph data
binary_activity <- readRDS("data/joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")
tf_df <- as_tibble(rownames(binary_activity)) #a dataframe that contains all the tf with 
# best representation of its identity in the binary_activity dataset


# ----------------------------------Pon data-------------------------------------------------------

pons_data <- read_tsv("data/joint_pons/Pons_join.2D.tsv") # for UMAP cluster

TF_active_pon <- as_tibble(read_rds("data/joint_pons/joint_pons.active_regulons.Rds"))

# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene_pon <- as_tibble(read_rds("data/joint_pons/joint_pons.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF_pon <- unique(TF_target_gene_pon[["TF"]])


TF_and_ext_pon <- identify_tf(TF_active_pon)

timeseries_input_meta_pons <- create_metadata_timeseries(pons_data,"pons") %>%
  filter(Cell != "___po_e12_TACGGGCGTCAAGCGA")


# metadata specific for each cell, corresponding to the activity data
#cell_metadata_pon_prep <- read_tsv("data/joint_pons/joint_pons.metadata.tsv")
# # 
#cell_metadata_pon <- create_cell_metadata_pon(cell_metadata_pon_prep) %>%
#    filter(Cell != "___po_e12_TACGGGCGTCAAGCGA")
# remove the extra line to make the number of cells the same as the binary activity pon data
# to correctly make the timeseires ribbon plot

# activity for cortex timeseries graph data
binary_activity_pon <- readRDS("data/joint_pons/joint_pons.binaryRegulonActivity_nonDupl.Rds")
tf_df_pon <- as_tibble(rownames(binary_activity_pon)) #a dataframe that contains all the tf with 
# best representation of its identity in the binary_activity dataset


# make two lists containing same name (will be assigned to a reactive list),
# then we can use the same name to code
data_cortex <- list(
  "cell_metadata"  = forebrain_data,
  "TF_and_ext" = TF_and_ext,
  "TF_target_gene_info" = TF_target_gene,
  "unique_active_TFs_bare" = unique_TF,
  "active_TFs" = TF_active,
  "binary_active_TFs" = tf_df,
  "timeseries_input_meta" = timeseries_input_meta_cortex,
  "binary_activity" = binary_activity

)

data_pons <- list(
  "cell_metadata" = pons_data,
  "TF_and_ext" = TF_and_ext_pon,
  "TF_target_gene_info" = TF_target_gene_pon,
  "unique_active_TFs_bare" = unique_TF_pon,
  "active_TFs" = TF_active_pon,
  "binary_active_TFs" = tf_df_pon,
  "timeseries_input_meta" = timeseries_input_meta_pons,
  "binary_activity" = binary_activity_pon
  
)


# ---------------------------cortex data-----------------------------
save(data_cortex, file = "data/joint_cortex/cortex_prep.Rda")

# -----------------------------pon data-----------------------------
save(data_pons, file = "data/joint_pons/pons_prep.Rda")

save(metadata,colour_palette_cluster,
     hm_anno,colour_palette, file = "data/joint_cortex/common_prep.Rda")





