library(tidyr)
library(dplyr)
library(readr)
TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))

TF_and_ext <- TF_active %>% 
  rename(name = value) %>%
  mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the gram
  mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
  mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
  mutate(type = str_replace(TF_type, "_ext", "")) %>%
  select(name, type, ext)

# make color palette
metadata <- read_tsv("data/joint_cortex/metadata_20190716.tsv")

colour_palette <- metadata %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  separate(Cluster, into = c("Prefix", "Cluster"), sep = " ") %>% 
  # Get two columns
  select(Cluster, Colour) %>% 
  distinct(Cluster, .keep_all = TRUE) %>% 
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe()
#head(colour_palette)
# metadata specific for each cell, corresponding to the activity data
cell_metadata_cortex <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")
cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)
# activity for cortex timeseries graph data
activity <- readRDS("data/joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")
tf_df <- as_tibble(rownames(activity))
save(TF_and_ext,TF_active,metadata,tf_df,
     cell_metadata_cortex,activity, file = "data/joint_cortex/cortex_prep.Rda")







