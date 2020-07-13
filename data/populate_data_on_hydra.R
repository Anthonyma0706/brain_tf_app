
# Script to populate the data directory for this app.
#   - Reads the data ontology in data.json
#   - Creates necessary subdirectories
#   - Identifies which files need to be copied over from elsewhere on the server
#   - Prints a command to copy the files to the right place
#
# Usage:
# 1. Run the script from within the data directory with:
# $ Rscript populate_data_on_hydra.R
#
# 2. Copy and execut the cp commands
#
# 3. Run the data_prep.R within the data directory with:
# $ Rscript data_prep.R

library(dplyr)
library(tidyr)
library(magrittr)
library(rjson)
library(purrr)
library(glue)

path_to_projects <- "/mnt/KLEINMAN_JBOD1/SCRATCH/projects/"

# Read in data.json, returning a list with one element per folder
app_data <- rjson::fromJSON(file = "data.json")

# Get the names of the directories
directories <- names(app_data)

# Create directories
walk(directories, dir.create, showWarnings = FALSE)

# For each list,
#     look for the files which are not produced by a script, and copy those over
process_file <- function(file, dir_name) {
  
  # If not, copy over the file
  if (is.null(file$script)) {
    
    # Get the full path to the file
    src  <- file.path(path_to_projects, file$path)
    dest <- file.path(dir_name, file$file)
    
    # Copy the file
    cmd <- glue("cp {src} {dest}")
    
    # Echo the command
    message(cmd)
    
  } else {
    
    NULL
    
  }
  
}

# Loop over the directories
iwalk(app_data, function(dir_data, dir_name) {
  
  # Loop over files, and get the copy commands, if it doesn't have a script property
  walk(dir_data, ~ process_file(.x, dir_name = dir_name))
  
})


