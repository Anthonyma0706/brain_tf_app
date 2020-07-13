library(shiny)
library(tidyr)
library(tibble)
library(stringr)
library(dplyr)
library(readr)
library(feather)
library(plotly)
library(ggplot2)
library(pheatmap)
library(DT)
library(rcytoscapejs2)
library(glue)
load("data/joint_cortex/cortex_prep.Rda") # a list, data_cortex
load("data/joint_pons/pons_prep.Rda")     # a list, data_pons
load("data/joint_cortex/common_prep.Rda") # metadata and colour_palettes
source("functions.R")