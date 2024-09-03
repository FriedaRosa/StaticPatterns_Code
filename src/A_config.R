## Config file for *A_PrepareData* for the Study on Static Patterns ##

# 1. Clean environment
rm(list=ls())
gc()








# 2. Libraries always needed
options(tidyverse.quiet = TRUE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE) # make planar
suppressPackageStartupMessages(library(rstatix))






# 3. Data paths
home_path <- "c:/Users/wolke/OneDrive - CZU v Praze/"

## Folders
source_atlas <- c(paste0(home_path, 
                         "Datasets/Processed/Atlases/Replicated/"))
source_predictors <- c(paste0(home_path, 
                              "Dokumenty/PhD_Projects/StaticPredictors/Data/"))
source_Git <- c(paste0(home_path, 
                       "Dokumenty/GitHub/StaticPredictors_Woelke_et_al_2024/"))

## Folder paths to atlas data
source_paths <- c(paste0(source_atlas, "Birds_Atlas_Czechia/"), 
                  paste0(source_atlas, "Birds_Atlas_New_York/"), 
                  paste0(source_atlas, "Birds_atlas_Japan/"), 
                  paste0(source_atlas, "Birds_atlas_EBBA/"))

## Folder path to output folder
out_path <- c(paste0(source_Git, "Results/"))

## Paths to data & grids
data_paths <- c(paste0(source_paths[1],"Birds_Atlas_Czechia_beast_data.rds"), 
                paste0(source_paths[2], "Birds_Atlas_New_York_beast_data.rds"), 
                paste0(source_paths[3], "Birds_atlas_Japan_beast_data.rds"),
                paste0(source_paths[4], "Birds_atlas_EBBA_beast_data_CHANGE.rds"))

grid_paths <- c(paste0(source_paths[1],"Birds_Atlas_Czechia_grid.gpkg"), 
                paste0(source_paths[2], "Birds_Atlas_New_York_grid.gpkg"), 
                paste0(source_paths[3], "Birds_atlas_Japan_grid.gpkg"),
                paste0(source_paths[4], "Birds_atlas_EBBA_grid.gpkg"))

# 4. Essential variables 
atlas_names <- c(
  "Birds_Atlas_Czechia", 
  "Birds_Atlas_New_York",
  "Birds_atlas_Japan", 
  "Birds_atlas_EBBA"
)

time_periods <- c(1,2)

# Define the desired order of factor levels
desired_levels <- factor(c("1", "2","4", "8", "16", "32", "64", "128"), 
                         ordered = T,  
                         levels = c("1", "2","4", "8", "16", "32", "64", "128")) 

crs_list <- c(CZ="epsg:5514", 
              NY="epsg:32118", 
              JP ="epsg:6684", 
              EBBA="epsg:3035") 









