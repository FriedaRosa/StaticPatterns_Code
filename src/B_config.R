## Config file for *B_MachineLearning* for the Study on Static Patterns ##

# 1. Clean environment
rm(list=ls())
gc()

source("../src/functions.R")

# 3. Data paths
home_path <- "c:/Users/wolke/OneDrive - CZU v Praze/"

## Folders
source_Git <- c(paste0(home_path, 
                       "Dokumenty/GitHub/StaticPredictors_Woelke_et_al_2024/"))

## Folder path to output folder
out_path <- c(paste0(source_Git, "Results/"))


atlas_names <- c(
  "Birds_Atlas_Czechia", 
  "Birds_Atlas_New_York",
  "Birds_atlas_Japan", 
  "Birds_atlas_EBBA"
)



## MachineLearning Packages
pckgs <- c("dplyr", "ggplot2", "reshape2", 
"ggcorrplot", 
"caret",  "recipes",   "caretEnsemble", 
"randomForest", "ranger", "gbm", "xgboost", 
"vegan", "pdp", 
"gridExtra", "kableExtra")

suppressMessages(install_and_load(pckgs))

# Hypotheses
H1_vars <- c(
  "sd_PC1", "sd_PC2", # Climatic Niche Breadth
  "GlobRangeSize_m2", "IUCN", "Mass", "Habitat", "Habitat.Density",
  "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle",
  "FP", # Phylogenetic Distinctness
  "Hand.Wing.Index") # Measure of dispersal ability
H2_vars <- c(
  "AOO", "rel_occ_Ncells", "mean_prob_cooccur", "D_AOO_a", 
  "moran", "x_intercept", "sp_centr_lon", "sp_centr_lat",
  "lengthMinRect", "widthMinRect", "elonMinRect", "bearingMinRect",
  "circ", "bearing", "Southernness", "Westernness",
  "rel_maxDist", "rel_ewDist", "rel_nsDist", "rel_elonRatio",
  "rel_relCirc", "rel_circNorm", "rel_lin", "Dist_centroid_to_COG",
  "maxDist_toBorder_border", "maxDist_toBorder_centr",
  "minDist_toBorder_centr")
H3_vars <- c("GammaSR", "AlphaSR_sp", "BetaSR_sp")
H4_vars <- c(
  "dataset", "mean_area", "Total_area_samp", "Total_Ncells_samp",
  "mean_cell_length", "atlas_lengthMinRect", "atlas_widthMinRect",
  "atlas_elonMinRect", "atlas_circ", "atlas_bearingMinRect",
  "atlas_bearing", "AtlasCOG_long", "AtlasCOG_lat")


# Labels
atlas_label <- c("Czech Republic", "New York State", "Japan", "Europe")
names(atlas_label) <- atlas_names
