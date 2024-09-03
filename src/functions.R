## 1. Functions for (A) Data Preparation.

 #                   _____            _               _____                        
 #     /\      _    |  __ \          | |             |  __ \                       
 #    /  \    (_)   | |  | |   __ _  | |_    __ _    | |__) |  _ __    ___   _ __  
 #   / /\ \         | |  | |  / _` | | __|  / _` |   |  ___/  | '__|  / _ \ | '_ \ 
 #  / ____ \   _    | |__| | | (_| | | |_  | (_| |   | |      | |    |  __/ | |_) |
 # /_/    \_\ (_)   |_____/   \__,_|  \__|  \__,_|   |_|      |_|     \___| | .__/ 
 #                                                                          | |    
 #                                                                          |_|    

## ======================================= ##
#                                           #
#         Function 1: Process Grids         #
#                                           #
## ======================================= ##

process_grids <- function(grid_paths, atlas_names, desired_levels) {
  
  # Step 1: Create the layers list
  layers_list <- list()
  
  for (i in seq_along(grid_paths)) {
    layers <- st_layers(grid_paths[i])$name
    layers_list[[i]] <- layers
  }
  
  names(layers_list) <- atlas_names

  # Step 2: Read and process grids to a list
  grid_list <- list()

  for (a in seq_along(grid_paths)) {
    grids <- sapply(layers_list[[a]], function(i) {
      st_read(grid_paths[[a]], paste(i), quiet = TRUE) %>% 
        reorder_levels(cell_grouping, order = desired_levels) %>% 
        mutate(atlas_cell_repeated = case_when(
          area1s != 0 & area2s != 0 ~ 2, # twice
          (area1s != 0 & area2s == 0) | (area1s == 0 & area2s != 0) ~ 1, # once
          area1s == 0 & area2s == 0 ~ 0, # never
          TRUE ~ NA)
          )
      }, simplify = FALSE)
    grid_list[[a]] <- grids
  }
  
  # Step 3: Combine the processed grids
  all_repeated_df <- do.call(rbind, lapply(seq_along(grid_list), function(j) {
    
    grids <- grid_list[[j]]
    do.call(rbind, lapply(grids, function(grid) {
      grid %>% 
        st_drop_geometry() %>% 
        distinct(cell_grouping, cell_label, atlas_cell_repeated, cell_lat, cell_long, area) %>% 
        mutate(dataset = atlas_names[j])
      }))
    }))

  rownames(all_repeated_df) <- NULL

  # Return the result
  return(list(all_repeated_df = all_repeated_df, grid_list = grid_list))
}








## ======================================= ##
#                                           #
#         Function 2: Process Data          #
#                                           #
## ======================================= ##


process_species_data <- function(data_paths, all_repeated_df, atlas_names, desired_levels) {
# Initialize the list
presence_data <- list()

# Loop through data paths
for (i in seq_along(data_paths)) {
  pres_dat <- readRDS(data_paths[i])
  sy <- sort(unique(pres_dat$start_year))
  
  # Read and process data
  pres_dat2 <- pres_dat %>%     
    ungroup() %>%
    mutate(tp = case_when(start_year == sy[1] ~ 1, 
                          start_year == sy[2] ~ 2)) %>% 
    filter(tp %in% c(1, 2)) %>%
    reorder_levels(cell_grouping, order = desired_levels) %>% 
    select(dataset, tp, cell_grouping, 
           cell_label, cell_lat, cell_long, area, 
           verbatim_name)
  
  # Combine with all_repeated_df
  pres_dat3 <- all_repeated_df %>% 
    filter(dataset == atlas_names[i]) %>% 
    group_by(cell_grouping, dataset) %>% 
    left_join(pres_dat2)
  
  # Store in list
  presence_data[[i]] <- pres_dat3
}

# Merge the list of atlases together
presence_data2 <- plyr::rbind.fill(presence_data, fill = TRUE) %>% 
  distinct(dataset, tp, verbatim_name, cell_grouping, cell_label, .keep_all = TRUE)

# Convert columns to factors
presence_data2$tp <- as.factor(presence_data2$tp)
presence_data2$dataset <- as.factor(presence_data2$dataset)
presence_data2$cell_grouping <- as.factor(presence_data2$cell_grouping)
presence_data2$atlas_cell_repeated <- as.factor(presence_data2$atlas_cell_repeated)

# Return the result
return(presence_data2)
}










## ======================================= ##
#                                           #
#   Function 3: Calculate grid information  #
#                                           #
## ======================================= ##

calculate_grid_info <- function(grid) {
  grid_info <- cbind(
    grid %>%
      st_drop_geometry() %>%
      summarise(
        Total_area = sum(area),
        Total_Ncells = n_distinct(cell_label)
      ),
    grid %>%
      filter(atlas_cell_repeated == 2) %>%
      st_drop_geometry() %>%
      summarise(
        Total_area_samp = sum(area),
        Total_Ncells_samp = n_distinct(cell_label)
      )
  )
  
  return(grid_info)
}










## ======================================= ##
#                                           #
#   Function 4: Jaccard Formula             #
#                                           #
## ======================================= ##

# Define the Jaccard function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}




## ======================================= ##
#                                           #
#  Function 5: Calculate Occupancy          #
#                                           #
## ======================================= ##

calculate_occupancy <- function(pres_dat_final, atlas_names) {
  
  occ_data_final <- pres_dat_final %>%
    ungroup() %>% 
    distinct(dataset, tp, cell_grouping, verbatim_name, cell_label, .keep_all = TRUE) %>%
    group_by(dataset, tp, cell_grouping, verbatim_name) %>%
    
    # Calculate AOO:
    dplyr::mutate(
      mean_area = mean(area),
      AOO = sum(area),
      occ_Ncells = n_distinct(cell_label)) %>%
    
    # Calculate relative Occupancy:
    dplyr::mutate(
      rel_AOO = AOO / Total_area_samp,
      rel_occ_Ncells = occ_Ncells / Total_Ncells_samp) %>%  # Prevalence
    
    # Remove duplicated rows:
    distinct() %>%
    
    # Add scale
    ungroup() %>%
    group_by(dataset) %>%
    dplyr::mutate(
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "1", 1/64, NA),
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "2", 1/32, scale), 
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "4", 1/16, scale), 
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "8", 1/8, scale),
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "16", 1/4, scale),
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "32", 1/2, scale),
      scale = ifelse(dataset %in% c(atlas_names[1]) & cell_grouping == "64", 1, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "1", 1/128, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "2", 1/16, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "4", 1/32, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "8", 1/16, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "16", 1/8, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "32", 1/4, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "64", 1/2, scale),
      scale = ifelse(dataset %in% c(atlas_names[2:4]) & cell_grouping == "128", 1, scale))
  

  return(occ_data_final)
}






## ======================================= ##
#                                           #
#   Function 6: OAR & Fractal Dimension     #
#                                           #
## ======================================= ##

process_and_calculate_OAR <- function(species_data, atlas_names, time_periods) {
  list_sp <- list()
  list_tp <- list()
  list_a <- list()
  
  for (a in seq_along(atlas_names)) {
    temp_df <- species_data %>% filter(dataset == atlas_names[a])
    
    for (t in seq_along(time_periods)) {
      temp_df_t <- temp_df %>% filter(tp == time_periods[t])
      species_names <- unique(temp_df_t$verbatim_name)
      
      for (s in seq_along(species_names)) {
        temp_df_s <- temp_df_t %>% filter(verbatim_name == species_names[s])
        
        # Exclude saturated scales
        temp_df_red <- temp_df_s %>% filter(rel_AOO < 1)
        
        if (nrow(temp_df_red) < 2 & nrow(temp_df_red) > 0) {
          out_df <- data.frame(verbatim_name = species_names[s],
                               dataset = atlas_names[a],
                               tp = time_periods[t],
                               exclude_sp_OAR = 1,
                               available_scales = nrow(temp_df_red),
                               mean_relAOO = mean(temp_df_red$rel_AOO, na.rm=T))
        } else if (nrow(temp_df_red) >= 2) {
          out_df <- data.frame(verbatim_name = species_names[s],
                               dataset = atlas_names[a],
                               tp = time_periods[t],
                               exclude_sp_OAR = 0,
                               available_scales = nrow(temp_df_red),
                               mean_relAOO = mean(temp_df_red$rel_AOO, na.rm=T))
        } else {
          out_df <- data.frame(verbatim_name = species_names[s],
                               dataset = atlas_names[a],
                               tp = time_periods[t],
                               exclude_sp_OAR = 1,
                               available_scales = nrow(temp_df_red),
                               mean_relAOO = 1)
        }
        list_sp[[s]] <- out_df
      }
      sp_df <- plyr::rbind.fill(list_sp)
      list_tp[[t]] <- sp_df
    }
    tp_df <- plyr::rbind.fill(list_tp)
    list_a[[a]] <- tp_df
  }
  
  atlas_df <- plyr::rbind.fill(list_a)
  atlas_df$tp <- as.factor(atlas_df$tp)
  
  sp_data_new <- full_join(species_data, atlas_df) %>% distinct(dataset, tp, cell_grouping, verbatim_name, .keep_all=T)
  
  dd <- sp_data_new %>% 
    filter(exclude_sp_OAR == 0) %>% 
    filter(rel_occ_Ncells < 1) %>% 
    unique() %>% 
    filter_at(vars(c(cell_grouping, AOO, mean_area)), any_vars(!is.na(.)))
  
  OAR_list_sp <- list()
  OAR_list_sp_tp <- list()
  OAR_list_sp_tp_atlas <- list()
  
  for (n_atlas in seq_along(atlas_names)) {
    atlas_name <- atlas_names[n_atlas]
    atlas <- dd %>% filter(dataset == atlas_name)
    
    for (time in seq_along(time_periods)) {
      period_nr <- time_periods[time]
      atlas_1time <- atlas %>% filter(tp == period_nr)
      sp <- unique(atlas_1time$verbatim_name)
      
      for (spec in seq_along(sp)) {
        species <- sp[spec]
        model_df <- atlas_1time %>% 
          filter(verbatim_name == species) %>% 
          distinct()
        
        OAR <- lm(log(AOO) ~ log(mean_area), data = model_df)
        
        OAR_df <- data.frame(
          verbatim_name = species,
          dataset = atlas_name,
          tp = period_nr,
          m_AOO_a = OAR$coefficients[2],
          b_AOO_a = OAR$coefficients[1])
        
        OAR_df$D_AOO_a <- -2 * OAR_df$m_AOO_a + 2
        OAR_df$dataset <- as.factor(OAR_df$dataset)
        OAR_list_sp[[spec]] <- OAR_df
      }
      OAR_df_sp <- plyr::rbind.fill(OAR_list_sp, fill = T)
      OAR_list_sp_tp[[time]] <- OAR_df_sp
    }
    OAR_df_sp_tp <- plyr::rbind.fill(OAR_list_sp_tp, fill = T)
    OAR_list_sp_tp_atlas[[n_atlas]] <- OAR_df_sp_tp
  }
  
  OAR_final <- plyr::rbind.fill(OAR_list_sp_tp_atlas, fill = T) %>% distinct()
  species_data_new <- merge(sp_data_new, OAR_final, by = intersect(names(sp_data_new), names(OAR_final)), all = T) %>% distinct()
  
  return(species_data_new)
}



## ======================================= ##
#                                           #
#     Function 7: Spatial Geometries        #
#     Function Credits: XXX ANONYMIZED      #
#                                           #
## ======================================= ##


# Function to get the max elongation per polygon, north-south length or east-west length, circularity, etc...
poly_attr <- function(x, type = NULL) {
  # Convert "sf" to SpatVector
  if (inherits(x, "sf") == T) {
    x <- vect(x)
  }
  # Reproject to WGS84. Currently, "terra" calculate distance in meters when given latitude and longitude points.
  x <- project(x, "epsg:4326")
  # Get East-West distance (longitudinal extension)
  if (type == "ewDist") {
    vert <- crds(as.points(ext(x)), df = T)
    dimNames <- list(c("point"),c("x","y"))
    sw <- matrix(c(min(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    se <- matrix(c(max(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    res <- distance(sw, se, lonlat = T)[[1]]
  }
  # Get South-North distance (latitudinal extension)
  if (type == "nsDist") {
    vert <- crds(as.points(ext(x)), df = T)
    dimNames <- list(c("point"),c("x","y"))
    sw <- matrix(c(min(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    nw <- matrix(c(min(vert[["x"]]), max(vert[["y"]])), ncol = 2, dimnames = dimNames)
    res <- distance(sw, nw, lonlat = T)[[1]]
  }
  # Get max distance between opposite vertices
  if (type == "maxDist") {
    vert <- crds(as.points(convHull(x)))
    res <- distance(vert, lonlat = T) %>% max()
  }
  # Get elongation ratio along longest axis
  if (type == "elonRatio") {
    convexHull <- convHull(x)
    vert <- crds(as.points(convexHull), df = T)
    dist <- as.data.frame.table(as.matrix(distance(vert, lonlat = T)), responseName = "distance") %>%
      slice_max(., distance)
    axisPoints <- vert[c(dist[[1]][[1]],dist[[2]][[1]]),]
    axisPoints <- arrange(axisPoints, desc(y))
    rotation <- -1*bearing(axisPoints[2,],axisPoints[1,])
    rotHull <- spin(convexHull, rotation)
    ext <- ext(convexHull)
    df <- as.vector(distance(crds(as.points(ext)), lonlat = T))
    df <- sort(df)
    length <- mean(df[[3]], df[[4]])
    width <- mean(df[[1]], df[[2]])
    res <- 1 - (width / length)
    # res <- list()
    # res[["vert"]] <- vert
    # res[["dist"]] <- dist
    # res[["axispoints"]] <- axisPoints
    # res[["bearing"]] <- rotation
    # res[["rotHull"]] <- rotHull
  }
  # Circularity
  if (type == "circ") {
    perimeter <- perim(x)
    area <- expanse(x)
    res <- (perimeter^2) / area
  }
  # Normalized circularity
  if (type == "circNorm") {
    perimeter <- perim(x)
    area <- expanse(x)
    res <- (perimeter^2) / (4 * pi * area)
  }
  # Major length of the minimum rectangle
  if (type == "lengthMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle), df = T), lonlat = T))
    df <- sort(df)
    res <- mean(df[[3]], df[[4]])
  }
  # Width of the minimum rectangle
  if (type == "widthMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle), df = T), lonlat = T))
    df <- sort(df)
    res <- mean(df[[1]], df[[2]])
  }
  # Elongation ratio of minimal encasing rectangle (from here: Dražić, S., Ralević, N., & Žunić, J. (2010). Shape elongation from optimal encasing rectangles. Computers & Mathematics with Applications, 60(7), 2035–2042. https://doi.org/10.1016/j.camwa.2010.07.043)
  if (type == "elonMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle)), lonlat = T))
    df <- sort(df)
    length <- mean(df[[3]], df[[4]])
    width <- mean(df[[1]], df[[2]])
    res <- 1 - (width / length)
  }
  # Related circumscribing circle
  if (type == "relCirc") {
    circle <- minCircle(x)
    areaCircle <- expanse(circle)
    area <- expanse(x)
    res <- 1-(area/areaCircle)
  }
  # Linearity index
  if (type == "lin") {
    hull <- convHull(x)
    df <- crds(as.points(hull), df = T)
    lm <- lm(y ~ x, data = df)
    res <- summary(lm)$r.squared
  }
  # North bearing of the minimum rectangle
  if (type == "bearingMinRect") {
    minRectangle <- minRect(x)
    df <- crds(as.points(minRectangle), df = T)
    cor <- cor(df[["x"]],df[["y"]])
    point1 <- slice_min(df, y)
    if (cor > 0){
      point2 <-slice_max(df, x)
    } else{
      point2 <-slice_min(df, x)
    }
    res <- bearing(point1, point2)
  }
  # Get bearing along longest axis
  if (type == "bearing") {
    convexHull <- convHull(x)
    vert <- crds(as.points(convexHull), df = T)
    dist <- as.data.frame.table(as.matrix(distance(vert, lonlat = T)), responseName = "distance") %>%
      slice_max(., distance)
    axisPoints <- vert[c(dist[[1]][[1]],dist[[2]][[1]]),]
    axisPoints <- arrange(axisPoints, desc(y))
    res <- bearing(axisPoints[2,],axisPoints[1,])
  }
  res <- as.numeric(res)
  return(res)
}








## 2. Functions for (B) Machine Learning =================================================== ####



#   ____          __  __                  _       _                    _                                      _                 
#  |  _ \   _    |  \/  |                | |     (_)                  | |                                    (_)                
#  | |_) | (_)   | \  / |   __ _    ___  | |__    _   _ __     ___    | |        ___    __ _   _ __   _ __    _   _ __     __ _ 
#  |  _ <        | |\/| |  / _` |  / __| | '_ \  | | | '_ \   / _ \   | |       / _ \  / _` | | '__| | '_ \  | | | '_ \   / _` |
#  | |_) |  _    | |  | | | (_| | | (__  | | | | | | | | | | |  __/   | |____  |  __/ | (_| | | |    | | | | | | | | | | | (_| |
#  |____/  (_)   |_|  |_|  \__,_|  \___| |_| |_| |_| |_| |_|  \___|   |______|  \___|  \__,_| |_|    |_| |_| |_| |_| |_|  \__, |
#                                                                                                                          __/ |
#                                                                                                                         |___/ 

# Function 1 : Install Packages function ============

install_and_load <- function(package_list) {
  for (pkg in package_list) {
    if (!require(pkg, character.only = TRUE)) {
      if (!is.element(pkg, installed.packages()[, "Package"])) {
        install.packages(pkg, dependencies = TRUE)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

## load dependencies of functions ##
suppressMessages(
install_and_load(
                 c(
                   "caret", "caretEnsemble",
                    "dplyr",
                   "randomForest", "recipes", "kableExtra")))


# Function 2: Pre-process the raw predictor data ============
process_data <- function(file_path, tp_value, response, vars) {
  dat <- readRDS(file_path) %>%
    filter(cell_grouping == 1 &  tp == tp_value) %>%
    select(all_of(c(response, vars))) %>%
    mutate(
      D_AOO_a = case_when(
        is.na(D_AOO_a) & rel_occ_Ncells > 0.97 ~ 2,
        TRUE ~ D_AOO_a
      ),
      mean_prob_cooccur = case_when(
        is.na(mean_prob_cooccur) & rel_occ_Ncells < 0.05 ~ 0,
        TRUE ~ mean_prob_cooccur
      )
    ) %>%
    filter(!is.na(moran)) %>%
    mutate_if(is.character, as.factor)

  # Convert specific columns to factors if needed
  dat$Habitat.Density <- as.factor(dat$Habitat.Density)
  dat$Migration <- as.factor(dat$Migration)

  return(dat)
}



# Function 3: NA summarise ============
summarize_NA <- function(dat) {
  result <- dat %>%
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "Variable",
                        values_to = "NA_Count") %>%
    filter(NA_Count > 0)

  # Print the table using kable
  kableExtra::kable(result)
}

# Funcion 4a: adjusted rfFunc helper with 5000 trees

library(caret)
library(caretEnsemble)
library(plyr)

## Default summary function (rfFuncs)
rfFuncs <- list(
  summary = # defaultSummary
    function(data, lev = NULL, model = NULL) {
      if (is.character(data$obs))
        data$obs <- factor(data$obs, levels = lev)
      postResample(data[, "pred"], data[, "obs"])
    },
  
  fit =
    function(x, y, first, last, ...) {
      loadNamespace("randomForest")
      randomForest::randomForest(x, y, importance = first, ntree = 1000, ...)
      # importance = TRUE: calculates importance for each subset of predictors
      # importance = first: calculates importance only for the model with all predictors
    },
  
  pred =
    function(object, x) {
      tmp <- predict(object, x)
      if (is.factor(object$y)) {
        out <- cbind(
          data.frame(pred = tmp),
          as.data.frame(predict(object, x, type = "prob"),
                        stringsAsFactors = TRUE))
      } else {
        out <- tmp
      }
      out
    },
  
  rank =
    function(object, x, y) {
      vimp <- caret::varImp(object, type = 1, scale = TRUE)
      # vimp <- varImp(object) # default setting
      if (is.factor(y)) {
        if (all(levels(y) %in% colnames(vimp))) {
          avImp <- apply(vimp[, levels(y), drop = TRUE], 1, mean)
          vimp$Overall <- avImp
        }
      }
      vimp <- vimp[order(vimp$importance$Overall, decreasing = TRUE), , drop = FALSE]
      if (ncol(x) == 1) {
        vimp$var <- colnames(x)
      } else {
        vimp$var <- rownames(vimp)
      }
      vimp
    },
  
  selectSize =
    function(x, metric, maximize) {
      best <- if (maximize)
        which.max(x[, metric])
      else which.min(x[, metric])
      min(x[best, "Variables"])
    },
  
  selectVar =
    function(y, size) {
      finalImp <- plyr::ddply(y[, c("Overall", "var")], .(var), function(x) mean(x$Overall, na.rm = TRUE))
      names(finalImp)[2] <- "Overall"
      finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE), ]
      as.character(finalImp$var[1:size])
    }
)


# Function 4b: rank importance in the set of all predictors
rank <- function(object, x, y) {
  vimp <- caret::varImp(object, type = 1, scale = TRUE)
  if (is.factor(y)) {
    if (all(levels(y) %in% colnames(vimp))) {
      avImp <- apply(vimp[, levels(y), drop = TRUE], 1,
                     mean)
      vimp$Overall <- avImp
    }
  }
  vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
  if (ncol(x) == 1) {
    vimp$var <- colnames(x)
  } else {
    vimp$var <- rownames(vimp)
  }
  vimp
}




## 3. Functions for plotting =================================================== ####


#  _____          ______   _                                     
#  |  __ \   _    |  ____| (_)                                    
#  | |  | | (_)   | |__     _    __ _   _   _   _ __    ___   ___ 
#  | |  | |       |  __|   | |  / _` | | | | | | '__|  / _ \ / __|
#  | |__| |  _    | |      | | | (_| | | |_| | | |    |  __/ \__ \
#  |_____/  (_)   |_|      |_|  \__, |  \__,_| |_|     \___| |___/
#                                __/ |                            
#                               |___/                             


# Function 1: Define a common theme for consistency ============
common_theme <- theme_classic() +
  theme(
    axis.line = element_line(colour = "azure4", linetype = "solid"), 
    axis.ticks = element_line(colour = "gray13"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(vjust = 0.9, hjust = 0.9),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 2),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"), 
    strip.text = element_text(size = 16, color = "black")
  )

common_theme_x_axis <- theme_classic() +
  theme(
    axis.line = element_line(colour = "azure4", linetype = "solid"), 
    axis.ticks = element_line(colour = "gray13"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 2),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"), 
    strip.text = element_text(size = 16, color = "black")
  )


## Color vector for Log Ratio =============================================== ##
my_cols_LR <- c("strong decrease (> halfing)" = "#d7191c",
                "weak decrease (< halfing)" = "#fdae61",
                "stable" = "#e0e0e0",
                "weak increase (< doubling)" = "#a6d96a",
                "strong increase (> doubling)" = "#1a9641")

## Color vector for Jaccard ================================================= ##

my_cols_J <- c("complete turnover (< 0.1)" = "#B35806",
               "strong turnover (0.1 - 0.25)" = "#F1A340",
               "strong intermediate turnover (0.25 - 0.5)" = "#FEE0B6",
               "weak intermediate turnover (0.5 - 0.75)" = "#D8DAEB",
               "weak turnover (> 0.75)" = "#998EC3",
               "stable (> 0.9)" = "#542788")


## Plot repeated cells ======================================================= ##

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(ggplot2))

plot_map <- function(grid_list, presence_data2, atlas_names, crs_list) {
  
  for (a in seq_along(grid_list)) {
    
    g1 <- grid_list[[a]]$cell1grid
      
      t_common_cells <- presence_data2 %>% 
        filter(dataset == atlas_names[[a]] & cell_grouping == "1") %>% 
        mutate(atlas_cell_repeated = as.factor(atlas_cell_repeated))
      
      t_grid <- g1 %>% st_transform(crs = crs_list[a]) %>% 
        mutate(atlas_cell_repeated = as.factor(atlas_cell_repeated))
      
      t_grid2 <- full_join(t_grid, t_common_cells) %>% 
        mutate(repeated = case_when(
          is.na(atlas_cell_repeated) ~ "never",
          atlas_cell_repeated == 1 ~ "once",
          atlas_cell_repeated == 0 ~ "never",
          atlas_cell_repeated == 2 ~ "both" )) %>% 
        mutate(repeated = as.factor(repeated)) %>% 
        reorder_levels(repeated, order = c("never", "once", "both"))
      
      # Mapping
      plot(t_grid2["repeated"], pal = hcl.colors(3, hcl.pals(type = "divergingx")[13], alpha = 0.7))
      
      
    }
  }





