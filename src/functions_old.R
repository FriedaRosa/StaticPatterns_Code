# Functions for Machine Learning Assignment #

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
install_and_load(
                 c(
                   "caret", "caretEnsemble",
                   "plyr", "dplyr",
                   "randomForest", "recipes", "kableExtra"))


# Function 2: Pre-process the raw predictor data ============
process_data <- function(file_path, tp_value, response, vars) {
  dat <- readRDS(file_path) %>%
    filter(cell_grouping == 1 & exclude == 0 & tp == tp_value) %>%
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



##############
library(caret)
library(caretEnsemble)
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
      randomForest::randomForest(x, y, importance = first, ntree = 5000, ...)
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
      finalImp <- ddply(y[, c("Overall", "var")], .(var), function(x) mean(x$Overall, na.rm = TRUE))
      names(finalImp)[2] <- "Overall"
      finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE), ]
      as.character(finalImp$var[1:size])
    }
)


## Custom summary function for randomForest (simplified)
rfRFE1 <- list(
  summary = defaultSummary,
  fit = function(x, y, first, last, ...) {
    library(randomForest)
    randomForest(x, y, importance = first, ...) # importance = first indicates that the first model is the complete one (with all predictors)
  },
  pred = function(object, x) predict(object, x),
  rank = function(object, x, y) {
    vimp <- caret::varImp(object, type = 1, scale = TRUE)
    vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
    vimp$var <- rownames(vimp)
    vimp
  },
  selectSize = pickSizeBest,
  selectVar = pickVars
)


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
