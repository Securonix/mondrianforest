## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE------------------------------------------------------
library(mondrianforest)
library(purrr)

## ------------------------------------------------------------------------
census_data <- readRDS("../sample_data/adult.3000.rds")

train <- census_data[1:2500, ]
test <- census_data[2501:3000, ]

head(train)

## ---- cache = TRUE-------------------------------------------------------
mf <- mondrian_forest(train[1:250, ], y_col_num = 15, lambda = 8)

# Class Probabilities:
predict(mf, test[1:5, ], type = "prob")

# Confusion Matrix:
table(test$X15, predict(mf, test, type = "class"), dnn = list("Actual", "Predicted"))

## ---- cache = TRUE-------------------------------------------------------
mf_extended <- extend_mondrian_forest(mf, train[251:500, ])

# Class Probabilities:
predict(mf_extended, test[1:5, ], type = "prob")

# Confusion Matrix:
table(test$X15, predict(mf_extended, test, type = "class"), dnn = list("Actual", "Predicted"))

## ------------------------------------------------------------------------
split_feature <- function(df, column_num, split_set, split_set_description, complement_description) {
  original <- df[[column_num]]
  split_a <- split_b <- original
  split_a[!(original %in% split_set)] <- complement_description # non-split_set variables aggregated into complement_description
  split_b[original %in% split_set] <- split_set_description
  df[[paste(names(df)[column_num], "a",sep = "_")]] <- split_a
  df[[paste(names(df)[column_num], "b",sep = "_")]] <- split_b
  df[[column_num]] <- NULL
  df <- dplyr::select(df, 1:(column_num-1), ncol(df)-1, ncol(df), column_num:(ncol(df)-2)) # Reorder so data_frame looks like original
  df
}

census_data <- split_feature(census_data, column_num = 4, split_set = c("10th", "11th", "12th", "1st-4th", "5th-6th", "7th-8th", "9th"
), "Grade_School", "Higher_Ed")

## ---- cache = TRUE-------------------------------------------------------
# Numeric variables are interpreted as-is
census_data <- map_df(census_data, function(x) if (is.numeric(x)) (x - min(x))/(max(x)-min(x)) else x)

# Categorical variables are defaulted so each dummy level has a range of 1
categorical_dummy_columns <- map(census_data[, -ncol(census_data)], function(x) if (!is.numeric(x)) length(unique(x)) - 1 else numeric(0)) %>% unlist()

# number of implicit dummy columns for each categorical variable
categorical_dummy_columns

## ---- cache = TRUE-------------------------------------------------------
train <- census_data[1:2500, ]
test <- census_data[2501:3000, ]

mf_scaled <- mondrian_forest(train[1:500, ], y_col_num = 16, lambda = 8, f_scale = 1/categorical_dummy_columns)

## ---- cache = TRUE-------------------------------------------------------
pryr::object_size(mf_extended)
pryr::object_size(mf_scaled)

table(test$X15, predict(mf_scaled, test, type = "class"))

