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

