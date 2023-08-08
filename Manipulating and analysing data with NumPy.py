#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

# similar to library() or require() in R

homelessness = pd.read_csv("~/Documents/Coding/Python/python_and_r_equivalents/homelessness.csv", index_col = 0)
sales = pd.read_csv("~/Documents/Coding/Python/python_and_r_equivalents/sales_subset.csv", index_col = 0)
temperatures = pd.read_csv("~/Documents/Coding/Python/python_and_r_equivalents/temperatures.csv", index_col = 0)

# similar to read.csv() or readr::read_csv() in R

# head of the data
homelessness.head()

# similar to head(homelessness) in R

# get information about the data
homelessness.info()

# similar to str(homelessness) in R

# dimensions of the data
homelessness.shape

# similar to dim(homelessness) in R

# summary statistics of the data
homelessness.describe()

# similar to summary(homelessness) in R

# see all data
homelessness.values

# R data frames are structured differently to NumPy arrays and I don't know why seeing the data this way would 
# be useful, but if you wanted row-wise chunks you could do this:
# unname(split(homelessness[-1], seq(nrow(homelessness))))

# get column names
homelessness.columns

# similar to names(homelessness) in R

# get the row indices
homelessness.index

# equivalent in R would be 1:nrow(homelessness) noting that the lowest index is 1 rather than 0

# sort the data by certain columns differentially by ascending and descending orders
homelessness.sort_values(["region", "family_members"], ascending = [True, False])

# equivalent in R (using dplyr) might be homelessness %>% arrange(region, desc(family_members))

# select certain columns
homelessness[["individuals", "state"]]

# in R: homelessness[, c("individuals", "state")]
# or using dplyr: homelessness %>% select(individuals, state)

# subsetting data
homelessness[(homelessness["family_members"] < 1000) & (homelessness["region"] == "Pacific")]

# equivalent in R:
# homelessness[homelessness$family_members < 1000 & homelessness$region == "Pacific", ]
# using dplyr:
# homelessness %>% filter(family_members < 1000, region = "Pacific") 

# subsetting using multiple values
homelessness[homelessness["region"].isin(["South Atlantic", "Mid-Atlantic"])]

# equivalent in R:
# homelessness[homelessness$region %in% c("South Atlantic", "Mid-Atlantic"), ]
# using dplyr:
# homelessness %>% filter(region %in% c("South Atlantic", "Mid-Atlantic") 

# creating a new column from other columns
homelessness["total"] = homelessness["individuals"] + homelessness["family_members"]
homelessness

# in R:
# homelessness$total <- homelessness$individuals + homelessness$family_members
# using dplyr:
# homelessness <- homelessness %>% mutate(total = individuals + family_members)

# Tying the above together in a "Combo Attack!" as DataCamp put it:
# Create indiv_per_10k col as homeless individuals per 10k state pop
homelessness["indiv_per_10k"] = 10000 * homelessness["individuals"] / homelessness["state_pop"] 

# Subset rows for indiv_per_10k greater than 20
high_homelessness = homelessness[homelessness["indiv_per_10k"] > 20]

# Sort high_homelessness by descending indiv_per_10k
high_homelessness_srt = high_homelessness.sort_values("indiv_per_10k", ascending = False)

# From high_homelessness_srt, select the state and indiv_per_10k cols
result = high_homelessness_srt[["state", "indiv_per_10k"]]

# See the result
print(result)

# a shorthand version in dplyr in R using multiple pipes and without creating extra objects:
# homelessness %>% 
#   mutate(indiv_per_10k = 10000 * individuals / state_pop) %>%
#   filter(indiv_per_10k > 20) %>%
#   arrange(desc(indiv_per_10k)) %>%
#   select(state, indiv_per_10k)

# calculate the mean or median of a column
print(sales["weekly_sales"].mean())
print(sales["weekly_sales"].median())

# equivalent in R:
# print(mean(sales$weekly_sales))
# print(median(sales$weekly_sales))

# calculate the mean or median of a column
print(sales["date"].max())
print(sales["date"].min())

# equivalent in R:
# print(max(sales$date))
# print(min(sales$date))

# find the interquartile range by defining a function and passing it to the agg() method
def iqr(column):
    return column.quantile(0.75) - column.quantile(0.25)
    
sales["temperature_c"].agg(iqr)

# in R there is an inbuilt function for this:
# IQR(sales$temperature_c)

# Sort sales_1_1 by date
sales_1_1 = sales[(sales["department"] == 1) & (sales["store"] == 1)]
sales_1_1 = sales_1_1.sort_values("date")

# Get the cumulative sum of weekly_sales, add as cum_weekly_sales col
sales_1_1["cum_weekly_sales"] = sales_1_1["weekly_sales"].cumsum()

# Get the cumulative max of weekly_sales, add as cum_max_sales col
sales_1_1["cum_max_sales"] = sales_1_1["weekly_sales"].cummax()

# See the columns you calculated
sales_1_1[["date", "weekly_sales", "cum_weekly_sales", "cum_max_sales"]]

# the equivalent using R and dplyr:
# sales %>% 
#  filter(store == 1, department == 1) %>%
#  arrange(date) %>%
#  mutate(cum_weekly_sales = cumsum(weekly_sales)) %>%
#  mutate(cum_max_sales = cummax(weekly_sales)) %>%
#  select(date, weekly_sales, cum_weekly_sales, cum_max_sales)

# dropping duplicates
sales.drop_duplicates(["store", "type"])

# using dplyr in R:
# sales %>% filter(!duplicated(store, type))

# using dplyr in R and also removing the introduced artifact of the first column:
# sales %>% filter(!duplicated(store, type)) %>% select(-X)

# dropping duplicates and subsetting
sales[sales["is_holiday"] == True].drop_duplicates("date")

# in R and dplyr, note that the column is_holiday will be treated as a character vector as logical values
# must be all caps and encoded as logical
# sales %>% filter(is_holiday == "True", !duplicated(date))
#
# to correct this:
# sales %>% 
#   mutate(is_holiday = as.logical(toupper(is_holiday))) %>%
#   filter(is_holiday, !duplicated(date))

# counting
store_types = sales.drop_duplicates(subset = ["store", "type"])
store_types["type"].value_counts()

# in R and dplyr:
# sales %>% 
#   filter(!duplicated(store, type)) %>%
#   count(type)
# however it does produce a warning because of two variables being in duplicated(). upon consulting
# about not getting the warning ChatGPT suggested the code could be written thus:
# sales[!duplicated(sales[c("store", "type")]), ] %>% count(type)

# counting, sorting and normalising
store_depts = sales.drop_duplicates(subset = ["store", "department"])
store_depts["department"].value_counts(sort = True, normalize = True)

# in R and dplyr without a warning and returned as a vector (which was returned in both this and the 
# above examples):
# sales %>%
#   distinct(store, department, .keep_all = TRUE) %>%
#   count(department, sort = TRUE) %>%
#   mutate(n = n / sum(n)) %>%
#   pull()

# Calc total weekly sales
sales_all = sales["weekly_sales"].sum()

# Subset for type A stores, calc total weekly sales
sales_A = sales[sales["type"] == "A"]["weekly_sales"].sum()

# Subset for type B stores, calc total weekly sales
sales_B = sales[sales["type"] == "B"]["weekly_sales"].sum()

# Subset for type C stores, calc total weekly sales
sales_C = sales[sales["type"] == "C"]["weekly_sales"].sum()

# Get proportion for each type
sales_propn_by_type = [sales_A, sales_B, sales_C] / sales_all
print(sales_propn_by_type)

# in R this can be done with remarkably similar code using = to assign and keeping some column calls
# as quoted names in square brackets
# sales_all = sum(sales["weekly_sales"])
# sales_A = sum(sales[sales$type == "A", ]["weekly_sales"])
# sales_B = sum(sales[sales$type == "B", ]["weekly_sales"])
# sales_C = sum(sales[sales$type == "C", ]["weekly_sales"])
# sales_propn_by_type = c(sales_A, sales_B, sales_C) / sales_all
# print(sales_propn_by_type)

# using shorthand dplyr code we achieve a slightly different result as it turns out there is not a "C"-type
# sales %>% 
#   group_by(type) %>%
#   mutate(type_sales = sum(weekly_sales)) %>%
#   distinct(type_sales) %>%
#   ungroup() %>% 
#   mutate(sales_propn_by_type = type_sales / sum(type_sales)) %>% 
#   pull(sales_propn_by_type)

# grouping
sales.groupby("type")["weekly_sales"].sum()

#with R and dplyr returning a data frame and vector:
# sales %>% 
#   group_by(type) %>% 
#   summarise(sum = sum(weekly_sales))
    
# sales %>% 
#   group_by(type) %>% 
#   summarise(sum = sum(weekly_sales)) %>%
#   pull()

# grouping and aggregate statistics
sales.groupby("type")[["unemployment", "fuel_price_usd_per_l"]].agg([np.min, np.max, np.mean, np.median])

# to my knowledge there is no real way to produce a similar output in R. the analysis would need to be run
# separately (as in the example), or simply name the output columns according to the group. Python:1, R:0?
# sales %>% 
#   group_by(type) %>%
#   summarise(amin = min(unemployment),
#             amax = max(unemployment),
#             mean = mean(unemployment),
#             median = median(unemployment))

# sales %>% 
#   group_by(type) %>%
#   summarise(amin = min(fuel_price_usd_per_l),
#             amax = max(fuel_price_usd_per_l),
#             mean = mean(fuel_price_usd_per_l),
#             median = median(fuel_price_usd_per_l))

# pivot tables can also be used to perform grouped calculations and this instance uses two variables
# which here are "department" and "type"
sales.pivot_table(values = "weekly_sales", index = "department", columns = "type", fill_value = 0, margins = True)

# in R, we might do this with tidyr::spread() or tidyr::pivot_wider() I couldn't do it off the cuff.
# ChatGPT also tried but failed twice, so I'll leave it for now!

# indexing
# Look at temperatures
print(temperatures)

# Set the index of temperatures to city
temperatures_ind = temperatures.set_index("city")

# Look at temperatures_ind
print(temperatures_ind)

# Reset the temperatures_ind index, keeping its contents
print(temperatures_ind.reset_index())

# Reset the temperatures_ind index, dropping its contents
print(temperatures_ind.reset_index(drop = True))

# in R, duplicate row names are not allowed so a command like the following fails:
# temperatures %>% tibble::column_to_rownames(var = "city")

