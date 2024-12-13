### DEPENDENT VARIABLE: MEDIAN RENT

# Libraries repository

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(googlesheets4)
library(sf)
library(httr)
library(tidycensus)
library(tigris)
library(purrr)
library(fixest)
library(stargazer)
library(sandwich)
library(lmtest)
library(broom)

# Clearing

rm(list = ls())

# Importing live data

gs4_auth()
sheet_url <- "https://docs.google.com/spreadsheets/d/183amxMrOQILhY3qLbOmf34fLknf0GYmpe7TsAr3Rasw/edit?gid=0#gid=0"
data <- read_sheet(sheet_url)

# Creating the histogram of change in betweenness centrality

ggplot(data, aes(x = diffincent)) +
  geom_histogram(bins = 25, fill = "white", color = "black") +
  labs(title = "Histogram of Difference in Betweenness Centrality Pre-Post Silver Line Phase I",
       x = "Difference in Betweenness Centrality Pre-Post", y = "Frequency") +
  theme_minimal()

# Function to get ACS data for a given year

get_acs_data_for_year <- function(year) {
  safely(function() {
    acs_data <- get_acs(
      geography = "tract",
      variables = c(
        median_rent = "B25064_001",          # Median gross rent
        housing_stock = "B25001_001",        # Total housing units
        white = "B02001_002",                # White population
        black = "B02001_003",                # Black or African American population
        native = "B02001_004",               # American Indian and Alaska Native population
        asian = "B02001_005",                # Asian population
        pacificislander = "B02001_006",      # Native Hawaiian and Other Pacific Islander population
        multiracial = "B02001_008",          # Two or more races
        in_poverty = "B17001_002",           # Population below the poverty level
        housing_cost_burden_30_34 = "B25070_007",  # Households paying 30-34.9% of income on housing
        housing_cost_burden_35_more = "B25070_008", # Households paying 35% or more on housing
        geographical_mobility = "B07003_003",      # Moved within the same county in the last year
        median_home_value = "B25077_001",     # Median home value
        vacant_housing = "B25002_003",        # Vacant housing units
        median_year_built = "B25035_001",     # Median year structure built
        homeownership_rate = "B25003_002",    # Owner-occupied housing units
        median_income = "B19013_001",         # Median household income
        unemployment_rate = "B23025_005",     # Unemployed population
        high_school_grads = "B15003_017",     # Population with high school diploma
        college_grads = "B15003_022",         # Population with Bachelor's degree or higher
        commute_time = "B08303_001",          # Mean travel time to work
        labor_force_participation = "B23025_002",  # Population in labor force
        housing_units_by_type = "B25024_001",  # Housing units by type
        population = "B01001_001"             # Total population (for population density)
      ),
      year = year,
      state = c("MD", "VA", "DC"),
      survey = "acs5",
      geometry = FALSE
    ) %>%
      select(GEOID, variable, estimate) %>%
      spread(variable, estimate)
    
    return(acs_data)
  })()
}

# Retrieving ACS data for 2014 (controls) and 2019 (dependent variable)

acs_2014 <- get_acs_data_for_year(2014)$result %>%
  rename_with(~ paste0(.x, "_control"), -GEOID)
acs_2019 <- get_acs_data_for_year(2019)$result %>%
  rename_with(~ paste0(.x, "_dep"), -GEOID)

# Data processing

if ("Coordinates" %in% names(data)) {
  data <- data %>%
    separate(Coordinates, into = c("latitude", "longitude"), sep = ",", convert = TRUE)
} else {
  stop("Coordinates column not found in data.")
}

data <- data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

tracts_md <- tracts(state = "MD", cb = TRUE, year = 2020)
tracts_va <- tracts(state = "VA", cb = TRUE, year = 2020)
tracts_dc <- tracts(state = "DC", cb = TRUE, year = 2020)

tracts <- rbind(tracts_md, tracts_va, tracts_dc)
data <- st_transform(data, crs = st_crs(tracts))
data_with_geoid <- st_join(data, tracts["GEOID"], join = st_within)

# Merging ACS data for controls (2014) with dependent variable (2019)

final_data <- data_with_geoid %>%
  left_join(acs_2014, by = "GEOID") %>%
  left_join(acs_2019, by = "GEOID") %>%
  mutate(
    median_rent_lag = median_rent_control,  # Lagged variable from 2014
    median_rent = median_rent_dep           # Dependent variable from 2019
  ) %>%
  mutate(
    controls_complete = if_else(
      !is.na(median_income_control) & !is.na(housing_cost_burden_35_more_control) &
        !is.na(unemployment_rate_control) & !is.na(college_grads_control) &
        !is.na(median_home_value_control),
      1, 0
    )
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Regression models

model1 <- lm(median_rent ~ diffincent + median_rent_lag + controls_complete, data = final_data)
model2 <- lm(median_rent ~ diffincent + median_rent_lag + housing_stock_control + controls_complete, data = final_data)
model3 <- lm(median_rent ~ diffincent + median_rent_lag + housing_stock_control + median_income_control + controls_complete, data = final_data)
model4 <- lm(median_rent ~ diffincent + median_rent_lag + housing_stock_control + median_income_control + controls_complete + asian_control + black_control +
               college_grads_control + commute_time_control + geographical_mobility_control + high_school_grads_control +
               homeownership_rate_control + housing_cost_burden_30_34_control + housing_cost_burden_35_more_control +
               in_poverty_control + labor_force_participation_control + median_home_value_control +
               median_year_built_control + multiracial_control + native_control + pacificislander_control + population_control +
               unemployment_rate_control + vacant_housing_control + white_control + LineRed + LineGreen +
               LineYellow + LineBlue + LineOrange + LineSilver + InBeltway + FurtherRail +
               Parking, data = final_data)


omit_vars <- c("asian_control", "black_control", "college_grads_control", "commute_time_control",
               "geographical_mobility_control", "high_school_grads_control", "homeownership_rate_control",
               "housing_cost_burden_30_34_control", "housing_cost_burden_35_more_control", "in_poverty_control",
               "labor_force_participation_control", "median_home_value_control",
               "median_year_built_control", "multiracial_control", "native_control", "pacificislander_control",
               "population_control", "unemployment_rate_control", "vacant_housing_control", "white_control", "LineRed", "LineGreen",
               "LineYellow", "LineBlue", "LineOrange", "LineSilver", "InBeltway", "FurtherRail", "Parking")

# Summary table

stargazer(
  model1, model2, model3, model4,
  type = "text",
  column.labels = c("OLS with lagged median rent", "+ Housing Stock", "+ Median Income", "+ Other Controls"),
  omit = omit_vars,
  keep.stat = c("n", "rsq", "adj.rsq", "f"),
  title = "Centrality Regression Results with Lagged Median Rent",
  add.lines = list(c("Other controls", "No", "No", "No", "Yes"))
)

# Full summary table

stargazer(
  model1, model2, model3, model4,
  type = "text",
  column.labels = c("OLS with lagged median rent", "+ Housing Stock", "+ Median Income", "+ Other Controls"),
  keep.stat = c("n", "rsq", "adj.rsq", "f"),
  title = "Centrality Regression Results with Lagged Median Rent",
  add.lines = list(c("Other controls", "No", "No", "No", "Yes"))
)


---
  
### DEPENDENT VARIABLE: HOUSING STOCK
  
# Libraries repository
  
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(googlesheets4)
library(sf)
library(httr)
library(tidycensus)
library(tigris)
library(purrr)
library(fixest)
library(stargazer)
library(sandwich)
library(lmtest)
library(broom)

# Clearing

rm(list = ls())

# Importing live data

gs4_auth()
sheet_url <- "https://docs.google.com/spreadsheets/d/183amxMrOQILhY3qLbOmf34fLknf0GYmpe7TsAr3Rasw/edit?gid=0#gid=0"
data <- read_sheet(sheet_url)

# Function to get ACS data for a given year

get_acs_data_for_year <- function(year) {
  safely(function() {
    acs_data <- get_acs(
      geography = "tract",
      variables = c(
        median_rent = "B25064_001",          # Median gross rent
        housing_stock = "B25001_001",        # Total housing units
        white = "B02001_002",                # White population
        black = "B02001_003",                # Black or African American population
        native = "B02001_004",               # American Indian and Alaska Native population
        asian = "B02001_005",                # Asian population
        pacificislander = "B02001_006",      # Native Hawaiian and Other Pacific Islander population
        multiracial = "B02001_008",          # Two or more races
        in_poverty = "B17001_002",           # Population below the poverty level
        housing_cost_burden_30_34 = "B25070_007",  # Households paying 30-34.9% of income on housing
        housing_cost_burden_35_more = "B25070_008", # Households paying 35% or more on housing
        geographical_mobility = "B07003_003",      # Moved within the same county in the last year
        median_home_value = "B25077_001",     # Median home value
        vacant_housing = "B25002_003",        # Vacant housing units
        median_year_built = "B25035_001",     # Median year structure built
        homeownership_rate = "B25003_002",    # Owner-occupied housing units
        median_income = "B19013_001",         # Median household income
        unemployment_rate = "B23025_005",     # Unemployed population
        high_school_grads = "B15003_017",     # Population with high school diploma
        college_grads = "B15003_022",         # Population with Bachelor's degree or higher
        commute_time = "B08303_001",          # Mean travel time to work
        labor_force_participation = "B23025_002",  # Population in labor force
        housing_units_by_type = "B25024_001",  # Housing units by type
        population = "B01001_001"             # Total population (for population density)
      ),
      year = year,
      state = c("MD", "VA", "DC"),
      survey = "acs5",
      geometry = FALSE
    ) %>%
      select(GEOID, variable, estimate) %>%
      spread(variable, estimate)
    
    return(acs_data)
  })()
}

# Retrieving ACS data for 2014 (controls) and 2019 (dependent variable)

acs_2014 <- get_acs_data_for_year(2014)$result %>%
  rename_with(~ paste0(.x, "_control"), -GEOID) 
acs_2019 <- get_acs_data_for_year(2019)$result %>%
  rename_with(~ paste0(.x, "_dep"), -GEOID)

# Data processing

data <- data %>%
  separate(Coordinates, into = c("latitude", "longitude"), sep = ",", convert = TRUE)

data <- data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

tracts_md <- tracts(state = "MD", cb = TRUE, year = 2020)
tracts_va <- tracts(state = "VA", cb = TRUE, year = 2020)
tracts_dc <- tracts(state = "DC", cb = TRUE, year = 2020)

tracts <- rbind(tracts_md, tracts_va, tracts_dc)
data <- st_transform(data, crs = st_crs(tracts))
data_with_geoid <- st_join(data, tracts["GEOID"], join = st_within)

# Merging ACS data for controls (2014) and dependent variable (2019)

final_data <- data_with_geoid %>%
  left_join(acs_2014, by = "GEOID") %>%
  left_join(acs_2019, by = "GEOID") %>%
  mutate(
    housing_stock_lag = housing_stock_control,
    housing_stock = housing_stock_dep 
  ) %>%
  mutate(
    controls_complete = if_else(
      !is.na(median_income_control) & !is.na(housing_cost_burden_35_more_control) &
        !is.na(unemployment_rate_control) & !is.na(college_grads_control) &
        !is.na(median_home_value_control),
      1, 0
    )
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Regression models

model1 <- lm(housing_stock ~ diffincent + housing_stock_lag + controls_complete, data = final_data)
model2 <- lm(housing_stock ~ diffincent + housing_stock_lag + median_income_control + controls_complete, data = final_data)
model3 <- lm(housing_stock ~ diffincent + housing_stock_lag  + median_income_control + population_control + controls_complete, data = final_data)
model4 <- lm(housing_stock ~ diffincent + housing_stock_lag + controls_complete + median_income_control + population_control + asian_control + black_control + college_grads_control + commute_time_control +
               geographical_mobility_control + high_school_grads_control + homeownership_rate_control +
               housing_cost_burden_30_34_control + housing_cost_burden_35_more_control + in_poverty_control +
               labor_force_participation_control + median_home_value_control +
               median_year_built_control + multiracial_control + native_control + pacificislander_control +
               unemployment_rate_control + vacant_housing_control + white_control + LineRed + LineGreen + LineYellow + LineBlue + LineOrange +
               LineSilver + InBeltway + FurtherRail + Parking, data = final_data)

omit_vars <- c("asian_control", "black_control", "college_grads_control", "commute_time_control",
               "geographical_mobility_control", "high_school_grads_control", "homeownership_rate_control",
               "housing_cost_burden_30_34_control", "housing_cost_burden_35_more_control", "in_poverty_control",
               "labor_force_participation_control", "median_home_value_control",
               "median_year_built_control", "multiracial_control", "native_control", "pacificislander_control",
               "unemployment_rate_control", "vacant_housing_control", "white_control", "AvgDailyTappedEntries",
               "AvgDailyNonTappedEntries", "LineRed", "LineGreen", "LineYellow", "LineBlue", "LineOrange",
               "LineSilver", "InBeltway", "FurtherRail", "Parking")

# Summary table

stargazer(
  model1, model2, model3, model4,
  type = "text",
  column.labels = c("OLS with lagged housing stock", "+ Median Income", "+ Population", "+ Other Controls"),
  omit = omit_vars,  # Hide control variables only in Model 4
  keep.stat = c("n", "rsq", "adj.rsq", "f"),
  title = "Housing Stock Regression Results with Lagged Variables",
  add.lines = list(c("Other controls", "No", "No", "No", "Yes"))
)

# Full summary table

stargazer(
  model1, model2, model3, model4,
  type = "text",
  column.labels = c("OLS with lagged housing stock", "+ Median Income", "+ Population", "+ Other Controls"),
  keep.stat = c("n", "rsq", "adj.rsq", "f"),
  title = "Housing Stock Regression Results with Lagged Variables",
  add.lines = list(c("Other controls", "No", "No", "No", "Yes"))
)

