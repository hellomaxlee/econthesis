### ALL DiD ANALYSIS

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
library(gt)
library(knitr)
library(kableExtra)
library(modelsummary)

# Clearing

rm(list = ls())

# Importing live data

gs4_auth()
sheet_url <- "https://docs.google.com/spreadsheets/d/1_z3MMwhvb9Y5W7JCCXpNasd27RRPwsIU7Jvz05HA8-Q/edit?gid=0#gid=0"
data <- read_sheet(sheet_url)
head(data)

data <- data %>% 
  separate(Coordinates, into = c("latitude", "longitude"), sep = ",", convert = TRUE)

# Spatial definition of treated vs. control groups

stations_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

treated_stations <- stations_sf %>% filter(Treated == 1)
control_stations <- stations_sf %>% filter(Treated == 0)

treated_buffer <- treated_stations %>%
  st_transform(crs = 32618) %>% 
  st_buffer(dist = 1609.34*.75) %>%
  st_transform(crs = 4326)

control_buffer <- control_stations %>%
  st_transform(crs = 32618) %>%
  st_buffer(dist = 1609.34*.75) %>% 
  st_transform(crs = 4326)

tracts_md <- tracts(state = "MD", cb = TRUE, year = 2020)
tracts_va <- tracts(state = "VA", cb = TRUE, year = 2020)
tracts_dc <- tracts(state = "DC", cb = TRUE, year = 2020)

tracts <- rbind(tracts_md, tracts_va, tracts_dc)

tracts <- st_transform(tracts, crs = st_crs(treated_buffer))

tracts_treated <- st_intersects(tracts, treated_buffer, sparse = FALSE)
tracts_control <- st_intersects(tracts, control_buffer, sparse = FALSE)

tracts <- tracts %>%
  mutate(treatment_status = case_when(
    apply(tracts_treated, 1, any) ~ "Treated",
    apply(tracts_control, 1, any) ~ "Control",
    TRUE ~ "Not in the Sample"
  ))


ggplot() +
  geom_sf(data = tracts, aes(fill = treatment_status), color = "black") +
  geom_sf(data = treated_buffer, fill = NA, color = "blue", linetype = "dashed") +
  geom_sf(data = control_buffer, fill = NA, color = "green", linetype = "dashed") +
  geom_sf(data = stations_sf, aes(color = factor(Treated)), size = 2) +
  labs(title = "Census Tracts within 1 Mile Buffer of Treated and Control Group 1 Stations") +
  coord_sf(xlim = c(-77.5, -76.6), ylim = c(38.6, 39.3), expand = FALSE)
  
# Data cleaning and transformation

tracts <- rbind(tracts_md, tracts_va, tracts_dc)

tracts <- st_transform(tracts, crs = st_crs(treated_buffer))

stations_buffer <- rbind(
  treated_buffer %>% mutate(Treated = 1),
  control_buffer %>% mutate(Treated = 0)
)

stations_buffer$station_id <- seq_len(nrow(stations_buffer))

tracts_stations_intersect <- st_intersection(tracts, stations_buffer)

tracts_stations_intersect <- tracts_stations_intersect %>%
  mutate(intersection_area = st_area(geometry))

tracts_stations <- tracts_stations_intersect %>%
  group_by(GEOID) %>%
  slice_max(intersection_area, n = 1) %>%
  ungroup() %>%
  select(GEOID, Treated, station_associated = station_id)

tracts_stations <- tracts_stations %>%
  mutate(Treated = ifelse(Treated == 1, 1, 0))

years <- data.frame(Year = 2009:2022)

tracts_years <- tracts_stations %>%
  st_drop_geometry() %>%   
  crossing(years)

tracts_years <- tracts_years %>%
  mutate(Post = ifelse(Year >= 2015, 1, 0))


# ACS API

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
      survey = "acs5",  # Use 5-year ACS data
      geometry = FALSE
    ) %>%
      # Select relevant columns and reshape the data
      select(GEOID, variable, estimate) %>%
      spread(variable, estimate) %>%
      # Rename columns to include the year in the variable names
      rename_with(~ paste0(., "_", year), -GEOID) %>%
      # Add a column indicating the data year
      mutate(Year = year)  # Add Year column here for the join
    
    return(acs_data)
  })()
}

years <- 2009:2022
acs_data_all_years <- map(years, ~ get_acs_data_for_year(.x))

successful_acs_data <- map(acs_data_all_years, "result") %>%
  discard(is.null)

acs_data_combined <- bind_rows(successful_acs_data)

merged_data <- tracts_years %>%
  left_join(acs_data_combined, by = c("GEOID", "Year"))

merged_data <- merged_data %>%
  mutate(year_num = Year - 2008)

# Dataset merging and transformation

did_data <- merged_data %>%
  mutate(
    Treated_Post = Treated * Post  # Create the interaction term in mutate()
  ) %>%
  select(
    GEOID,                    # Census Tract ID
    Treated,                  # Treated status
    year_num,                 # Year number (numeric: 1 to 14)
    Year,                     # Actual year
    Post,                     # Post-treatment indicator
    Treated_Post,             # Interaction term (Treated * Post)
    station_associated,       # Associated station
    starts_with("median_rent"), # Median rent by year
    starts_with("housing_stock"),
    starts_with("white"),
    starts_with("black"),
    starts_with("native"),
    starts_with("asian"),
    starts_with("pacificislander"),
    starts_with("multiracial"),
    starts_with("in_poverty"),
    starts_with("housing_cost_burden_30_34"),
    starts_with("housing_cost_burden_35_more"),
    starts_with("geographical_mobility"),
    starts_with("median_home_value"),
    starts_with("vacant_housing"),
    starts_with("median_year_built"),
    starts_with("homeownership_rate"),
    starts_with("median_income"),
    starts_with("unemployment_rate"),
    starts_with("high_school_grads"),
    starts_with("college_grads"),
    starts_with("commute_time"),
    starts_with("labor_force_participation"),
    starts_with("housing_units_by_type"),
    starts_with("population")
  )

acs_data_long <- acs_data_combined %>%
  pivot_longer(
    cols = -c(GEOID, Year),
    names_to = "variable",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  )

merged_data <- tracts_years %>%
  left_join(acs_data_long, by = c("GEOID", "Year"))

variables <- c("median_rent", "housing_stock", "white", "black", "native", "asian",
               "pacificislander", "multiracial", "in_poverty", "housing_cost_burden_30_34",
               "housing_cost_burden_35_more", "geographical_mobility", "median_home_value",
               "vacant_housing", "median_year_built", "homeownership_rate", "median_income",
               "unemployment_rate", "high_school_grads", "college_grads", "commute_time",
               "labor_force_participation", "housing_units_by_type", "population")

for (variable in variables) {
  merged_data[[variable]] <- merged_data %>%
    select(starts_with(variable)) %>%
    reduce(coalesce)  # Coalesce to create a single column for each variable
}

merged_data <- merged_data %>%
  select(-matches("_(2009|2010|2011|2012|2013|2014|2015|2016|2017|2018|2019|2020|2021|2022)$"))

did_data_final <- merged_data %>%
  mutate(
    year_num = Year - 2008,
    Post = ifelse(Year >= 2015, 1, 0),
    Treated_Post = Treated * Post
  ) %>%
  select(GEOID, Treated, year_num, Year, Post, Treated_Post, station_associated, everything())

---

# Summary statistics with clustered standard errors
  
rent_summary <- did_data_final %>%
  group_by(year_num, Treated) %>%
  summarise(
    mean_median_rent = mean(median_rent, na.rm = TRUE),
    se_median_rent = {
      grouped_data <- cur_data()
      num_clusters <- n_distinct(grouped_data$station_associated)
      n <- n()
      cluster_variance <- grouped_data %>%
        group_by(station_associated) %>%
        summarise(cluster_mean = mean(median_rent, na.rm = TRUE)) %>%
        summarise(var = sum((cluster_mean - mean_median_rent)^2) / (num_clusters - 1))
      sqrt(cluster_variance$var / n)
    }
  ) %>%
  ungroup()

rent_summary_wide <- rent_summary %>%
  pivot_wider(
    names_from = Treated,
    values_from = c(mean_median_rent, se_median_rent),
    names_prefix = "Treated_"
  )

names(rent_summary_wide) <- c(
  "year_num", "Treated_0_mean_median_rent", "Treated_1_mean_median_rent",
  "Treated_0_se_median_rent", "Treated_1_se_median_rent"
)

diff_summary <- rent_summary_wide %>%
  filter(year_num >= 4) %>% 
  mutate(
    diff = Treated_1_mean_median_rent - Treated_0_mean_median_rent, 
    se_diff = sqrt(Treated_1_se_median_rent^2 + Treated_0_se_median_rent^2)
  )

year_labels <- seq(-2, 8, length.out = nrow(diff_summary))

reference_diff <- diff_summary$diff[diff_summary$year_num == 6]

ggplot(diff_summary, aes(x = year_num, y = diff_centered)) +
    geom_line(color = "black", size = 1.2) +
    geom_point(color = "black", size = 2) +
    geom_ribbon(
      aes(
        ymin = diff_centered - se_diff,
        ymax = diff_centered + se_diff
      ),
      fill = "gray", alpha = 0.3
    ) +
    geom_vline(xintercept = 6, linetype = "dashed", color = "red") +
    labs(
      title = "Median Rent: Treated Minus Control Group",
      x = "Period",
      y = "Difference in Median Rent (Treated - Control)"
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = diff_summary$year_num, labels = year_labels) +
    geom_hline(yintercept = 0, linetype = "dotted") + 
    theme(legend.position = "none")

---
  
# Housing stock event study

housing_summary <- did_data_final %>%
  group_by(year_num, Treated) %>%
  summarise(
    mean_housing_stock = mean(housing_stock, na.rm = TRUE),
    se_housing_stock = {
      grouped_data <- cur_data()
      num_clusters <- n_distinct(grouped_data$station_associated)
      n <- n()
      cluster_variance <- grouped_data %>%
        group_by(station_associated) %>%
        summarise(cluster_mean = mean(housing_stock, na.rm = TRUE)) %>%
        summarise(var = sum((cluster_mean - mean_housing_stock)^2) / (num_clusters - 1))
      sqrt(cluster_variance$var / n)
    }
  ) %>%
  ungroup()

housing_summary_wide <- housing_summary %>%
  pivot_wider(
    names_from = Treated,
    values_from = c(mean_housing_stock, se_housing_stock),
    names_prefix = "Treated_"
  )

names(housing_summary_wide) <- c(
  "year_num", "Treated_0_mean_housing_stock", "Treated_1_mean_housing_stock",
  "Treated_0_se_housing_stock", "Treated_1_se_housing_stock"
)

diff_summary <- housing_summary_wide %>%
  filter(year_num >= 4) %>%
  mutate(
    diff = Treated_1_mean_housing_stock - Treated_0_mean_housing_stock,
    se_diff = sqrt(Treated_1_se_housing_stock^2 + Treated_0_se_housing_stock^2)
  )

year_labels <- seq(-2, 8, length.out = nrow(diff_summary))

ggplot(diff_summary, aes(x = year_num, y = diff_centered)) +
  geom_line(color = "black", size = 1.2) +
  geom_point(color = "black", size = 2) +
  geom_ribbon(
    aes(
      ymin = diff_centered - se_diff,
      ymax = diff_centered + se_diff
    ),
    fill = "gray", alpha = 0.3
  ) +  # Add standard error ribbon
  geom_vline(xintercept = 6, linetype = "dashed", color = "red") +
  labs(
    title = "Median Rent: Treated Minus Control Group",
    x = "Period",
    y = "Difference in Median Rent (Treated - Control)"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = diff_summary$year_num, labels = year_labels) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(legend.position = "none")

---
  
# Separated post periods - median rent
  
did_data_final <- did_data_final %>%
  mutate(
    Post_Period_1 = ifelse(year_num >= 7 & year_num <= 10, 1, 0),
    Post_Period_2 = ifelse(year_num >= 11 & year_num <= 14, 1, 0),
    Treated_Post_Period_1 = Treated * Post_Period_1,
    Treated_Post_Period_2 = Treated * Post_Period_2
  )

model_comb_summary <- summary(model_comb)
model_comb_lm <- lm(median_rent ~ Treated + Post_Period_1 + Post_Period_2 + Treated_Post_Period_1 + Treated_Post_Period_2, 
                    data = did_data_final)

model_comb_lm$coefficients <- coef(model_comb_summary)
model_comb_lm$vcov <- vcov(model_comb_summary)

stargazer(model_comb_lm, type = "text",
          title = "Difference-in-Differences Estimates for Post Periods 1 and 2",
          dep.var.labels = "Median Rent",
          covariate.labels = c("Treated", "Post Period 1", "Post Period 2", 
                               "Treated:Post Period 1", "Treated:Post Period 2"),
          column.labels = c("Model 1"),
          omit.stat = c("f", "ser"),
          no.space = TRUE,
          digits = 2)

---

# Different model specifications - median rent

# Model 1: Difference-in-Differences (DiD) without controls
model_1 <- feols(
  median_rent ~ Treated + Post + Treated_Post,
  data = did_data_final,
  cluster = ~station_associated
  )

# Model 2: Fixed Effects with no controls
model_2 <- feols(
  median_rent ~ Treated + Post + Treated_Post | GEOID + Year,
  data = did_data_final,
  cluster = ~station_associated
)

# Model 3: Difference-in-Differences (DiD) with controls
model_3 <- feols(
  median_rent ~ Treated + Post + Treated_Post +
    housing_stock + white + black + native + asian +
    pacificislander + multiracial + in_poverty + housing_cost_burden_30_34 +
    housing_cost_burden_35_more + geographical_mobility + median_home_value +
    vacant_housing + median_year_built + homeownership_rate + median_income +
    unemployment_rate + high_school_grads + college_grads + commute_time +
    labor_force_participation + population,
  data = did_data_final,
  cluster = ~station_associated
)

# Model 4: Fixed Effects with controls
model_4 <- feols(
  median_rent ~ Treated + Post + Treated_Post +
    housing_stock + white + black + native + asian +
    pacificislander + multiracial + in_poverty + housing_cost_burden_30_34 +
    housing_cost_burden_35_more + geographical_mobility + median_home_value +
    vacant_housing + median_year_built + homeownership_rate + median_income +
    unemployment_rate + high_school_grads + college_grads + commute_time +
    labor_force_participation + population | GEOID + Year,
  data = did_data_final,
  cluster = ~station_associated
)

model_list <- list(
  "DiD (No Controls)" = model_1,
  "Fixed Effects (No Controls)" = model_2,
  "DiD (With Controls)" = model_3,
  "Fixed Effects (With Controls)" = model_4
)

custom_rows <- tribble(
  ~term,                        ~"DiD (No Controls)", ~"Fixed Effects (No Controls)", ~"DiD (With Controls)", ~"Fixed Effects (With Controls)",
  "Difference-in-differences",  "Yes",               "No",                          "Yes",                  "No",
  "Treated",                    "Shown",             "Omitted (FE)",                "Shown",                "Omitted (FE)",
  "Post",                       "Shown",             "Omitted (FE)",                "Shown",                "Omitted (FE)",
  "Fixed Effects: Census Tract","No",               "Yes",                         "No",                   "Yes",
  "Fixed Effects: Year",        "No",               "Yes",                         "No",                   "Yes",
  "Controls",                   "No",               "No",                          "Yes",                  "Yes"
)

modelsummary(
  model_list,
  gof_omit = "R2|RMSE|IC|Log|Std", # Omit unnecessary goodness-of-fit stats
  add_rows = custom_rows,          # Add custom rows to the table
  coef_map = c(
    "Treated_Post" = "Treated * Post",
    "Treated" = "Treated",
    "Post" = "Post"
  ), # Consistent names for coefficients
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01) # Custom significance levels
)

---
  
# Different model specifications - housing stock
  
# Model 5: Difference-in-Differences (DiD) without controls
model_5 <- feols(
  housing_stock ~ Treated + Post + Treated_Post,
  data = did_data_final,
  cluster = ~station_associated
)

# Model 6: Fixed Effects with no controls
model_6 <- feols(
  housing_stock ~ Treated + Post + Treated_Post | GEOID + Year,
  data = did_data_final,
  cluster = ~station_associated
)

# Model 7: Difference-in-Differences (DiD) with controls
model_7 <- feols(
  housing_stock ~ Treated + Post + Treated_Post +
    white + black + native + asian +
    pacificislander + multiracial + in_poverty + housing_cost_burden_30_34 +
    housing_cost_burden_35_more + geographical_mobility + median_home_value +
    vacant_housing + median_year_built + homeownership_rate + median_income +
    unemployment_rate + high_school_grads + college_grads + commute_time +
    labor_force_participation + population,
  data = did_data_final,
  cluster = ~station_associated
)

# Model 8: Fixed Effects with controls
model_8 <- feols(
  housing_stock ~ Treated + Post + Treated_Post +
    white + black + native + asian +
    pacificislander + multiracial + in_poverty + housing_cost_burden_30_34 +
    housing_cost_burden_35_more + geographical_mobility + median_home_value +
    vacant_housing + median_year_built + homeownership_rate + median_income +
    unemployment_rate + high_school_grads + college_grads + commute_time +
    labor_force_participation + population | GEOID + Year,
  data = did_data_final,
  cluster = ~station_associated
)

model_list <- list(
  "DiD (No Controls)" = model_5,
  "Fixed Effects (No Controls)" = model_6,
  "DiD (With Controls)" = model_7,
  "Fixed Effects (With Controls)" = model_8
)

custom_rows <- tribble(
  ~term,                        ~"DiD (No Controls)", ~"Fixed Effects (No Controls)", ~"DiD (With Controls)", ~"Fixed Effects (With Controls)",
  "Difference-in-differences",  "Yes",                       "No",                                   "Yes",                           "No",
  "Treated",                    "Shown",                     "Omitted (FE)",                        "Shown",                         "Omitted (FE)",
  "Post",                       "Shown",                     "Omitted (FE)",                        "Shown",                         "Omitted (FE)",
  "Fixed Effects: Census Tract","No",                       "Yes",                                 "No",                           "Yes",
  "Fixed Effects: Year",        "No",                       "Yes",                                 "No",                           "Yes",
  "Controls",                   "No",                       "No",                                  "Yes",                          "Yes"
)

# Generate model summary with significance levels
modelsummary(
  model_list,
  gof_omit = "R2|RMSE|IC|Log|Std",
  add_rows = custom_rows,
  coef_map = c(
    "Treated_Post" = "Treated * Post",
    "Treated" = "Treated",
    "Post" = "Post"
  ), # Consistent names for coefficients
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01)
)

---
  
# New summary stats
  
variables_to_summarize <- c(
  "median_rent", "housing_stock", "white", "black", "native", "asian",
  "pacificislander", "multiracial", "in_poverty", "housing_cost_burden_30_34",
  "housing_cost_burden_35_more", "geographical_mobility", "median_home_value",
  "vacant_housing", "median_year_built", "homeownership_rate", "median_income",
  "unemployment_rate", "high_school_grads", "college_grads", "commute_time",
  "labor_force_participation", "population"
)

summary_stats <- did_data_final %>%
  filter(year_num == 6) %>%
  group_by(Treated) %>%
  summarise(across(
    all_of(variables_to_summarize),
    list(mean = ~ round(mean(.x, na.rm = TRUE), 2), 
         sd = ~ round(sd(.x, na.rm = TRUE), 2)),
    .names = "{.col}.{.fn}"
  ))

summary_stats_long <- summary_stats %>%
  pivot_longer(
    cols = -Treated,
    names_to = c("Variable", "Statistic"),
    names_sep = "\\."
  ) %>%
  pivot_wider(
    names_from = c("Statistic", "Treated"),
    values_from = value,
    names_glue = "Treated_{Treated}_{Statistic}"
  )

summary_stats_long %>%
  kable("html", caption = "Summary Statistics by Treated Status") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed", "responsive"))

---
# Buffer zone playground, sandbox, whatever you want to call it
# Basically I switch a few parameters and rerun my DiD's 

# Fit the model
lm_model <- feols(median_rent ~ Treated + Post + Treated_Post, 
                  data = did_data_final, 
                  cluster = ~station_associated)

# Customize the table
modelsummary(lm_model,
             title = "Base Difference-in-Differences Model",
             estimate = "{estimate}{stars}", # Add significance stars
             statistic = "({std.error})",    # Standard errors in parentheses below coefficients
             stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01), # Define significance levels
             gof_omit = "RMSE|Log|Adj",     # Omit RMSE and other unnecessary statistics
             coef_map = c("Treated" = "Treated", "Post" = "Post", "Treated_Post" = "Treated * Post"))






