# Jihoon Lim
# SRMA
# Load packages
library(readxl); library(dplyr); library(tidyr); library(openxlsx); library(lubridate)
library(sandwich); library(betareg); library(MASS); library(meta); library(forestplot)
library(grid)

#### Directory Setup ####
folder_path <- "C:/Users/limji/Desktop/Research Associate/SRMA/Tables"
file_list <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
data_list <- lapply(file_list, read_excel)
names(data_list) <- tools::file_path_sans_ext(basename(file_list))

# | Exclude observations with < 8 time points ####
list_filtered <- Filter(function(x) {is.data.frame(x) && nrow(x) >= 8}, data_list)
list_gen <- Filter(function(x) "pain_type" %in% colnames(x) 
                   && any(x$pain_type == "General"), list_filtered)
list_cnc <- Filter(function(x) "pain_type" %in% colnames(x) 
                   && any(x$pain_type == "CNCP"), list_filtered)
list_chr <- Filter(function(x) "pain_type" %in% colnames(x) 
                   && any(x$pain_type == "Chronic"), list_filtered)
list_acu <- Filter(function(x) "pain_type" %in% colnames(x) 
                   && any(x$pain_type %in% c("Acute", "Surgical", "ED", "Trauma")), list_filtered)
list_oth <- Filter(function(x) "pain_type" %in% colnames(x) 
                   && any(x$pain_type == "Other"), list_filtered)

#### Part A: Create Lists ####
# | 1. Create a list containing each outcome measure for each pain type ####
a1_gen <- pain_category(list_gen) # General
a1_cnc <- pain_category(list_cnc) # CNCP
a1_chr <- pain_category(list_chr) # Chronic
a1_acu <- pain_category(list_acu) # Acute

# | 2. Sort lists by outcome types ####
# General
a2_gen_prop <- sorted_list(a1_gen, "Percentage") # N = 2
a2_gen_count <- sorted_list(a1_gen, "Count") # N = 3
a2_gen_cont <- sorted_list(a1_gen, "Continuous") # N = 3
# CNCP
a2_cnc_prop <- sorted_list(a1_cnc, "Percentage") # N = 1, n = 11
a2_cnc_count <- sorted_list(a1_cnc, "Count") # N = 0
a2_cnc_cont <- sorted_list(a1_cnc, "Continuous") # N = 2, n = 5, 3
# Chronic
a2_chr_prop <- sorted_list(a1_chr, "Percentage") # N = 2, n = 14, 4
a2_chr_count <- sorted_list(a1_chr, "Count") # N = 0
a2_chr_cont <- sorted_list(a1_chr, "Continuous") # N = 0
# Acute
a2_acu_prop <- sorted_list(a1_acu, "Percentage") # N = 0
a2_acu_count <- sorted_list(a1_acu, "Count") # N = 0
a2_acu_cont <- sorted_list(a1_acu, "Continuous") # N = 0

#### Part B: ITS Regression ####
# | 1. Run regression for all outcomes ####
# Note: This removes studies that had too few time points.
# General
b1_gen_cont <- reg_log_linear(a2_gen_cont)
b1_gen_count <- reg_count(a2_gen_count)
b1_gen_prop <- reg_proportion(a2_gen_prop)
# CNCP
b1_cnc_cont <- reg_log_linear(a2_cnc_cont)
b1_cnc_prop <- reg_proportion(a2_cnc_prop)
# Chronic
b1_chr_prop <- reg_proportion(a2_chr_prop)

#### Part C: Meta-Analysis ####
# | 1. Extract ITS results ####
# || a. Percent persons dispensed (%) ####
# General (N = 14)
c1_gen_prev_level <- its_extraction(b1_gen_prop[[1]], 3, a2_gen_prop[[1]])
c1_gen_prev_trend <- its_extraction(b1_gen_prop[[1]], 4, a2_gen_prop[[1]])
# CNCP (N = 10)
c1_cnc_prev_level <- its_extraction(b1_cnc_prop[[1]], 3, a2_cnc_prop[[1]])
c1_cnc_prev_trend <- its_extraction(b1_cnc_prop[[1]], 4, a2_cnc_prop[[1]])
# Chronic (N = 13)
c1_chr_prev_level <- its_extraction(b1_chr_prop[[1]], 3, a2_chr_prop[[1]])
c1_chr_prev_trend <- its_extraction(b1_chr_prop[[1]], 4, a2_chr_prop[[1]])

# || b. Percent incident persons dispensed ####
# General (N = 3)
c1_gen_inc_level <- its_extraction(b1_gen_prop[[2]], 3, a2_gen_prop[[2]])
c1_gen_inc_trend <- its_extraction(b1_gen_prop[[2]], 4, a2_gen_prop[[2]])

# || c. MME per day ####
# General (N = 6)
c1_gen_mme_day_level <- its_extraction(b1_gen_cont[[1]], 3, a2_gen_cont[[1]])
c1_gen_mme_day_trend <- its_extraction(b1_gen_cont[[1]], 4, a2_gen_cont[[1]])

# || d. MME per person ####
# General (N = 7)
c1_gen_mme_pop_level <- its_extraction(b1_gen_cont[[2]], 3, a2_gen_cont[[2]])
c1_gen_mme_pop_trend <- its_extraction(b1_gen_cont[[2]], 4, a2_gen_cont[[2]])
# CNCP (N = 3, 5)
c1_cnc_mme_pop_level <- its_extraction(b1_cnc_cont[[1]], 3, a2_cnc_cont[[1]])
c1_cnc_mme_pop_trend <- its_extraction(b1_cnc_cont[[1]], 4, a2_cnc_cont[[1]])

# || e. Persons dispensed ####
# General (N = 4)
c1_gen_prevn_level <- its_extraction(b1_gen_count[[1]], 3, a2_gen_count[[1]])
c1_gen_prevn_trend <- its_extraction(b1_gen_count[[1]], 4, a2_gen_count[[1]])

# || f. Prescriptions ####
# General (N = 12)
c1_gen_rx_level <- its_extraction(b1_gen_count[[3]], 3, a2_gen_count[[3]])
c1_gen_rx_trend <- its_extraction(b1_gen_count[[3]], 4, a2_gen_count[[3]])

# || g. Duration ####
# General (N = 4)
c1_gen_dur_level <- its_extraction(b1_gen_cont[[3]], 3, a2_gen_cont[[3]])
c1_gen_dur_trend <- its_extraction(b1_gen_cont[[3]], 4, a2_gen_cont[[3]])
# CNCP (N = 3)
c1_cnc_dur_level <- its_extraction(b1_cnc_cont[[2]], 3, a2_cnc_cont[[2]])
c1_cnc_dur_trend <- its_extraction(b1_cnc_cont[[2]], 4, a2_cnc_cont[[2]])

# | 2. Meta-analysis ####
# || a. Percent persons dispensed (%) ####
# General (lower level, lower trend)
re_meta(c1_gen_prev_level, "Level", "Prevalence")
re_meta(c1_gen_prev_trend, "Trend", "Prevalence")
# CNCP (no changes in level or trend)
re_meta(c1_cnc_prev_level, "Level", "Prevalence")
re_meta(c1_cnc_prev_trend, "Trend", "Prevalence")
# Chronic (non-significant lower level, lower trend)
re_meta(c1_chr_prev_level, "Level", "Prevalence")
re_meta(c1_chr_prev_trend, "Trend", "Prevalence")

# || b. Percent incident persons dispensed ####
# General (no changes in level or trend)
re_meta(c1_gen_inc_level, "Level", "Incidence")
re_meta(c1_gen_inc_trend, "Trend", "Incidence")

# || c. MME per day ####
# General (no level change, declining trend)
re_meta(c1_gen_mme_day_level, "Level", "MME per Day")
re_meta(c1_gen_mme_day_trend, "Trend", "MME per Day")

# || d. MME per person ####
# General (no level change, declining trend)
re_meta(c1_gen_mme_pop_level, "Level", "MME per Person")
re_meta(c1_gen_mme_pop_trend, "Trend", "MME per Person")
# CNCP (no level change, declining trend)
re_meta(c1_cnc_mme_pop_level, "Level", "MME per Person")
re_meta(c1_cnc_mme_pop_trend, "Trend", "MME per Person")

# || e. Persons dispensed ####
# General (no level change, declining trend)
re_meta(c1_gen_prevn_level, "Level", "Persons Dispensed")
re_meta(c1_gen_prevn_trend, "Trend", "Persons Dispensed")

# || f. Prescriptions ####
# General (marginally insignificant level change, declining trend)
re_meta(c1_gen_rx_level, "Level", "Number of Prescriptions")
re_meta(c1_gen_rx_trend, "Trend", "Number of Prescriptions")

# || g. Duration ####
# General (no changes in level or trend)
re_meta(c1_gen_dur_level, "Level", "Prescription Duration")
re_meta(c1_gen_dur_trend, "Trend", "Prescription Duration")
# CNCP (no changes in level, declining trend)
re_meta(c1_cnc_dur_level, "Level", "Prescription Duration")
re_meta(c1_cnc_dur_trend, "Trend", "Prescription Duration")
