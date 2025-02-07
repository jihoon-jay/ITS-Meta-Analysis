# Jihoon Lim
# SRMA Functions
# December 3, 2024

#### Part A #1: List containing each outcome measure for each pain type ####
pain_category <- function(list_name) {
  # Outcome vector
  outcomes <- c("persons_dispensed", "percent_persons_dispensed", # Prevalence
                "incident_persons_dispensed", "percent_incident_persons_dispensed", # Incidence
                "prescriptions", "prescriptions_per_population", # Prescriptions
                "mme", "mme_per_day", "mme_per_population", # MME
                "mean_days_supplied") # Duration
  list_outcome <- list()
  for (outcome in outcomes) {
    list_outcome[[outcome]] <- Filter(function(x) outcome %in% colnames(x), list_name)
  }
  # Include only the outcomes with at least 3 studies
  list_outcome <- Filter(function(x) length(x) >= 3, list_outcome)
  return(list_outcome)
}

#### Part A #2: Sorting list by outcome types ####
sorted_list <- function(list_filtered, outcome_measure) {
  if (outcome_measure == "Percentage") {
    list_sorted <- subset(list_filtered, names(list_filtered) %in% 
                            c("percent_persons_dispensed","percent_incident_persons_dispensed"))
  } else if (outcome_measure == "Count") {
    list_sorted <- subset(list_filtered, names(list_filtered) %in% 
                            c("prescriptions", "prescriptions_per_population", 
                              "persons_dispensed", "incident_persons_dispensed"))
  } else {
    list_sorted <- subset(list_filtered, names(list_filtered) %in% 
                            c("mme", "mme_per_day", "mme_per_population", "mean_days_supplied"))
  }
  return(list_sorted)
}

#### Part B: ITS Regression ####
# | Newey-West Outputs ####
nw_reg <- function(m_df, lag_num) {
  # m_df: data frame (m or m2)
  # lag_num: lag number
  ar1.vcov = NeweyWest(m_df, lag=lag_num)
  
  # Find std.err for parameter on 'intercept'
  se.int = sqrt(ar1.vcov[1,1])
  int.est <- round(coef(m_df)[1], 4)
  lcl.int <- round(coef(m_df)[1]-se.int*1.96, 4)
  ucl.int <- round(coef(m_df)[1]+se.int*1.96, 4)
  
  # Find std.err for parameter on 'time'
  se.time = sqrt(ar1.vcov[2,2])
  time.est <- round(coef(m_df)[2], 4)
  lcl.time <- round(coef(m_df)[2]-se.time*1.96, 4)
  ucl.time <- round(coef(m_df)[2]+se.time*1.96, 4)
  
  # Find std.err for parameter on 'level'
  se.level = sqrt(ar1.vcov[3,3])
  level.est <- round(coef(m_df)[3], 4)
  lcl.level <- round(coef(m_df)[3]-se.level*1.96, 4)
  ucl.level <- round(coef(m_df)[3]+se.level*1.96, 4)
  
  # Find std.err for parameter on 'time:level'
  time.level.est <- round(coef(m_df)[4], 4)
  se.time.level = sqrt(ar1.vcov[4,4])
  time.level.lcl <- round(coef(m_df)[4]-se.time.level*1.96, 4)
  time.level.ucl <- round(coef(m_df)[4]+se.time.level*1.96, 4)
  
  ests <- c(int.est, time.est, level.est, time.level.est)
  ses <- c(round(se.int, 4), round(se.time, 4), round(se.level, 4), round(se.time.level, 4))
  lcls <- c(lcl.int, lcl.time, lcl.level, time.level.lcl)
  ucls <- c(ucl.int, ucl.time, ucl.level, time.level.ucl)
  
  df <- as.data.frame(list(ests, ses, lcls, ucls))
  names(df) <- c("Estimate", "SE", "LCL", "UCL")
  df <- as.data.frame(t(df))
  return(df)
}

# | Log-linear regression ####
reg_log_linear <- function(list_filtered) {
  # Initialize outer list
  reg_list <- list()
  for (i in names(list_filtered)) {
    # Initialize inner list
    reg_list[[i]] <- list()
    for (j in names(list_filtered[[i]])) {
      # Use tryCatch to handle errors
      result <- tryCatch({
        # Response variable and formula setup
        response_variable <- i
        reg_formula <- as.formula(paste("log(",response_variable,") ~ time + level + time * level"))
        m <- glm(reg_formula, data = list_filtered[[i]][[j]])
        nw_reg(m_df = m, lag_num = 2) # Newey-West 
      }, error = function(e) {
        # On error, print a message and return NULL
        warning(paste("Error in i =", i, "j =", j, ":", e$message))
        NULL
      })
      
      # Store the result in the list
      reg_list[[i]][[j]] <- result
    }
  }
  reg_list_cleaned <- lapply(reg_list, function(sublist) Filter(Negate(is.null), sublist))
  return(reg_list_cleaned)
}

# | Negative binomial regression ####
reg_count <- function(list_filtered) {
  # Initialize outer list
  reg_list <- list()
  for (i in names(list_filtered)) {
    # Initialize inner list
    reg_list[[i]] <- list()
    for (j in names(list_filtered[[i]])) {
      # Use tryCatch to handle errors
      result <- tryCatch({
        # Response variable and formula setup
        response_variable <- i
        reg_formula <- as.formula(paste(response_variable, "~ time + level + time*level"))
        m <- glm.nb(reg_formula, data = list_filtered[[i]][[j]])
        nw_reg(m_df = m, lag_num = 2) # Newey-West 
      }, error = function(e) {
        # On error, print a message and return NULL
        warning(paste("Error in i =", i, "j =", j, ":", e$message))
        NULL
      })
      
      # Store the result in the list
      reg_list[[i]][[j]] <- result
    }
  }
  reg_list_cleaned <- lapply(reg_list, function(sublist) Filter(Negate(is.null), sublist))
  return(reg_list_cleaned)
}

# | Beta regression ####
reg_proportion <- function(list_filtered) {
  # Initialize outer list
  reg_list <- list()
  for (i in names(list_filtered)) {
    # Initialize inner list
    reg_list[[i]] <- list()
    for (j in names(list_filtered[[i]])) {
      # Use tryCatch to handle errors
      result <- tryCatch({
        # Response variable and formula setup
        response_variable <- i
        reg_formula <- as.formula(paste(response_variable, "~ time + level + time*level"))
        m <- betareg(reg_formula, data = list_filtered[[i]][[j]])
        nw_reg(m_df = m, lag_num = 2) # Newey-West 
      }, error = function(e) {
        # On error, print a message and return NULL
        warning(paste("Error in i =", i, "j =", j, ":", e$message))
        NULL
      })
      
      # Store the result in the list
      reg_list[[i]][[j]] <- result
    }
  }
  reg_list_cleaned <- lapply(reg_list, function(sublist) Filter(Negate(is.null), sublist))
  return(reg_list_cleaned)
}

#### Part C: Meta-Analysis ####
# | Extract ITS results ####
its_extraction <- function(list_name, col_num, ref_list) {
  # ITS regression output
  m_output <- matrix(NA, nrow = length(list_name), ncol = 4) # row = num_row for list_name[[i]], col = 4
  for (i in seq_along(list_name)) {
    m_output[i,] <- t(list_name[[i]][, col_num])
  }
  m_output <- cbind(m_output, names(list_name))
  colnames(m_output) <- c("Estimate", "SE", "LCL", "UCL", "Study")
  m_output <- as.data.frame(m_output)
  
  # Guideline names
  m_guideline <- matrix(NA, nrow = length(ref_list), ncol = 2)
  m_guideline[, 1] <- names(ref_list)
  for (i in seq_along(ref_list)) {
    m_guideline[i, 2] <- ref_list[[i]]$guideline[1]
  }
  colnames(m_guideline) <- c("Study", "Guideline")
  m_guideline <- as.data.frame(m_guideline)
  
  # Left join
  m_example <- left_join(m_output, m_guideline, by = "Study")
  return(m_example)
}

# | Percent Changes ####
rr_to_percent_change <- function(rr) {
  percent_change <- (exp(rr) - 1) * 100
  return(percent_change)
}

# | Meta-Analysis ####
re_meta <- function(reg_dataset, coef_text1, coef_text2) {
  # Initial meta-analysis of regression coefficients
  ma_time_level <- metagen(TE = as.numeric(Estimate), seTE = as.numeric(SE), data = reg_dataset,
                           sm = "RR", common = FALSE, random = TRUE, method.tau = "REML",
                           studlab = reg_dataset$Study)
  
  # Transform RR and 95% CI to percent change
  ## TE: study-specific
  ## TE.random: random-effects TE from meta-MA
  ma_time_level$TE <- rr_to_percent_change(ma_time_level$TE)
  ma_time_level$TE.random <- rr_to_percent_change(ma_time_level$TE.random)
  tes <- c(ma_time_level$TE, ma_time_level$TE.random)
  ## lower: study-specific
  ## lower.random: random-effects LCL from MA
  ma_time_level$lower <- rr_to_percent_change(ma_time_level$lower)
  ma_time_level$lower.random <- rr_to_percent_change(ma_time_level$lower.random)
  lcls <- c(ma_time_level$lower, ma_time_level$lower.random)
  ## upper: study-specific
  ## upper.random: random-effects UCL from MA
  ma_time_level$upper <- rr_to_percent_change(ma_time_level$upper)
  ma_time_level$upper.random <- rr_to_percent_change(ma_time_level$upper.random)
  ucls <-c(ma_time_level$upper, ma_time_level$upper.random)
  
  # Study Characteristics
  studies <- c(ma_time_level$studlab, "Pooled Estimate")
  weights <- c(round(ma_time_level$w.random/sum(ma_time_level$w.random) * 100, 2), NA)
  pooled = c(rep(FALSE, length(reg_dataset[,1])), TRUE)
  guidelines <- c(reg_dataset$Guideline, NA)
  
  # Combine data into a data frame for forest plot
  forest_data <- data.frame(Study = studies, Effect = round(tes, 1),
                            LCI = round(lcls, 1), UCI = round(ucls, 1),
                            Weight = weights, Guideline = guidelines)
  tau_sq <- ma_time_level$tau2; print(paste("Tau^2: ", tau_sq)) # Tau^2
  i_sq <- ma_time_level$I2; print(paste("I^2: ", i_sq)) # I^2
  
  # Forest Plot
  forest <- forestplot(forest_data, labeltext = c(Study, Effect, LCI, UCI, Weight, Guideline),
                       mean = Effect, lower = LCI, upper = UCI,
                       is.summary = c(rep(FALSE, length(reg_dataset[, 1])), TRUE),
                       col = fpColors(box = "black", lines = "black", summary = "darkblue"),
                       colgap = unit(2, "mm"),
                       title = paste("Pooled", coef_text1, "Change: ", coef_text2),
                       txt_gp = fpTxtGp(label = gpar(fontsize = 8), ticks = gpar(fontsize = 16), xlab  = gpar(fontsize = 16)),
                       xticks = seq(-2.5, 2.5, 0.5), boxsize = 0.15,
                       xlab = "% Change") %>%
    fp_add_header(Study = "Study", Effect = "% Change", LCI = "LCL", UCI = "UCL", 
                  Weight = "Weight (%)", Guideline = "Guideline") %>%
    fp_add_lines() %>% 
    fp_set_style(box = "royalblue", line = "royalblue4", align = "lccccc")
  grid.text(paste("I² =", round(i_sq*100, 1), "%"), x = 0.9, y = 0.05, gp = gpar(fontsize = 9))
  print(forest)
}
