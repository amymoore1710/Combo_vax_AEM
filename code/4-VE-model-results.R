#Adding it to the simulation

pacman::p_load(dplyr,
               tidyr,
               future,
               furrr,
               purrr,
               scales,
               stringr,
               ggplot2,
               survival,
               here)

today_date <- "2026-09-15"


dataset_names <- c("0.1", "0.2", "0.3", "0.4", 
                   "0.5", "0.6", "0.7", "0.8", "0.9")

file_names <- paste0("ve_0_", dataset_names, "_0_0.12")


dataset_paths <- paste0("sim data/sim_results_", 
                        file_names,
                        ".RDS")

any_shig_estimates <- data.frame(scenario = NA,
                                 country = NA,
                                 true_ve = NA,
                                 estimated_ve = NA,
                                 SE_ve = NA,
                                 lowCI_ve = NA,
                                 highCI_ve = NA,
                                 CI_contains = NA,
                                 num_events = NA,
                                 total_follow_up = NA,
                                 num_coinf = NA)

shig_CU_estimates <- data.frame(scenario = NA,
                                  country = NA,
                                 true_ve = NA,
                                 estimated_ve = NA,
                                 SE_ve = NA,
                                 lowCI_ve = NA,
                                 highCI_ve = NA,
                                 CI_contains = NA,
                                 num_events = NA,
                                 total_follow_up = NA,
                                 num_coinf = NA)

weighted_estimates <- data.frame(scenario = NA,
                                 country = NA,
                                 true_ve = NA,
                                 estimated_ve = NA,
                                 SE_ve = NA,
                                 lowCI_ve = NA,
                                 highCI_ve = NA,
                                 CI_contains = NA,
                                 num_events = NA,
                                 total_follow_up = NA,
                                 num_coinf = NA)

# weighted_adj_estimates <- data.frame(scenario = NA,
#                                  country = NA,
#                                  true_ve = NA,
#                                  estimated_ve = NA,
#                                  SE_ve = NA,
#                                  lowCI_ve = NA,
#                                  highCI_ve = NA,
#                                  CI_contains = NA,
#                                  num_events = NA,
#                                  total_follow_up = NA,
#                                  num_coinf = NA)
# 
# weighted_alt_adj_estimates <- data.frame(scenario = NA,
#                                      country = NA,
#                                      true_ve = NA,
#                                      estimated_ve = NA,
#                                      SE_ve = NA,
#                                      lowCI_ve = NA,
#                                      highCI_ve = NA,
#                                      CI_contains = NA,
#                                      num_events = NA,
#                                      total_follow_up = NA,
#                                      num_coinf = NA)

for (i in 1:length(file_names)) {
  
  scenario <- dataset_names[i]
  dataset_path <- here::here(dataset_paths[i])
  file_name <- file_names[i]
  
  # true VE against disease
  super_true_ve <- sub("^[^_]*_[^_]*_([^_]*).*", "\\1", file_name)
  
  # Load dataset
  raw_data <- readRDS(dataset_path)
  
  n_simulations <- length(unique(raw_data$sim))
  
  # Separate by simulation
  for (j in 1:n_simulations) {
    data_by_country <- raw_data %>% filter(sim == j)
    
    countries <- unique(data_by_country$country_id)
    # for (k in 1:length(countries)) {
      # country <- countries[k]
      country <- "PK"
      
      data <- data_by_country %>% filter(country_id == country)
    
      # Defining Events
        # Any Shigella: event = diarrhea with shigella AFE > 0
        # Highest AFE: event = diarrhea with shigella AFE > all other AFEs
        # High Shigella: event = diarrhea with shigella AFE > 0.5
        # Weighted model: event = any shigella event, weight = shigella AFE
      
      data <- data %>%
        mutate(any_shig_event = ifelse(diarrhea == 1 & shigella_afe > 0, 
                                       1, 0),
               event_weight = ifelse(any_shig_event == 1, shigella_afe, 1),
               # event_adj_weight = ifelse(any_shig_event == 1, shigella_adj_afe, 1) ,
               # event_alt_adj_weight = ifelse(any_shig_event == 1, shigella_alt_adj_afe, 1),
               shig_CU_event = ifelse(diarrhea == 1 & shigella_ct < 28.8,
                                        1, 0)
               )
      
      vax_positives <- sum(data %>% filter(vax == 1) %>% select(shigella_diarrhea))
      vax_negatives <- nrow(data %>% filter(vax == 1) %>% select(shigella_diarrhea)) - vax_positives
      unvax_positives <- sum(data %>% filter(vax == 0) %>% select(shigella_diarrhea))
      unvax_negatives <- nrow(data %>% filter(vax == 0) %>% select(shigella_diarrhea)) - unvax_positives
      
      true_ve <- 1 - (vax_positives/vax_negatives)/(unvax_positives/unvax_negatives)
      
      # checking <- data %>% summarize(weights = sum(event_weight))
      
      
      data <- data %>%
        mutate(unique_id = paste0(sim,"_",country_id, "_",sprintf("%04d", child)))
      
      
      #Any Shigella Model
      
      any_shig_data <- data %>%
        group_by(unique_id) %>%
        arrange(month) %>%
        # Keep rows up to and including the first time 'any_shig_event' is 1
        filter(cumsum(lag(any_shig_event, default = 0)) < 1) %>%
        arrange(unique_id) %>% 
        mutate(prev_month = month - 1)
      
      
      any_shig_model <- coxph(Surv(time = prev_month,
                                   time2 = month,
                                   event = any_shig_event) ~ vax,
                              data = any_shig_data)
      summary(any_shig_model)
      
      logHR <- summary(any_shig_model)$coefficients[1]
      SE_ve <- summary(any_shig_model)$coefficients[3]
      
      estimated_ve <- 1 - exp(logHR)
      lowCI_logHR <- logHR + qnorm(p = 0.025) * SE_ve 
      highCI_logHR <- logHR + qnorm(p = 0.975) * SE_ve 
      
      lowCI_ve <- 1 - exp(highCI_logHR)
      highCI_ve <- 1 - exp(lowCI_logHR)
      
      CI_contains <- ifelse(lowCI_ve <= true_ve & highCI_ve >= true_ve, 1, 0)
      
      # print(paste0("True VE: ", super_true_ve))
      # print(paste0("Estimated VE: ", estimated_ve))
      
      
        #Number of events 
      num_events <- sum(any_shig_data$any_shig_event)
      
        #Total Follow up time --> the total number of person months before censoring
      total_follow_up <- nrow(any_shig_data)
      
        #Coinfections are events that have a nonzero AFe for another pathogen
      num_coinf <- sum(ifelse(any_shig_data$any_shig_event == 1 & (any_shig_data$ETEC_afe > 0 | any_shig_data$other_afe > 0),
                              1, 0))
      
      any_shig_estimates <- rbind(any_shig_estimates, 
                                  c(i,
                                    country,
                                    true_ve,
                                    estimated_ve,
                                    SE_ve,
                                    lowCI_ve,
                                    highCI_ve,
                                    CI_contains,
                                    num_events,
                                    total_follow_up,
                                    num_coinf))
      
      #Shigella Cut-Off Model
      
      shig_CU_data <- data %>%
        group_by(unique_id) %>%
        arrange(month) %>%
        # Keep rows up to and including the first time 'any_shig_event' is 1
        filter(cumsum(lag(shig_CU_event, default = 0)) < 1) %>%
        arrange(unique_id) %>% 
        mutate(prev_month = month - 1)
      
      
      shig_CU_model <- coxph(Surv(time = prev_month,
                                   time2 = month,
                                   event = shig_CU_event) ~ vax,
                              data = shig_CU_data)
      summary(shig_CU_model)
      
      logHR <- summary(shig_CU_model)$coefficients[1]
      SE_ve <- summary(shig_CU_model)$coefficients[3]
      
      estimated_ve <- 1 - exp(logHR)
      lowCI_logHR <- logHR + qnorm(p = 0.025) * SE_ve 
      highCI_logHR <- logHR + qnorm(p = 0.975) * SE_ve 
      
      lowCI_ve <- 1 - exp(highCI_logHR)
      highCI_ve <- 1 - exp(lowCI_logHR)
      
      CI_contains <- ifelse(lowCI_ve <= true_ve & highCI_ve >= true_ve, 1, 0)
      
      # print(paste0("True VE: ", super_true_ve))
      # print(paste0("Estimated VE: ", estimated_ve))
      
      
      #Number of events 
      num_events <- sum(shig_CU_data$shig_CU_event)
      
      #Total Follow up time --> the total number of person months before censoring
      total_follow_up <- nrow(shig_CU_data)
      
      #Coinfections are events that have a nonzero AFe for another pathogen
      num_coinf <- sum(ifelse(shig_CU_data$shig_CU_event == 1 & (shig_CU_data$ETEC_afe > 0 | shig_CU_data$other_afe > 0),
                              1, 0))
      
      shig_CU_estimates <- rbind(shig_CU_estimates, 
                                  c(i,
                                    country,
                                    true_ve,
                                    estimated_ve,
                                    SE_ve,
                                    lowCI_ve,
                                    highCI_ve,
                                    CI_contains,
                                    num_events,
                                    total_follow_up,
                                    num_coinf))
  
      
      
      
      
      #Weighted Model
      
      weighted_data <- data %>%
        group_by(unique_id) %>%
        arrange(month) %>%
        # Keep rows up to and including the first time 'any_shig_event' is 1
        filter(cumsum(lag(any_shig_event, default = 0)) < 1) %>%
        arrange(unique_id) %>% 
        mutate(prev_month = month - 1)
      
      extra_events <- weighted_data %>%
        filter(any_shig_event == 1) %>%
        mutate(any_shig_event = 0,
               event_weight = 1 - event_weight)
      
      weighted_data <- rbind(weighted_data, extra_events) %>%
        arrange(unique_id, month)
      
      
      weighted_model <- coxph(Surv(time = prev_month,
                                   time2 = month,
                                   event = any_shig_event) ~ vax,
                              weights = event_weight,
                              data = weighted_data)
      # summary(weighted_model)
      
      logHR <- summary(weighted_model)$coefficients[1]
      SE_ve <- summary(weighted_model)$coefficients[3]
      
      estimated_ve <- 1 - exp(logHR)
      lowCI_logHR <- logHR + qnorm(p = 0.025) * SE_ve 
      highCI_logHR <- logHR + qnorm(p = 0.975) * SE_ve 
      
      lowCI_ve <- 1 - exp(highCI_logHR)
      highCI_ve <- 1 - exp(lowCI_logHR)
      
      CI_contains <- ifelse(lowCI_ve <= true_ve & highCI_ve >= true_ve, 1, 0)
      
      # print(paste0("True VE: ", super_true_ve))
      # print(paste0("Estimated VE: ", estimated_ve))
      
        #Number of events 
      num_events <- sum(weighted_data$any_shig_event)
      
        #Total Follow up time --> the total number of person months before censoring
      total_follow_up <- nrow(weighted_data)
      
        #Coinfections are events that have a nonzero AFe for another pathogen
      num_coinf <- sum(ifelse(weighted_data$any_shig_event == 1 & (weighted_data$ETEC_afe > 0 | weighted_data$other_afe > 0),
                              1, 0))
      
      
      weighted_estimates <- rbind(weighted_estimates, 
                                  c(i,
                                    country,
                                    true_ve,
                                    estimated_ve,
                                    SE_ve,
                                    lowCI_ve,
                                    highCI_ve,
                                    CI_contains,
                                    num_events,
                                    total_follow_up,
                                    num_coinf))
      
      
      # #Weighted Model - Sum to 1 adjustment
      # 
      # weighted_adj_data <- data %>%
      #   group_by(unique_id) %>%
      #   arrange(month) %>%
      #   # Keep rows up to and including the first time 'any_shig_event' is 1
      #   filter(cumsum(lag(any_shig_event, default = 0)) < 1) %>%
      #   arrange(unique_id) %>% 
      #   mutate(prev_month = month - 1)
      # 
      # 
      # weighted_adj_model <- coxph(Surv(time = prev_month,
      #                              time2 = month,
      #                              event = any_shig_event) ~ vax,
      #                         weights = event_adj_weight,
      #                         data = weighted_adj_data)
      # # summary(weighted_model)
      # 
      # logHR <- summary(weighted_adj_model)$coefficients[1]
      # SE_ve <- summary(weighted_adj_model)$coefficients[3]
      # 
      # estimated_ve <- 1 - exp(logHR)
      # lowCI_logHR <- logHR + qnorm(p = 0.025) * SE_ve 
      # highCI_logHR <- logHR + qnorm(p = 0.975) * SE_ve 
      # 
      # lowCI_ve <- 1 - exp(highCI_logHR)
      # highCI_ve <- 1 - exp(lowCI_logHR)
      # 
      # CI_contains <- ifelse(lowCI_ve <= true_ve & highCI_ve >= true_ve, 1, 0)
      # 
      # # print(paste0("True VE: ", super_true_ve))
      # # print(paste0("Estimated VE: ", estimated_ve))
      # 
      # #Number of events 
      # num_events <- sum(weighted_adj_data$any_shig_event)
      # 
      # #Total Follow up time --> the total number of person months before censoring
      # total_follow_up <- nrow(weighted_adj_data)
      # 
      # #Coinfections are events that have a nonzero AFe for another pathogen
      # num_coinf <- sum(ifelse(weighted_adj_data$any_shig_event == 1 & (weighted_adj_data$ETEC_afe > 0 | weighted_adj_data$other_afe > 0),
      #                         1, 0))
      # 
      # 
      # weighted_adj_estimates <- rbind(weighted_adj_estimates, 
      #                             c(i,
      #                               country,
      #                               true_ve,
      #                               estimated_ve,
      #                               SE_ve,
      #                               lowCI_ve,
      #                               highCI_ve,
      #                               CI_contains,
      #                               num_events,
      #                               total_follow_up,
      #                               num_coinf))
      # 
      # 
      # 
      # #Weighted Model - Multiplicative Adjustment (alternate)
      # 
      # weighted_alt_adj_data <- data %>%
      #   group_by(unique_id) %>%
      #   arrange(month) %>%
      #   # Keep rows up to and including the first time 'any_shig_event' is 1
      #   filter(cumsum(lag(any_shig_event, default = 0)) < 1) %>%
      #   arrange(unique_id) %>% 
      #   mutate(prev_month = month - 1)
      # 
      # 
      # weighted_alt_adj_model <- coxph(Surv(time = prev_month,
      #                                  time2 = month,
      #                                  event = any_shig_event) ~ vax,
      #                             weights = event_alt_adj_weight,
      #                             data = weighted_alt_adj_data)
      # # summary(weighted_model)
      # 
      # logHR <- summary(weighted_alt_adj_model)$coefficients[1]
      # SE_ve <- summary(weighted_alt_adj_model)$coefficients[3]
      # 
      # estimated_ve <- 1 - exp(logHR)
      # lowCI_logHR <- logHR + qnorm(p = 0.025) * SE_ve 
      # highCI_logHR <- logHR + qnorm(p = 0.975) * SE_ve 
      # 
      # lowCI_ve <- 1 - exp(highCI_logHR)
      # highCI_ve <- 1 - exp(lowCI_logHR)
      # 
      # CI_contains <- ifelse(lowCI_ve <= true_ve & highCI_ve >= true_ve, 1, 0)
      # 
      # # print(paste0("True VE: ", super_true_ve))
      # # print(paste0("Estimated VE: ", estimated_ve))
      # 
      # #Number of events 
      # num_events <- sum(weighted_alt_adj_data$any_shig_event)
      # 
      # #Total Follow up time --> the total number of person months before censoring
      # total_follow_up <- nrow(weighted_alt_adj_data)
      # 
      # #Coinfections are events that have a nonzero AFe for another pathogen
      # num_coinf <- sum(ifelse(weighted_alt_adj_data$any_shig_event == 1 & (weighted_alt_adj_data$ETEC_afe > 0 | weighted_alt_adj_data$other_afe > 0),
      #                         1, 0))
      # 
      # 
      # weighted_alt_adj_estimates <- rbind(weighted_adj_estimates, 
      #                                 c(i,
      #                                   country,
      #                                   true_ve,
      #                                   estimated_ve,
      #                                   SE_ve,
      #                                   lowCI_ve,
      #                                   highCI_ve,
      #                                   CI_contains,
      #                                   num_events,
      #                                   total_follow_up,
      #                                   num_coinf))
    
    # }
  }
}

any_shig_estimates <- any_shig_estimates %>% filter(!is.na(true_ve))

any_shig_estimates <- any_shig_estimates %>%
  mutate(CI_contains = as.numeric(CI_contains),
         true_ve = as.numeric(true_ve),
         estimated_ve = as.numeric(estimated_ve),
         SE_ve = as.numeric(SE_ve),
         method = "Any Detection")

shig_CU_estimates <- shig_CU_estimates %>% filter(!is.na(true_ve))

shig_CU_estimates <- shig_CU_estimates %>%
  mutate(CI_contains = as.numeric(CI_contains),
         true_ve = as.numeric(true_ve),
         estimated_ve = as.numeric(estimated_ve),
         SE_ve = as.numeric(SE_ve),
         method = "Detection Cut Off")

weighted_estimates <- weighted_estimates %>% filter(!is.na(true_ve))

weighted_estimates <- weighted_estimates %>%
  mutate(CI_contains = as.numeric(CI_contains),
         true_ve = as.numeric(true_ve),
         estimated_ve = as.numeric(estimated_ve),
         SE_ve = as.numeric(SE_ve),
         method = "Weighted")

# weighted_adj_estimates <- weighted_adj_estimates %>% filter(!is.na(true_ve))
# 
# weighted_adj_estimates <- weighted_adj_estimates %>%
#   mutate(CI_contains = as.numeric(CI_contains),
#          true_ve = as.numeric(true_ve),
#          estimated_ve = as.numeric(estimated_ve),
#          SE_ve = as.numeric(SE_ve),
#          method = "Weighted Adjusted")
# 
# weighted_alt_adj_estimates <- weighted_alt_adj_estimates %>% filter(!is.na(true_ve))
# 
# weighted_alt_adj_estimates <- weighted_alt_adj_estimates %>%
#   mutate(CI_contains = as.numeric(CI_contains),
#          true_ve = as.numeric(true_ve),
#          estimated_ve = as.numeric(estimated_ve),
#          SE_ve = as.numeric(SE_ve),
#          method = "Weighted Alternate Adjusted")



estimates <- rbind(any_shig_estimates,
                   shig_CU_estimates,
                   weighted_estimates) #,
                   # weighted_adj_estimates,
                   # weighted_alt_adj_estimates)
estimates <- estimates %>%
  mutate(VE = as.numeric(scenario) *0.1)

# estimates_BG <- estimates %>% filter(country == "BG") 
# 
# true_ve_BG <- estimates_BG %>%
#   group_by(scenario) %>%
#   summarize(true_ve = mean(true_ve)) %>% 
#   ungroup()
# 
# true_ve_plot <- true_ve_BG$true_ve
# 
# 
# boxplot_BG <- estimates_BG %>% 
#   ggplot(aes(x=as.factor(VE), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("True Simulated VE") +
#   ylab("Estimated VE") +
#   geom_segment(aes(x = 0.5, y = true_ve_plot[1], xend = 1.5, yend = true_ve_plot[1]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)  +
#   geom_segment(aes(x = 1.5, y = true_ve_plot[2], xend = 2.5, yend = true_ve_plot[2]), 
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) + 
#   geom_segment(aes(x = 2.5, y = true_ve_plot[3], xend = 3.5, yend = true_ve_plot[3]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 3.5, y = true_ve_plot[4], xend = 4.5, yend = true_ve_plot[4]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 4.5, y = true_ve_plot[5], xend = 5.5, yend = true_ve_plot[5]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 5.5, y = true_ve_plot[6], xend = 6.5, yend = true_ve_plot[6]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 6.5, y = true_ve_plot[7], xend = 7.5, yend = true_ve_plot[7]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 7.5, y = true_ve_plot[8], xend = 8.5, yend = true_ve_plot[8]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 8.5, y = true_ve_plot[9], xend = 9.5, yend = true_ve_plot[9]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)
# 
# boxplot_BG
# 
# ggsave(here("figures","aem", today_date,"results_boxplot_VE_scale_BG.png"), plot = boxplot_BG, width = 6, height = 4, dpi = 300)
# 
# 
# estimates_IN <- estimates %>% filter(country == "IN") 
# 
# true_ve_IN <- estimates_IN %>%
#   group_by(scenario) %>%
#   summarize(true_ve = mean(true_ve)) %>% 
#   ungroup()
# 
# true_ve_plot <- true_ve_IN$true_ve
# 
# 
# boxplot_IN <- estimates_IN %>% 
#   ggplot(aes(x=as.factor(VE), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("True Simulated VE") +
#   ylab("Estimated VE") +
#   geom_segment(aes(x = 0.5, y = true_ve_plot[1], xend = 1.5, yend = true_ve_plot[1]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)  +
#   geom_segment(aes(x = 1.5, y = true_ve_plot[2], xend = 2.5, yend = true_ve_plot[2]), 
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) + 
#   geom_segment(aes(x = 2.5, y = true_ve_plot[3], xend = 3.5, yend = true_ve_plot[3]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 3.5, y = true_ve_plot[4], xend = 4.5, yend = true_ve_plot[4]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 4.5, y = true_ve_plot[5], xend = 5.5, yend = true_ve_plot[5]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 5.5, y = true_ve_plot[6], xend = 6.5, yend = true_ve_plot[6]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 6.5, y = true_ve_plot[7], xend = 7.5, yend = true_ve_plot[7]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 7.5, y = true_ve_plot[8], xend = 8.5, yend = true_ve_plot[8]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 8.5, y = true_ve_plot[9], xend = 9.5, yend = true_ve_plot[9]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)
# 
# boxplot_IN
# 
# ggsave(here("figures","aem", today_date,"results_boxplot_VE_scale_IN.png"), plot = boxplot_IN, width = 6, height = 4, dpi = 300)
# 
# 
# 
# estimates_PE <- estimates %>% filter(country == "PE") 
# 
# true_ve_PE <- estimates_PE %>%
#   group_by(scenario) %>%
#   summarize(true_ve = mean(true_ve)) %>% 
#   ungroup()
# 
# true_ve_plot <- true_ve_PE$true_ve
# 
# 
# boxplot_PE <- estimates_PE %>% 
#   ggplot(aes(x=as.factor(VE), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("True Simulated VE") +
#   ylab("Estimated VE") +
#   geom_segment(aes(x = 0.5, y = true_ve_plot[1], xend = 1.5, yend = true_ve_plot[1]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)  +
#   geom_segment(aes(x = 1.5, y = true_ve_plot[2], xend = 2.5, yend = true_ve_plot[2]), 
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) + 
#   geom_segment(aes(x = 2.5, y = true_ve_plot[3], xend = 3.5, yend = true_ve_plot[3]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 3.5, y = true_ve_plot[4], xend = 4.5, yend = true_ve_plot[4]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 4.5, y = true_ve_plot[5], xend = 5.5, yend = true_ve_plot[5]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 5.5, y = true_ve_plot[6], xend = 6.5, yend = true_ve_plot[6]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 6.5, y = true_ve_plot[7], xend = 7.5, yend = true_ve_plot[7]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 7.5, y = true_ve_plot[8], xend = 8.5, yend = true_ve_plot[8]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
#   geom_segment(aes(x = 8.5, y = true_ve_plot[9], xend = 9.5, yend = true_ve_plot[9]),
#                color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)
# 
# boxplot_PE
# 
# ggsave(here("figures","aem", today_date,"results_boxplot_VE_scale_PE.png"), plot = boxplot_PE, width = 6, height = 4, dpi = 300)
# 
# 
# 
estimates_PK <- estimates %>% filter(country == "PK") 

true_ve_PK <- estimates_PK %>%
  group_by(scenario) %>%
  summarize(true_ve = mean(true_ve)) %>%
  ungroup()

true_ve_plot <- true_ve_PK$true_ve


boxplot_PK <- estimates_PK %>% 
  ggplot(aes(x=as.factor(VE), y=estimated_ve, fill=method)) + 
  geom_boxplot(alpha=0.2) + 
  xlab("True Simulated VE") +
  ylab("Estimated VE") +
  geom_segment(aes(x = 0.5, y = true_ve_plot[1], xend = 1.5, yend = true_ve_plot[1]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)  +
  geom_segment(aes(x = 1.5, y = true_ve_plot[2], xend = 2.5, yend = true_ve_plot[2]), 
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) + 
  geom_segment(aes(x = 2.5, y = true_ve_plot[3], xend = 3.5, yend = true_ve_plot[3]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 3.5, y = true_ve_plot[4], xend = 4.5, yend = true_ve_plot[4]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 4.5, y = true_ve_plot[5], xend = 5.5, yend = true_ve_plot[5]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 5.5, y = true_ve_plot[6], xend = 6.5, yend = true_ve_plot[6]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 6.5, y = true_ve_plot[7], xend = 7.5, yend = true_ve_plot[7]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 7.5, y = true_ve_plot[8], xend = 8.5, yend = true_ve_plot[8]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 8.5, y = true_ve_plot[9], xend = 9.5, yend = true_ve_plot[9]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)

boxplot_PK

ggsave(here("figures","aem", today_date,"results_boxplot_VE_scale_PK.png"), plot = boxplot_PK, width = 6, height = 4, dpi = 300)


  #Simplified Version

estimates_PK <- estimates %>% filter(country == "PK", method != "Weighted")

true_ve_PK <- estimates_PK %>%
  group_by(scenario) %>%
  summarize(true_ve = mean(true_ve)) %>%
  ungroup()

true_ve_plot <- true_ve_PK$true_ve


boxplot_PK <- estimates_PK %>% 
  ggplot(aes(x=as.factor(VE), y=estimated_ve, fill=method)) + 
  geom_boxplot(alpha=0.2) + 
  xlab("True Simulated VE") +
  ylab("Estimated VE") +
  geom_segment(aes(x = 0.5, y = true_ve_plot[1], xend = 1.5, yend = true_ve_plot[1]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)  +
  geom_segment(aes(x = 1.5, y = true_ve_plot[2], xend = 2.5, yend = true_ve_plot[2]), 
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) + 
  geom_segment(aes(x = 2.5, y = true_ve_plot[3], xend = 3.5, yend = true_ve_plot[3]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 3.5, y = true_ve_plot[4], xend = 4.5, yend = true_ve_plot[4]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 4.5, y = true_ve_plot[5], xend = 5.5, yend = true_ve_plot[5]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 5.5, y = true_ve_plot[6], xend = 6.5, yend = true_ve_plot[6]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 6.5, y = true_ve_plot[7], xend = 7.5, yend = true_ve_plot[7]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 7.5, y = true_ve_plot[8], xend = 8.5, yend = true_ve_plot[8]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_segment(aes(x = 8.5, y = true_ve_plot[9], xend = 9.5, yend = true_ve_plot[9]),
               color = "red", linetype = "solid", linewidth = 1, alpha = 0.5)

boxplot_PK

ggsave(here("figures","aem", today_date,"results_boxplot_VE_scale_PK_og_models.png"), plot = boxplot_PK, width = 6, height = 4, dpi = 300)




estimates_BG <- estimates_BG %>%
  mutate(MSE = (true_ve - estimated_ve)^2)

estimates_BG_by_VE <- estimates_BG %>%
  group_by(VE, method) %>%
  summarize(Average_MSE = mean(MSE))


BG_MSE_plot <- estimates_BG_by_VE %>% ggplot(aes(x=VE, y=Average_MSE, color = method)) +
  geom_line() +
  geom_point()

BG_MSE_plot

ggsave(here("figures","aem", today_date,"results_MSE_plot_VE_scale_BG.png"), plot = BG_MSE_plot, width = 6, height = 4, dpi = 300)



estimates_IN <- estimates_IN %>%
  mutate(MSE = (true_ve - estimated_ve)^2)

estimates_IN_by_VE <- estimates_IN %>%
  group_by(VE, method) %>%
  summarize(Average_MSE = mean(MSE))


IN_MSE_plot <- estimates_IN_by_VE %>% ggplot(aes(x=VE, y=Average_MSE, color = method)) +
  geom_line() +
  geom_point()

IN_MSE_plot

ggsave(here("figures","aem", today_date,"results_MSE_plot_VE_scale_IN.png"), plot = IN_MSE_plot, width = 6, height = 4, dpi = 300)



estimates_PE <- estimates_PE %>%
  mutate(MSE = (true_ve - estimated_ve)^2)

estimates_PE_by_VE <- estimates_PE %>%
  group_by(VE, method) %>%
  summarize(Average_MSE = mean(MSE))


PE_MSE_plot <- estimates_PE_by_VE %>% ggplot(aes(x=VE, y=Average_MSE, color = method)) +
  geom_line() +
  geom_point()

PE_MSE_plot

ggsave(here("figures","aem", today_date,"results_MSE_plot_VE_scale_PE.png"), plot = PE_MSE_plot, width = 6, height = 4, dpi = 300)


estimates_PK <- estimates %>% filter(country == "PK") 

estimates_PK <- estimates_PK %>%
  mutate(MSE = (true_ve - estimated_ve)^2)

estimates_PK_by_VE <- estimates_PK %>%
  group_by(VE, method) %>%
  summarize(Average_MSE = mean(MSE))


PK_MSE_plot <- estimates_PK_by_VE %>% ggplot(aes(x=VE, y=Average_MSE, color = method)) +
  geom_line() +
  geom_point() + 
  ylab("Average MSE") +
  xlab("True Simulated VE")

PK_MSE_plot

ggsave(here("figures","aem", today_date,"results_MSE_plot_VE_scale_PK.png"), plot = PK_MSE_plot, width = 6, height = 4, dpi = 300)







# coinfections_summary <- estimates %>% 
#   mutate(num_coinf = as.numeric(num_coinf),
#          num_events = as.numeric(num_events)) %>%
#   group_by(country, scenario) %>%
#   summarize(avg_co_inf = mean(num_coinf),
#             avg_event_num = mean(num_events),
#             co_inf_rate = avg_co_inf/avg_event_num)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# metrics <- estimates %>% 
#   mutate(num_events = as.numeric(num_events),
#          total_follow_up = as.numeric(total_follow_up),
#          num_coinf = as.numeric(num_coinf)) %>%
#   group_by(scenario, method) %>%
#   summarize(avg_num_events = mean(num_events),
#             avg_follow_up = mean(total_follow_up),
#             avg_coinfections = mean(num_coinf),
#             coinfection_prop = avg_coinfections/avg_num_events) %>%
#   filter(method != "Weighted")
# 
# 
# write.csv(metrics, file = here("figures","aem", "together_table1.csv"))
# 
# 
# 
# 
# 
# boxplot1 <- estimates %>% 
#   filter(scenario %in% c(1,2,3)) %>%
#   ggplot(aes(x=as.factor(scenario), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("Scenario") +
#   ylab("Estimated VE") +
#   geom_hline(yintercept = 0.4, color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)
# 
# boxplot1
# 
# boxplot2 <- estimates %>% 
#   filter(scenario %in% c(1,4,5)) %>%
#   ggplot(aes(x=as.factor(scenario), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("Scenario") +
#   ylab("Estimated VE") +
#   geom_segment(aes(x = 0.5, y = 0.4, xend = 1.5, yend = 0.4), 
#                color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)  +
#   geom_segment(aes(x = 1.5, y = 0.3, xend = 2.5, yend = 0.3), 
#                color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)  +
#   geom_segment(aes(x = 2.5, y = 0.5, xend = 3.5, yend = 0.5), 
#                color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)
# 
# boxplot2
# 
# boxplot3 <- estimates %>% 
#   filter(scenario %in% c(1,6,7)) %>%
#   ggplot(aes(x=as.factor(scenario), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("Scenario") +
#   ylab("Estimated VE") +
#   geom_hline(yintercept = 0.4, color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)
# 
# boxplot3
# 
# boxplot4 <- estimates %>% 
#   filter(scenario %in% c(1,8,9)) %>%
#   ggplot(aes(x=as.factor(scenario), y=estimated_ve, fill=method)) + 
#   geom_boxplot(alpha=0.2) + 
#   xlab("Scenario") +
#   ylab("Estimated VE") +
#   geom_hline(yintercept = 0.4, color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5)
# 
# boxplot4
# 
# 
# ggsave(here("figures","aem","results_boxplot_123.png"), plot = boxplot1, width = 6, height = 4, dpi = 300)
# ggsave(here("figures","aem","results_boxplot_145.png"), plot = boxplot2, width = 6, height = 4, dpi = 300)
# ggsave(here("figures","aem","results_boxplot_167.png"), plot = boxplot3, width = 6, height = 4, dpi = 300)
# ggsave(here("figures","aem","results_boxplot_189.png"), plot = boxplot4, width = 6, height = 4, dpi = 300)
# 
# 
# 
# 
# 
# 
# table(data$shigella_diarrhea + data$ETEC_diarrhea + data$other_diarrhea)
