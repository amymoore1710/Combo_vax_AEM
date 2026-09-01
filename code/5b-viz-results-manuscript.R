# libraries ----
pacman::p_load(dplyr,
               tidyr,
               ggplot2,
               scales,
               purrr,
               openxlsx,
               stringr,
               ggokabeito)

# data ----
ss.df <- readRDS(here::here("clean data","sens_spec_results.RDS"))
sample.df <- readRDS(here::here("clean data","sample_sizes3.RDS"))
sample.total.df <- readRDS(here::here("clean data", "sample_sizes3_comparison.RDS"))

# main data ----
main_results <- readRDS(here::here("clean data","scenario ve results","Main_v2.RDS"))
ve_main   <- main_results$ve

# update severity variable label for severity score ----
# List of dataset names as strings
datasets <- c("ss.df", "sample.df", "ve_main", "sample.total.df")

for (d in datasets) {
  # Get the dataset from the global environment
  temp <- get(d)
  
  # Replace the string in the 'severity' column (assuming it's a character vector)
  temp$severity[temp$severity == "MAL-ED Score >=6"] <- "Score-based MSD"
  temp$severity[temp$severity == "GEMS MSD"] <- "Criteria-based MSD"
  
  temp <- temp |> 
    filter(!(endpoint %in% c("Endpoint 2B:\nShigella Ct<30.4", "Endpoint 4B:\nShigella Ct<30.4 +\nno other pathogen Ct<30"))) |> 
    mutate(endpoint = as.character(endpoint),
           endpoint = factor(endpoint,
                           levels = c("Endpoint 1:\nAny Shigella",
                                      "Endpoint 2A:\nShigella Ct<28.8",
                                      "Endpoint 3:\nAny Shigella +\nno other pathogen",
                                      "Endpoint 4A:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
                                      "Endpoint Culture:\nCulture positive"),
                           labels = c("Endpoint 1:\nAny Shigella",
                                      "Endpoint 2:\nShigella Ct<28.8",
                                      "Endpoint 3:\nAny Shigella +\nno other pathogen",
                                      "Endpoint 4:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
                                      "Endpoint 5:\nCulture positive")))

  # Assign the modified dataset back to the global environment
  assign(d, temp)
}

################################################################################
########### Figures for manuscript #############################################
################################################################################
# Figure: Sensitivity, specificity, and VE bias for all sites ----
fig2.ss <- ss.df |>
  filter(scenario == "Main" & country_id == "All Sites")

fig2.ve.bias <- ve_main |> 
  filter(scenario == "Main" & 
    country_id == "All Sites" & 
    measure == "VE Bias") |>
  mutate(measure = "Bias in observed VE")

fig2.df <- fig2.ss |> 
  add_row(fig2.ve.bias) |> 
  mutate(measure = factor(measure, levels = c("Sensitivity",
"Specificity","Bias in observed VE")))

ggplot(fig2.df,
       aes(x = endpoint, y = mean, color = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  geom_hline(data = filter(fig2.df, measure == "Bias in observed VE"),
             aes(yintercept = 0), linetype = "dashed", color = "grey40", 
             inherit.aes = FALSE) +
  geom_blank(data = data.frame(measure = factor("Sensitivity", 
                                               levels = levels(fig2.df$measure)),
                              mean = c(0, 1)),
           aes(y = mean), inherit.aes = FALSE) +
  facet_grid(measure ~ ., scales = "free", switch = "y") +  
  labs(x = "", y = NULL, 
       color = "Severity") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_linedraw() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black")) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1))

ggsave(here::here("figures","manuscript","Figure2.jpg"),width = 6, height = 6)

# Figure: Site-specific VE ----
obs_ve <- ve_main |> filter(measure == "Observed VE") |> 
  mutate(endpoint = substr(endpoint,1,10))
true_ve <- ve_main |> filter(measure == "True VE")

ggplot(obs_ve, aes(x = endpoint, y = mean)) +
  geom_hline(data = true_ve, aes(yintercept = mean),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.35) +
  facet_grid(severity ~ country_id) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "", y = "Observed VE") +
  theme_linedraw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","manuscript","Figure3.jpg"),width = 8, height = 5)

# Figure: VE bias vs sample size ----
bias.sample.total <- ve_main |>
  filter(measure == "VE Bias", country_id == "All Sites") |>
  left_join(sample.total.df %>% 
              filter(scenario == "Main") |> 
    select(endpoint,severity,scenario,ntotal_sum_of_sites), 
            by = c("endpoint","severity","scenario")) |>
  mutate(endpoint = factor(endpoint,
                           levels = c("Endpoint 1:\nAny Shigella",
                                      "Endpoint 2:\nShigella Ct<28.8",
                                      "Endpoint 3:\nAny Shigella +\nno other pathogen",
                                      "Endpoint 4:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
                                      "Endpoint 5:\nCulture positive"),
                           labels = c("Endpoint 1: Any Shigella",
                                      "Endpoint 2: Shigella Ct<28.8",
                                      "Endpoint 3: Any Shigella + no other pathogen",
                                      "Endpoint 4: Shigella Ct<28.8 + no other pathogen Ct<30",
                                      "Endpoint 5: Culture positive"))) |> 
#   group_by(severity) |> 
  mutate(min.sample = min(ntotal_sum_of_sites),
         p.diff.sample = (ntotal_sum_of_sites - min.sample)/((ntotal_sum_of_sites+min.sample)/2))

ggplot(data = bias.sample.total, aes(x = p.diff.sample,y = mean, color = endpoint))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.05,
                linewidth = 0.35, alpha = 0.8) + 
  facet_wrap(vars(severity),ncol = 1)+
  theme_linedraw() +  
  labs(x = "Percent Difference in Sample size", y = "Bias in observed VE", 
       color = "Endpoint") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_x_continuous(labels = scales::percent) +
#   ggokabeito::scale_color_okabe_ito() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"),
guide = "none")+
  theme_linedraw() + 
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","manuscript","Figure4.jpg"),height = 6, width = 7)

################################################################################
########### Calculations for manuscript text ##################################
################################################################################
# Difference in mean sensitivity/specificity between each MSD-type outcome and
# the corresponding "all diarrhea" outcome, by country and endpoint (Main
# scenario, matching Figure 2).
calc_diff <- ss.df |>
  filter(scenario == "Main") |>
  select(country_id, endpoint, severity, measure, mean) |>
  pivot_wider(names_from = severity, values_from = mean) |>
  mutate(
    `Criteria-based MSD minus All diarrhea`          = `Criteria-based MSD`           - `All diarrhea`,
    `Score-based MSD minus All diarrhea` = `Score-based MSD` - `All diarrhea`
  )

sens_rows <- calc_diff |> filter(str_detect(measure, regex("sens", ignore_case = TRUE)))
spec_rows <- calc_diff |> filter(str_detect(measure, regex("spec", ignore_case = TRUE)))

message("\nTable 1: Sensitivity, Criteria-based MSD minus All diarrhea (rows = country, cols = endpoint)")
print(sens_rows |>
  select(country_id, endpoint, `Criteria-based MSD minus All diarrhea`) |>
  arrange(endpoint) |> 
  mutate(endpoint = substr(endpoint,1,10),
         `Criteria-based MSD minus All diarrhea` = 100 * `Criteria-based MSD minus All diarrhea`) |> 
  pivot_wider(names_from = endpoint, values_from = `Criteria-based MSD minus All diarrhea`))

message("\nTable 2: Sensitivity, Score-based MSD minus All diarrhea (rows = country, cols = endpoint)")
print(sens_rows |>
  select(country_id, endpoint, `Score-based MSD minus All diarrhea`) |>
  arrange(endpoint) |> 
  mutate(endpoint = substr(endpoint,1,10),
         `Score-based MSD minus All diarrhea` = 100 * `Score-based MSD minus All diarrhea`) |> 
  pivot_wider(names_from = endpoint, values_from = `Score-based MSD minus All diarrhea`))

message("\nTable 3: Specificity, Criteria-based MSD minus All diarrhea (rows = country, cols = endpoint)")
print(spec_rows |>
  select(country_id, endpoint, `Criteria-based MSD minus All diarrhea`) |>
  arrange(endpoint) |> 
  mutate(endpoint = substr(endpoint,1,10),
         `Criteria-based MSD minus All diarrhea` = 100 * `Criteria-based MSD minus All diarrhea`) |> 
  pivot_wider(names_from = endpoint, values_from = `Criteria-based MSD minus All diarrhea`))

message("\nTable 4: Specificity, Score-based MSD minus All diarrhea (rows = country, cols = endpoint)")
print(spec_rows |>
  select(country_id, endpoint, `Score-based MSD minus All diarrhea`) |>
  arrange(endpoint) |> 
  mutate(endpoint = substr(endpoint,1,10),
         `Score-based MSD minus All diarrhea` = 100 * `Score-based MSD minus All diarrhea`) |> 
  pivot_wider(names_from = endpoint, values_from = `Score-based MSD minus All diarrhea`))

# Supplementary Figures #######################################################

# Figure: Sensitivity and specificity ----
ggplot(ss.df |>
         filter(scenario == "Main"),
       aes(x = country_id, y = mean, color = country_id, shape = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free", switch = "y") +  
  labs(x = "Country", y = NULL, 
       shape = "Severity", color = "Country") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_linedraw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black")) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1))

ggsave(here::here("figures","manuscript","Supp_Figure1.jpg"),width = 9, height = 6)

# Figure: Site-specific VE ----
obs_ve <- ve_main |> filter(measure == "Observed VE") |> 
  mutate(endpoint = substr(endpoint,1,10))
true_ve <- ve_main |> filter(measure == "True VE")

ggplot(obs_ve, aes(x = endpoint, y = mean)) +
  geom_hline(data = true_ve, aes(yintercept = mean),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.35) +
  facet_grid(severity ~ country_id) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Endpoint", y = "Observed VE") +
  theme_linedraw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","manuscript","Figure2.jpg"),width = 8, height = 5)

# Figure: VE bias ----
ggplot(ve_main |>
         filter(measure == "VE Bias"), 
       aes(x = country_id, y = mean, color = country_id, shape = severity)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2,
                linewidth = 0.35) + 
  facet_grid(. ~ endpoint) +  
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     # limits = c(-1, 1),
                     breaks = seq(-0.5,0.5,by = 0.1)) + 
  labs(x = "Country", y = "Bias in observed VE",
       shape = NULL) +
  scale_color_brewer(palette = "Dark2") +  
  theme_linedraw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black")) +
  guides(color = "none") # Hide the legend for color

ggsave(here::here("figures","manuscript","Supp_Figure3.jpg"),width = 9, height = 3.5)

# Figure: VE bias vs sample size ----
bias.sample <- ve_main |>
  filter(measure == "VE Bias" & country_id != "All Sites") |>
  left_join(sample.df %>% 
              filter(scenario == "Main"), 
            by = c("country_id","endpoint","severity","scenario")) |>
  mutate(endpoint = factor(endpoint,
                           levels = c("Endpoint 1:\nAny Shigella",
                                      "Endpoint 2:\nShigella Ct<28.8",
                                      "Endpoint 3:\nAny Shigella +\nno other pathogen",
                                      "Endpoint 4:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
                                      "Endpoint 5:\nCulture positive"),
                           labels = c("Endpoint 1: Any Shigella",
                                      "Endpoint 2: Shigella Ct<28.8",
                                      "Endpoint 3: Any Shigella + no other pathogen",
                                      "Endpoint 4: Shigella Ct<28.8 + no other pathogen Ct<30",
                                      "Endpoint 5: Culture positive"))) |> 
#   group_by(severity) |> 
  mutate(min.sample = min(ntotal),
         p.diff.sample = (ntotal - min.sample)/((ntotal+min.sample)/2))

hulls <- bias.sample %>%
  group_by(endpoint, severity) %>%
  slice(chull(p.diff.sample, mean))

ggplot(data = bias.sample, aes(x = p.diff.sample,y = mean, shape = country_id,color = endpoint))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_polygon(data = hulls, aes(group = endpoint, fill = endpoint),
             alpha = 0.2, color = NA) +
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.05,
                linewidth = 0.35, alpha = 0.8) + 
  facet_wrap(vars(severity),ncol = 1)+
  theme_linedraw() +  
  labs(x = "Percent Difference in Sample size", y = "Bias in observed VE", 
       shape = "Country", color = "Endpoint") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_x_continuous(labels = scales::percent) +
#   ggokabeito::scale_color_okabe_ito() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"),
guide = "none")+
  scale_shape_manual(values = c(16,17,15,18))+
  theme_linedraw() + 
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","manuscript","Supp_Figure4.jpg"),height = 6, width = 7)
