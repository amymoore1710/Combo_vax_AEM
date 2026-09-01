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
inc_df <- readRDS(here::here("clean data","inc_results3.RDS"))
risk_df <- readRDS(here::here("clean data","risk_results3.RDS"))
ve_df <- readRDS(here::here("clean data","ve_results3.RDS"))
sample.df <- readRDS(here::here("clean data","sample_sizes3.RDS"))

# main data ----
main_results <- readRDS(here::here("clean data","scenario ve results","Main_v2.RDS"))
inc_main  <- main_results$incidence
risk_main <- main_results$risk
ve_main   <- main_results$ve

# update severity variable label for severity score ----
# List of dataset names as strings
datasets <- c("ss.df", "inc_df", "risk_df", "ve_df", "sample.df", "inc_main", "risk_main", "ve_main")

for (d in datasets) {
  # Get the dataset from the global environment
  temp <- get(d)
  
  # Replace the string in the 'severity' column (assuming it's a character vector)
  temp$severity[temp$severity == "MAL-ED Score >=6"] <- "Severity score >=6"
  
  # Assign the modified dataset back to the global environment
  assign(d, temp)
}

# function: VE figure ----
f.ve.plot <- function(my.scenario) {
  ggplot(ve_df |>
           filter(scenario == my.scenario) |>
           filter(measure %in% c("True VE","Observed VE")) |>
           mutate(lower = ifelse(lower < 0,0,lower)), 
         aes(x = country_id, y = mean, color = country_id, shape = measure)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +  
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(width = 0.5), width = 0.2) + 
    facet_grid(severity ~ endpoint, switch = "y") +  
    labs(x = "Country", y = "Vaccine Efficacy (VE)", 
         shape = NULL, color = "Country") +
    scale_y_continuous(labels = percent_format(scale = 100),limits = c(0,1),
                       breaks = c(0,0.2,0.4,0.6,0.8,1)) + 
    scale_color_brewer(palette = "Dark2") +  
    theme_minimal() + 
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside") 
}
################################################################################
########### Figures for main scenario ##########################################
################################################################################
# Figure: Sensitivity and specificity ----
ggplot(ss.df |>
         filter(scenario == "Main"),
       aes(x = country_id, y = mean, color = country_id, shape = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free", switch = "y") +  
  labs(x = "Country", y = NULL, 
       shape = "Severity", color = "Country") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1))

ggsave(here::here("figures","sens_spec.jpg"),width = 14, height = 7)

# Figure for defense: Sensitivity and specificity ----
ggplot(ss.df |> 
         filter(scenario == "Main"),
       aes(x = country_id, y = mean, color = country_id, shape = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free", switch = "y") +  
  labs(x = "Country", y = NULL, 
       shape = "Severity", color = "Country") +
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_linedraw() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black")) 

ggsave(here::here("figures","sens_spec_defense.jpg"),width = 16, height = 5.5)

# Figure for PDVAC: sensitivity and specificity for Peru GEMS MSD ----

ggplot(ss.df |> 
         filter(scenario == "Main",
                country_id == "Peru",
                severity == "GEMS MSD",
                !str_starts(endpoint, "Endpoint 2B")) |>
         mutate(endpoint = str_replace(endpoint, "^Endpoint 2A", "Endpoint 2")),
       aes(x = endpoint, y = mean)) +
  geom_point(size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2) + 
  facet_grid(. ~ measure, scales = "free", switch = "y") +  
  labs(x = "Endpoint", y = NULL) +
  scale_y_continuous(labels = percent_format(scale = 100),
                     limits = c(0,1), breaks = seq(0,1,0.2)) + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","sens_spec_PDVAC.jpg"),width = 6, height = 4)

# Figure: All diarrhea incidence ----

ggplot(inc_main |>
         filter(severity == "All diarrhea",
                        measure %in% c("True","Observed")), 
       aes(x = country_id, y = mean, color = country_id, shape = measure)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  scale_y_continuous(breaks = seq(0,100,by=20))+
  facet_grid(vax ~ endpoint, switch = "y") +  
  labs(x = "Country", y = "Incidence per 100 child-years", 
       shape = NULL, color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","inc_all_diarrhea.jpg"),width = 14, height = 5)

# Figure: Severe diarrhea (Severity score >=6) incidence ----

ggplot(inc_main |>
         filter(severity == "Severity score >=6",
                        measure %in% c("True","Observed")), 
       aes(x = country_id, y = mean, color = country_id, shape = measure)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(vax ~ endpoint, switch = "y") +  
  labs(x = "Country", y = "Incidence per 100 child-years", 
       shape = NULL, color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","inc_sev_score_diarrhea.jpg"),width = 14, height = 5)

# Figure: Severe diarrhea (GEMS MSD) incidence ----

ggplot(inc_main |>
         filter(severity == "GEMS MSD",
                        measure %in% c("True","Observed")), 
       aes(x = country_id, y = mean, color = country_id, shape = measure)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(vax ~ endpoint, switch = "y") +  
  labs(x = "Country", y = "Incidence per 100 child-years", 
       shape = NULL, color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","inc_sev_gems_diarrhea.jpg"),width = 14, height = 5)

# Figure: VE ----
ggplot(ve_main |>
         filter(measure %in% c("True VE","Observed VE")) |>
         mutate(lower = ifelse(lower < 0,0,lower)), 
       aes(x = country_id, y = mean, color = country_id, shape = measure)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(severity ~ endpoint, switch = "y") +  
  labs(x = "Country", y = "Vaccine Efficacy (VE)", 
       shape = NULL, color = "Country") +
  scale_y_continuous(labels = percent_format(scale = 100),limits = c(0,1),
                     breaks = c(0,0.2,0.4,0.6,0.8,1)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","ve_obs_true_main.jpg"),width = 14, height = 8)


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
  theme_minimal() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") +
  guides(color = "none") # Hide the legend for color

ggsave(here::here("figures","ve_bias.jpg"),width = 16, height = 4)

# Figure for PDVAC: VE bias for Peru GEMS MSD ----

ggplot(ve_main |>
         filter(measure == "VE Bias",
                country_id == "Peru",
                severity == "GEMS MSD",
                !str_starts(endpoint, "Endpoint 2B")) |>
         mutate(endpoint = str_replace(endpoint, "^Endpoint 2A", "Endpoint 2")), 
       aes(x = endpoint, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_point(size = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2,
                linewidth = 0.35) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     limits = c(-0.5, 0.5),
                     breaks = seq(-0.5,0.5,by = 0.1)) + 
  labs(x = "Endpoint", y = "Bias in observed VE") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside")

ggsave(here::here("figures","ve_bias_PDVAC.jpg"),width = 3, height = 4)

# Figure: sample size ----
ggplot(sample.df %>% 
         filter(scenario == "Main"), aes(x = country_id, y = ntotal, fill = country_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Side-by-side bars
  geom_text(aes(label = ntotal), 
            position = position_dodge(width = 0.8), 
            vjust = -0.3, size = 2.5) +  # Text at the top of bars
  scale_y_continuous(breaks = seq(0,15000,by = 2000))+
  facet_grid(severity ~ endpoint) +  # Facet by severity
  scale_fill_brewer(palette = "Dark2") +  
  labs(x = "Country", y = "Total Trial Size", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Rotate x-axis labels

ggsave(here::here("figures","main_sample_size.jpg"),width = 10,height = 7)

# Figure: VE bias vs sample size ----
bias.sample <- ve_main |>
  filter(measure == "VE Bias") |>
  left_join(sample.df %>% 
              filter(scenario == "Main"), 
            by = c("country_id","endpoint","severity","scenario")) |>
 filter(!(endpoint %in% c("Endpoint 2A:\nShigella Ct<28.8", "Endpoint 4A:\nShigella Ct<28.8 +\nno other pathogen Ct<30"))) |> 
  mutate(endpoint = factor(endpoint,
                           levels = c("Endpoint 1:\nAny Shigella",
                                      "Endpoint 2A:\nShigella Ct<28.8",
                                      "Endpoint 2B:\nShigella Ct<30.4",
                                      "Endpoint 3:\nAny Shigella +\nno other pathogen",
                                      "Endpoint 4A:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
                                      "Endpoint 4B:\nShigella Ct<30.4 +\nno other pathogen Ct<30",
                                      "Endpoint Culture:\nCulture positive"),
                           labels = c("Endpoint 1: Any Shigella",
                                      "Endpoint 2A: Shigella Ct<28.8",
                                      "Endpoint 2B: Shigella Ct<30.4",
                                      "Endpoint 3: Any Shigella + no other pathogen",
                                      "Endpoint 4A: Shigella Ct<28.8 + no other pathogen Ct<30",
                                      "Endpoint 4B: Shigella Ct<30.4 + no other pathogen Ct<30",
                                      "Endpoint Culture: Culture positive"))) |> 
#   group_by(severity) |> 
  mutate(min.sample = min(ntotal),
         p.diff.sample = (ntotal - min.sample)/((ntotal+min.sample)/2))

ggplot(data = bias.sample, aes(x = p.diff.sample,y = mean, shape = country_id,color = endpoint))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
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
  scale_shape_manual(values = c(16,17,15,18))+
  theme_linedraw() + 
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black"))

ggsave(here::here("figures","ve_bias_vs_sample.jpg"),height = 5, width = 7)
################################################################################
########### Summarize across scenarios #########################################
################################################################################

# format scenarios across datasets ----
dfs <- list(ss.df, inc_df, risk_df, ve_df, sample.df)

dfs <- map(dfs, ~ .x %>% 
             mutate(scenario = factor(scenario,
                                      levels = c("Main", 
                                                 "VE.inf0.0", "VE.inf0.2", 
                                                 "VE.dis0.3", "VE.dis0.5", 
                                                 "VE.sev0.0", "VE.sev0.3", 
                                                 "VE.ct0.0", "VE.ct0.2"),
                                      labels = c("Main",
                                                 "VEinf = 0%","VEinf = 20%",
                                                 "VEdis = 30%","VEdis = 50%",
                                                 "VEmsd = 40%",
                                                 "VEmsd = 80%",
                                                 "0% reduction\nin pathogen quantity",
                                                 "20% reduction\nin pathogen quantity"))))

# Unpack the list back into separate objects
ss.df <- dfs[[1]]
inc_df <- dfs[[2]]
risk_df <- dfs[[3]]
ve_df <- dfs[[4]]
sample.df <- dfs[[5]]

# Figures: VE for varying VE against infection and disease ----

f.ve.plot("VEinf = 0%")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_veinf0.jpg"),width = 14, height = 8)

f.ve.plot("VEinf = 20%")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_veinf20.jpg"),width = 14, height = 8)

f.ve.plot("VEdis = 30%")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_vedis30.jpg"),width = 14, height = 8)

f.ve.plot("VEdis = 50%")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_vedis50.jpg"),width = 14, height = 8)

f.ve.plot("0% reduction\nin severity")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_vemsd40.jpg"),width = 14, height = 8)

f.ve.plot("30% reduction\nin severity")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_vemsd80.jpg"),width = 14, height = 8)

f.ve.plot("0% reduction\nin pathogen quantity")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_ct0.jpg"),width = 14, height = 8)

f.ve.plot("20% reduction\nin pathogen quantity")
ggsave(here::here("figures","ve across scenarios","ve_obs_true_ct20.jpg"),width = 14, height = 8)

# Tables: Main incidence and risk ----
inc.tbl <- inc_df %>% 
  filter(scenario == "Main" &
           measure == "True") %>% 
  mutate(val = paste0(sprintf("%.1f",mean),"\n(",
                      sprintf("%.1f",lower),", ",
                      sprintf("%.1f",upper),")")) %>% 
  select(severity,country_id,endpoint,vax,val) %>% 
  distinct() %>% 
  arrange(severity,country_id,endpoint,vax) %>% 
  pivot_wider(names_from = c("endpoint","vax"),values_from = "val")

risk.tbl <- risk_df %>% 
  filter(scenario == "Main" &
           measure == "True") %>% 
  mutate(val = paste0(sprintf("%.1f",100*mean),"%\n(",
                      sprintf("%.1f",100*lower),", ",
                      sprintf("%.1f",100*upper),")")) %>% 
  select(severity,country_id,endpoint,vax,val) %>% 
  distinct() %>% 
  arrange(severity,country_id,endpoint,vax) %>% 
  pivot_wider(names_from = c("endpoint","vax"),values_from = "val")

# Tables: Sensitivity and specificity ----
main.ss.tbl <- ss.df %>% 
  filter(scenario == "Main") %>% 
  mutate(val = paste0(sprintf("%.1f",100*mean),"%\n(",
                        sprintf("%.1f",100*lower),", ",
                      sprintf("%.1f",100*upper),")")) %>% 
  select(severity,country_id,endpoint,measure,val) %>% 
  distinct() %>% 
  arrange(severity,country_id,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "val")

ss.tbl <- ss.df %>% 
  group_by(severity,scenario,endpoint,measure) %>% 
  mutate(range = paste0(sprintf("%.1f",100*min(mean))," - ",
                        sprintf("%.1f",100*max(mean)))) %>% 
  select(severity,scenario,endpoint,measure,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "range")

ss.tbl2 <- ss.df %>% 
  group_by(severity,endpoint,measure) %>% 
  mutate(range = paste0(sprintf("%.1f",100*min(mean))," - ",
                        sprintf("%.1f",100*max(mean)))) %>% 
  select(severity,endpoint,measure,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "range")

# Figures: Sensitivity and specificity ----
ggplot(data = ss.df %>% filter(severity == "All diarrhea"),
       aes(x = scenario, y = mean, color = country_id))+
  geom_point(position = position_dodge(width = 0), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free") +  
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  labs(x = "Vaccine scenario", y = "Bias in observed VE",
       color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here::here("figures","all_diar_sens_spec_scenarios.jpg"),width = 18,height = 8)

ggplot(data = ss.df %>% filter(severity == "Severity score >=6"),
       aes(x = scenario, y = mean, color = country_id))+
  geom_point(position = position_dodge(width = 0), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free") +  
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  labs(x = "Vaccine scenario", y = "Bias in observed VE",
       color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here::here("figures","score_msd_sens_spec_scenarios.jpg"),width = 18,height = 8)

ggplot(data = ss.df %>% filter(severity == "GEMS MSD"),
       aes(x = scenario, y = mean, color = country_id))+
  geom_point(position = position_dodge(width = 0), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free") +  
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  labs(x = "Vaccine scenario", y = "Bias in observed VE",
       color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here::here("figures","gems_msd_sens_spec_scenarios.jpg"),width = 18,height = 8)

# Tables: VE ----
main.ve.tbl <- ve_df %>% 
  filter(scenario == "Main") %>% 
  filter(measure == "Observed VE" | 
           (measure == "True VE" & endpoint == "Endpoint 1:\nAny Shigella")) %>% 
  mutate(measure = factor(measure,levels = c("True VE","Observed VE"))) %>% 
  mutate(val = ifelse(measure == "True VE",
                      paste0(sprintf("%.1f",100*mean),"%"),
                      paste0(sprintf("%.1f",100*mean),"%\n(",
                             sprintf("%.1f",100*lower),", ",
                             sprintf("%.1f",100*upper),")"))) %>% 
  select(severity,country_id,endpoint,measure,val) %>% 
  distinct() %>% 
  arrange(severity,country_id,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "val")

ve.tbl <- ve_df %>% 
  filter(measure == "Observed VE" | 
           (measure == "True VE" & endpoint == "Endpoint 1:\nAny Shigella")) %>% 
  mutate(measure = factor(measure,levels = c("True VE","Observed VE"))) %>% 
  group_by(severity,scenario,endpoint,measure) %>% 
  mutate(range = ifelse(measure == "True VE",
                        sprintf("%.1f",100*min(mean)),
                        paste0(sprintf("%.1f",100*min(mean))," - ",
                               sprintf("%.1f",100*max(mean))))) %>% 
  select(severity,scenario,endpoint,measure,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "range")

ve.tbl2 <- ve_df %>% 
  filter(measure == "Observed VE" | 
           (measure == "True VE" & endpoint == "Endpoint 1:\nAny Shigella")) %>% 
  mutate(measure = factor(measure,levels = c("True VE","Observed VE"))) %>% 
  group_by(severity,endpoint,measure) %>% 
  mutate(range = paste0(sprintf("%.1f",100*min(mean))," - ",
                        sprintf("%.1f",100*max(mean)))) %>% 
  select(severity,endpoint,measure,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint,measure) %>% 
  pivot_wider(names_from = c("endpoint","measure"),values_from = "range")

# Tables: VE bias ----
main.ve.bias.tbl <- ve_df %>% 
  filter(scenario == "Main") %>% 
  filter(measure == "VE Bias") %>% 
  mutate(val = paste0(sprintf("%.1f",100*mean),"%\n(",
                      sprintf("%.1f",100*lower),", ",
                      sprintf("%.1f",100*upper),")")) %>% 
  select(severity,country_id,endpoint,val) %>% 
  distinct() %>% 
  arrange(severity,country_id,endpoint) %>% 
  pivot_wider(names_from = "endpoint",values_from = "val")

ve.bias.tbl <- ve_df %>% 
  filter(measure == "VE Bias") %>% 
  group_by(severity,scenario,endpoint) %>% 
  mutate(range = paste0(sprintf("%.1f",100*min(mean))," - ",
                        sprintf("%.1f",100*max(mean)))) %>% 
  select(severity,scenario,endpoint,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint) %>% 
  pivot_wider(names_from = "endpoint",values_from = "range")

ve.bias.tbl2 <- ve_df %>% 
  filter(measure == "VE Bias") %>% 
  group_by(severity,endpoint,) %>% 
  mutate(range = paste0(sprintf("%.1f",100*min(mean))," - ",
                        sprintf("%.1f",100*max(mean)))) %>% 
  select(severity,endpoint,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint) %>% 
  pivot_wider(names_from = "endpoint",values_from = "range")

# Table: sample size ----
sample.tbl <- sample.df %>% 
  group_by(severity,scenario,endpoint) %>% 
  mutate(range = paste0(min(ntotal)," - ",
                        max(ntotal))) %>% 
  select(severity,scenario,endpoint,range) %>% 
  distinct() %>% 
  arrange(severity,endpoint,scenario) %>% 
  pivot_wider(names_from = "endpoint",values_from = "range")

# Figure: summarize ve bias ----
ggplot(data = ve_df %>% filter(measure == "VE Bias"),
       aes(x = scenario, y = mean, color = country_id))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(severity ~ endpoint) +  
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     breaks = seq(-0.6,0.6,by=0.2)) + 
  labs(x = "Vaccine scenario", y = "Bias in observed VE",
       color = "Country") +
  scale_color_brewer(palette = "Dark2") +  
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here::here("figures","ve_bias_scenarios.jpg"),width = 18,height = 8)

# save all tables ----
wb <- createWorkbook()
addWorksheet(wb, "Main true inc")
writeData(wb, "Main true inc", inc.tbl)
addWorksheet(wb, "Main true risk")
writeData(wb, "Main true risk", risk.tbl)
addWorksheet(wb, "Main sens spec")
writeData(wb, "Main sens spec", main.ss.tbl)
addWorksheet(wb, "Sens spec")
writeData(wb, "Sens spec", ss.tbl)
addWorksheet(wb, "Sens spec sum.")
writeData(wb, "Sens spec sum.", ss.tbl2)
addWorksheet(wb, "Main VE")
writeData(wb, "Main VE", main.ve.tbl)
addWorksheet(wb, "VE")
writeData(wb, "VE", ve.tbl)
addWorksheet(wb, "VE sum.")
writeData(wb, "VE sum.", ve.tbl2)
addWorksheet(wb, "Main VE bias")
writeData(wb, "Main VE bias", main.ve.bias.tbl)
addWorksheet(wb, "VE bias")
writeData(wb, "VE bias", ve.bias.tbl)
addWorksheet(wb, "VE bias sum.")
writeData(wb, "VE bias sum.", ve.bias.tbl2)
addWorksheet(wb, "Sample size")
writeData(wb, "Sample size", sample.tbl)

saveWorkbook(wb, here::here("tables","vax_sim_tables.xlsx"), 
             overwrite = TRUE)