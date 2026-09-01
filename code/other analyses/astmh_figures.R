# libraries ----
pacman::p_load(dplyr,ggplot2,scales,pwr)

# bangladesh data ----
bang_ss <- readRDS(here::here("sim data","sens_spec_results.RDS")) %>% 
  filter(country_id == "Bangladesh")
bang_ve <- readRDS(here::here("sim data","ve_results.RDS")) %>% 
  filter(country_id == "Bangladesh")
bang_inc <- readRDS(here::here("sim data","inc_results.RDS")) %>% 
  filter(country_id == "Bangladesh")

# sens/spec figure ----
ggplot(bang_ss, aes(x = severity, y = mean, color = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(measure ~ endpoint, scales = "free", switch = "y") +  
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("All diarrhea" = "black",
                                "Severe diarrhea" = "red"))+
  scale_y_continuous(labels = percent_format(scale = 100)) + 
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","astmh","bangladesh_sens_spec.jpg"),width = 7, height = 5)

# VE figure ----
ggplot(bang_ve |>
         mutate(lower = ifelse(lower < 0,0,lower)) |>
         filter(measure %in% c("True VE","Observed VE")) |>
         mutate(measure = factor(measure,
                                 levels = c("True VE","Observed VE"))), 
       aes(x = measure, y = mean, color = measure)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(severity ~ endpoint, switch = "y") +  
  labs(x = NULL, y = "Vaccine Efficacy (VE)") +
  scale_y_continuous(labels = percent_format(scale = 100),limits = c(0,1),
                     breaks = c(0,0.2,0.4,0.6,0.8,1)) + 
  scale_color_brewer(palette = "Dark2") +  
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 

ggsave(here::here("figures","astmh","bangladesh_ve.jpg"),width = 7, height = 5)

# VE bias figure ----
ggplot(bang_ve |>
         filter(measure == "VE Bias"), 
       aes(x = severity, y = mean, color = severity)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.2) + 
  facet_grid(. ~ endpoint) +  
  scale_y_continuous(labels = percent_format(scale = 100),limits = c(-0.9,0.9),
                     breaks = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8)) + 
  labs(x = NULL, y = "Bias in observed VE") +
  scale_color_manual(values = c("All diarrhea" = "black",
                                "Severe diarrhea" = "red"))+
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside") 
  
ggsave(here::here("figures","astmh","bangladesh_ve_bias.jpg"),width = 7, height = 5)

# power calculation ----

bang_power <- bang_inc %>% 
  #filter(measure == "Observed") %>% 
  mutate(mean = mean/1200,
         risk = 1-exp(-mean*12)) %>% 
  select(-lower,-upper,-mean) %>% 
  pivot_wider(names_from = "vax", values_from = "risk") %>% 
  mutate(measure = factor(measure,
                          levels = c("True","Observed")))

bang_true <- bang_power %>% 
  filter(measure == "True") %>% 
  mutate(endpoint = "True Shigella\ndiarrhea") %>% 
  select(-measure) %>% 
  distinct()

bang_power <- bang_power %>% 
  filter(measure == "Observed") %>% 
  select(-measure) %>% 
  add_row(bang_true)
  

# Apply power.prop.test to each row using mapply
bang_power$sample_size <- mapply(function(p1, p2) {
  power.prop.test(p1 = p1, p2 = p2, power = 0.8, sig.level = 0.05)$n
}, bang_power$Placebo, bang_power$Vaccine)

# From pwr gives very similar results
# https://rpubs.com/mbounthavong/sample_size_power_analysis_R
bang_power$sample_size2 <- mapply(function(p1, p2) {
  pwr.2p.test(h = ES.h(p1 = p1, p2 = p2), sig.level = 0.05, power = .80)$n
}, bang_power$Placebo, bang_power$Vaccine)

bang_power$total_sample <- bang_power$sample_size*2

bang_power$endpoint <- factor(bang_power$endpoint,
                              levels = c("True Shigella\ndiarrhea",
                                         "Endpoint 1:\nAny Shigella",
                                         "Endpoint 2A:\nShigella Ct<28.8",
                                         "Endpoint 2B:\nShigella Ct<30.4",
                                         "Endpoint 3:\nAny Shigella +\nno other pathogen"))

# From epiR gives different results
library(epiR)
epi.sscohortt(irexp1 = 0.0493, irexp0 = 0.0805, FT = 12, n = NA, power = 0.80, r = 1,
              design = 1, sided.test = 1, nfractional = FALSE, conf.level = 0.95)$n.total

epi.sscohortt(irexp1 = 0.0481159174, irexp0 = 0.07732515, FT = 12, n = NA, power = 0.80, r = 1,
              design = 1, sided.test = 1, nfractional = FALSE, conf.level = 0.95)$n.total

epi.sscohortc(irexp1 = 0.4466, irexp0 = 0.6193, pexp = 0.5, n = NA, power = 0.80, r = 1,
              design = 1, sided.test = 1, nfractional = FALSE, conf.level = 0.95)$n.total

# sample size figure ----
ggplot(bang_power, aes(x = endpoint, y = total_sample)) +
  geom_bar(stat = "identity", width = 0.6) +  # Skinny bars with dodge
  geom_text(aes(label = round(total_sample)), 
            vjust = -0.5, size = 3) +  # Text on top of the bars, slightly above the bar
  facet_grid(severity ~ .) +  # Separate top and bottom panels by severity
  labs(x = NULL, y = "Number of participants",
       fill = "Incidence") +
  ylim(0,3300)+
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "bottom")  # Adjust legend position

ggsave(here::here("figures","astmh","bangladesh_power.jpg"),width = 6.5, height = 5)

# parameters ----
## IR of shigella and other diarrhea per child month
diarrhea_ir <- readRDS(here::here("sim param data","diarrhea_incidence.RDS"))
## Probability of shigella subclinical infection per child month
shigella_sub <- readRDS(here::here("sim param data","shigella_subclinical_prob.RDS"))
## Probability of other pathogen infection per child month
other_sub <- readRDS(here::here("sim param data","other_subclinical_prob.RDS"))
## Mean and SD of severity score for shigella and other diarrhea
severity <- readRDS(here::here("sim param data","severity.RDS"))
## Mean and SD of CT for shigella diarrhea and subclincial infections
shigella_ct <- readRDS(here::here("sim param data","shigella_ct.RDS"))
## Mean and SD of CT for other infections 
other_ct <- readRDS(here::here("sim param data","other_ct.RDS"))

# Vaccine efficacy assumptions 
ve_infection <- 0.10  # 10% VE against infection
ve_disease <- 0.40  # 40% VE against disease
ve_severe_dis <- 0.60 # 60% VE against severe disease
ve_severity <- 1-((1-ve_severe_dis)/((1-ve_infection)*(1-ve_disease))) # VE against severe symptoms, given diarrhea (from Avnika)
ve_ct <- 0.10  # reduction in pathogen quantity among breakthrough cases

# Vaccinated parameters + filter for bangladesh
bg_diarrhea_ir <- diarrhea_ir %>% 
  mutate(vax_IR_shigella = IR_shigella * (1-ve_disease)) %>% 
  filter(country_id == "BG")
bg_shigella_sub <- shigella_sub %>% 
  left_join(bg_diarrhea_ir,by=c("country_id","agegrp")) %>% 
  mutate(vax_predicted_prob = (predicted_prob + (IR_shigella - vax_IR_shigella)) * (1-ve_infection)) %>% 
  filter(country_id == "BG")
bg_severity <- severity %>% 
  mutate(vax_mean_sev = ifelse(type == "Shigella",
                               mean_sev * (1-ve_severity),
                               mean_sev)) %>% 
  filter(country_id == "BG")
bg_shigella_ct <- shigella_ct %>% 
  mutate(vax_mean_ct = mean_ct * (1+ve_ct)) %>% # does this make sense???
  mutate(vax_mean_ct = ifelse(vax_mean_ct > 35, 34.9,vax_mean_ct)) %>% 
  filter(country_id == "BG")
bg_other_sub <- other_sub %>% 
  filter(country_id == "BG")
bg_other_ct <- other_ct %>% 
  filter(country_id == "BG")
