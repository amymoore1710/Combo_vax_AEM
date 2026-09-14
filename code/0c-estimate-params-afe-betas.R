#Analyzing the relationship between Ct and AFE
pacman::p_load(dplyr, 
               here,
               tidyr,
               dplyr,
               readr)

# data -----
tac <- read_csv(here::here("maled data/Tac_Sep2018.csv"))

# shigella -----

  #subset data to just shigella related variables
shigella_dataset <- tac %>% select("sid","country_id", "shigella_eiec", "shigella_eiec_afe") %>% filter(!is.na(shigella_eiec) & !is.na(shigella_eiec_afe))

  #plot relationship between ct and AFE
ggplot(shigella_dataset, aes(x=shigella_eiec, y=shigella_eiec_afe)) + 
  geom_point()
ggplot(shigella_dataset, aes(x=shigella_eiec, y=shigella_eiec_afe, color = country_id)) + 
  geom_point()

#Subset to just countries in Simulation
shigella_dataset <- shigella_dataset %>% filter(country_id %in% c("BG", "IN", "PE", "PK")) 
ggplot(shigella_dataset, aes(x=shigella_eiec, y=shigella_eiec_afe)) + 
  geom_point()
ggplot(shigella_dataset, aes(x=shigella_eiec, y=shigella_eiec_afe, color = country_id)) + 
  geom_point() +
  xlab("Shigella Ct") +
  ylab("Shigella AFe") +
  labs(title = "Shigella Pathogen Load versus AFe by country") +
  facet_wrap(~ country_id, ncol = 2)

  # Using the definition,
  # AFE_i = 1 - exp(-beta * (35 - ct_i))
  #So, -log(1 - AFE_i) = beta * (35 - ct_i)

  #Pathogen quantity = 35 - shig_ct
shigella_dataset$pathogen_quantity <- (35 -shigella_dataset$shigella_eiec)
  #transformed AFE is log(1 - AFE)
shigella_dataset$transformed_AFE <- log(1 - shigella_dataset$shigella_eiec_afe)

  #Plot Pathogen Quantity versus transformed AFE
ggplot(shigella_dataset, aes(x=pathogen_quantity, y=transformed_AFE)) + 
  geom_point() 

  #Estimate the slope of the relationship (assume no intercept - if no pathogen present, AFE = 0)
shigella.model <- glm(transformed_AFE ~ 0 + pathogen_quantity:country_id, 
                      data = shigella_dataset)
summary(shigella.model)

  #Extract the beta coefficient
shigella_afe_beta <- summary(shigella.model)$coefficients[,1]

shigella_dataset <- shigella_dataset %>%
  mutate(beta_afe = case_when(country_id == "BG" ~ shigella_afe_beta[1],
                              country_id == "IN" ~ shigella_afe_beta[2],
                              country_id == "PE" ~ shigella_afe_beta[3],
                              country_id == "PK" ~ shigella_afe_beta[4]),
         shigella_afe_predictions = 1 - exp(beta_afe * pathogen_quantity))

  #Check for goodness of fit
# shigella_dataset$shigella_afe_predictions <-  1 - exp(shigella_afe_beta*shigella_dataset$pathogen_quantity)
shighella_pos_subset <- shigella_dataset %>% filter(pathogen_quantity != 0) %>% select("sid", "shigella_eiec", "shigella_eiec_afe", "shigella_afe_predictions")
MSE <- mean( (shighella_pos_subset$shigella_afe_predictions - shighella_pos_subset$shigella_eiec_afe)^2 )
ggplot(shigella_dataset, aes(x = shigella_eiec_afe, y = shigella_afe_predictions, color = country_id)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "True versus estimated shigella AFe values") +
  xlab("True AFe") +
  ylab("Estimated AFe") +
  facet_wrap(~ country_id, ncol = 2)

  #Save the beta estimate to transform in the simulation
saveRDS(shigella_afe_beta, here::here("sim param data", "shigella_afe_beta.RDS"))


# ETEC -----

  #Subset to just ETEC related variables
ETEC_dataset <- tac %>% 
  select("sid","country_id", "ST_ETEC","ST_ETEC_afe") %>% 
  filter(!is.na(ST_ETEC)) %>%
  filter(!is.na(ST_ETEC_afe))%>%
  filter(ST_ETEC != 0 | ST_ETEC_afe != 0)

  #plot relationship between Ct and AFE
ggplot(ETEC_dataset, aes(x=ST_ETEC, y=ST_ETEC_afe)) + 
  geom_point()
ggplot(ETEC_dataset, aes(x=ST_ETEC, y=ST_ETEC_afe, color = country_id)) + 
  geom_point()

  #Subset to just countries of interest
ETEC_dataset <- ETEC_dataset %>% filter(country_id %in% c("BG", "IN", "PE", "PK")) 
ggplot(ETEC_dataset, aes(x=ST_ETEC, y=ST_ETEC_afe)) + 
  geom_point()
ggplot(ETEC_dataset, aes(x=ST_ETEC, y=ST_ETEC_afe, color = country_id)) + 
  geom_point() +
  xlab("ETEC Ct") +
  ylab("ETEC AFe") +
  labs(title = "ETEC Pathogen Load versus AFe by country") +
  facet_wrap(~ country_id, ncol = 2)

  # Using the definition,
  # AFE_i = 1 - exp(-beta * (35 - ct_i))
  #So, -log(1 - AFE_i) = beta * (35 - ct_i)

  #Pathogen quantity = 35 - ETEC_ct
ETEC_dataset$pathogen_quantity <- (35 -ETEC_dataset$ST_ETEC)
  #Transformed_AFE = log( 1 - ETEC_AFE)
ETEC_dataset$transformed_AFE <- log(1 - ETEC_dataset$ST_ETEC_afe)

  #Plot Pathogen Quantity versus transformed AFE
ggplot(ETEC_dataset, aes(x=pathogen_quantity, y=transformed_AFE, color = country_id)) + 
  geom_point() 

  #Estimate the slope of the relationship (assume no intercept - if no pathogen present, AFE = 0)
ETEC.model <- glm(transformed_AFE ~ 0 + pathogen_quantity:country_id, 
                      data = ETEC_dataset)
summary(ETEC.model)

  #Extract the beta coefficient
ETEC_afe_beta <- summary(ETEC.model)$coefficients[,1]

ETEC_dataset <- ETEC_dataset %>%
  mutate(beta_afe = case_when(country_id == "BG" ~ ETEC_afe_beta[1],
                              country_id == "IN" ~ ETEC_afe_beta[2],
                              country_id == "PE" ~ ETEC_afe_beta[3],
                              country_id == "PK" ~ ETEC_afe_beta[4]),
         ETEC_afe_predictions = 1 - exp(beta_afe * pathogen_quantity))


  #Check for goodness of fit
ETEC_pos_subset <- ETEC_dataset %>% 
  filter(pathogen_quantity != 0) %>% 
  select("sid", "ST_ETEC", "ST_ETEC_afe", "ETEC_afe_predictions")
MSE <- mean( (ETEC_pos_subset$ETEC_afe_predictions - ETEC_pos_subset$ST_ETEC_afe)^2 )
ggplot(ETEC_dataset, aes(x = ST_ETEC_afe, y = ETEC_afe_predictions, color = country_id)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  xlab("True AFe") +
  ylab("Estimated AFe") +
  labs(title = "True versus estimated ETEC AFe values") + 
  facet_wrap(~ country_id, ncol = 2)

  #Save the beta estimate to transform in the simulation
saveRDS(ETEC_afe_beta, here::here("sim param data", "ETEC_afe_beta.RDS"))


# rotavirus (other) -----

  #Subset to just rotavirus related variables
rotavirus_dataset <- tac %>% select("sid","country_id", "rotavirus", "rotavirus_afe") %>% filter(!is.na(rotavirus) & !is.na(rotavirus_afe))

  #plot relationship between Ct and AFE
ggplot(rotavirus_dataset, aes(x=rotavirus, y=rotavirus_afe)) + 
  geom_point()
ggplot(rotavirus_dataset, aes(x=rotavirus, y=rotavirus_afe, color = country_id)) + 
  geom_point()

  #Subset to just 1 country
rotavirus_dataset <- rotavirus_dataset %>% filter(country_id %in% c("BG", "IN", "PE", "PK")) 
ggplot(rotavirus_dataset, aes(x=rotavirus, y=rotavirus_afe)) + 
  geom_point()
ggplot(rotavirus_dataset, aes(x=rotavirus, y=rotavirus_afe, color = country_id)) + 
  geom_point() +
  xlab("Rotavirus Ct") +
  ylab("Rotavirus AFe") +
  labs(title = "Rotavirus pathogen load versus AFe by country") +
  facet_wrap(~ country_id, ncol = 2)

  # Using the definition,
  # AFE_i = 1 - exp(-beta * (35 - ct_i))
  #So, -log(1 - AFE_i) = beta * (35 - ct_i)

  #Pathogen quantity = 35 - ETEC_ct
rotavirus_dataset$pathogen_quantity <- (35 -rotavirus_dataset$rotavirus)
  #Transformed_AFE = log( 1 - ETEC_AFE)
rotavirus_dataset$transformed_AFE <- log(1 - rotavirus_dataset$rotavirus_afe)

  #Plot Pathogen Quantity versus transformed AFE
ggplot(rotavirus_dataset, aes(x=pathogen_quantity, y=transformed_AFE)) + 
  geom_point() 

  #Estimate the slope of the relationship (assume no intercept - if no pathogen present, AFE = 0)
rotavirus.model <- glm(transformed_AFE ~ 0 + pathogen_quantity:country_id, 
                  data = rotavirus_dataset)
summary(rotavirus.model)

  #Extract the beta coefficient
other_afe_beta <- summary(rotavirus.model)$coefficients[,1]

rotavirus_dataset <- rotavirus_dataset %>%
  mutate(beta_afe = case_when(country_id == "BG" ~ other_afe_beta[1],
                              country_id == "IN" ~ other_afe_beta[2],
                              country_id == "PE" ~ other_afe_beta[3],
                              country_id == "PK" ~ other_afe_beta[4]),
         rotavirus_afe_predictions = 1 - exp(beta_afe * pathogen_quantity))

  #Check for goodness of fit
# rotavirus_dataset$rotavirus_afe_predictions <-  1 - exp(other_afe_beta*rotavirus_dataset$pathogen_quantity)
rotavirus_pos_subset <- rotavirus_dataset %>% filter(pathogen_quantity != 0) %>% select("sid", "rotavirus", "rotavirus_afe", "rotavirus_afe_predictions")
MSE <- mean( (rotavirus_pos_subset$rotavirus_afe_predictions - rotavirus_pos_subset$rotavirus_afe)^2 )
ggplot(rotavirus_dataset, aes(x = rotavirus_afe, y = rotavirus_afe_predictions, color = country_id)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  xlab("True AFe") +
  ylab("Estimated AFe") +
  labs(title = "True versus estimated rotavirus AFe values") +
  facet_wrap(~ country_id)

  #Save the beta estimate to transform in the simulation
saveRDS(other_afe_beta, here::here("sim param data", "other_afe_beta.RDS"))
