##############################################################################
##############################################################################
##############################################################################
########################### CODE FOR PART 3 ASSGN. ###########################
##############################################################################
##############################################################################
##############################################################################

# Preliminary --------------------------------------------------------------

# import the libraries
library(tidyverse)
library(haven)
library(labelled)
library(fixest)

# set the environment
setwd("~/Documents/Tilburg/classes/Labor Econ./part3_empirical")

# import data
imm   <- read_dta("origin_data.dta")
crime <- read_dta("Immigration_Crime_European_Regions_2002_2017.dta")



# Data Issues -------------------------------------------------------------

visdat::vis_miss(select(crime, empl, migr, ends_with('rate')))


# Data Preparation: dta ---------------------------------------------------

crime1 <- 
  crime %>% 
  to_factor() %>% 
  mutate(region = str_squish(region)) %>% 
  rename(nuts = id) %>% 
  relocate(migr, .before = population)

imm1 <- 
  imm %>% 
  select(-c(country, nameNUTS)) %>% 
  separate_rows(nuts, sep = "-") 

imm2 <- 
  imm1 %>% 
  filter(year == 2000) %>% 
  group_by(nuts, native) %>% 
  summarise(sum_value = sum(value)) %>% 
  pivot_wider(
    names_from   = "native", 
    names_prefix = "native", values_from = sum_value) %>% 
  ungroup() %>% 
  transmute(nuts, 
            migr_01 = native0 / (native0 + native1),
            migr_01 = migr_01 * 100)

imm2a <- imm2 %>% rename(migr = migr_01)

imm3 <- 
  imm1 %>% 
  filter(year == 1990) %>% 
  group_by(nuts, native) %>% 
  summarise(sum_value = sum(value)) %>% 
  pivot_wider(
    names_from   = "native", 
    names_prefix = "native", values_from = sum_value) %>% 
  ungroup() %>% 
  transmute(nuts,
            migr_90 = native0 / (native0 + native1),
            migr_90 = migr_90 * 100)

diffed <- 
  bind_rows(crime1, imm2a) %>% 
  arrange(nuts, year) %>% 
  mutate(
    d_MIGR = migr               - lag(migr),
    d_lpop = (log(population)    - lag(log(population)))*100,
    d_hom  = homicide_rate      - lag(homicide_rate),
    d_vehc = vehicle_theft_rate - lag(vehicle_theft_rate),
    d_imgr = d_MIGR - d_lpop,
    nuts = as.factor(nuts)
  ) %>% 
  filter(!is.na(d_hom),!is.na(d_vehc),
         !is.na(d_hom),!is.na(d_imgr),
         !is.na(d_MIGR), !is.na(d_lpop),
         !is.na(empl), !is.na(GDP3))

dta <- 
  diffed %>% 
  filter(year != 2002) %>% 
  left_join(x=.,y=imm2,by="nuts") %>% 
  left_join(x=.,y=imm3,by="nuts") 

# Data Preparation: dta2 --------------------------------------------------

tmp_id_keep <- 
  imm1 %>% 
  mutate(tmp_id = str_c(nuts,cob_english,sep="-")) %>% 
  count(tmp_id) %>% 
  filter(n==3) %>% 
  pull(tmp_id)

imm5 <- 
  imm1 %>% 
  mutate(tmp_id = str_c(nuts,cob_english,sep="-")) %>% 
  filter(tmp_id %in% tmp_id_keep) %>% 
  select(-tmp_id)

ivdt <- 
  imm5 %>%
  filter(year == 1990) %>% 
  mutate(value = value + 1) %>% 
  group_by(nuts) %>% 
  summarise(sum_val = sum(value))

imm6 <- 
  imm5 %>% 
  filter(year == 1990) %>% 
  left_join(x=.,y=ivdt,by='nuts') %>% 
  transmute(nuts, cob_english,theta = value / sum_val)

ivdt2 <- 
  imm5 %>%
  filter(year != 1990) %>% 
  group_by(nuts, year) %>% 
  summarise(sum_val = sum(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = year, names_prefix = "y", values_from = sum_val) %>% 
  transmute(nuts,lpop_diff = log(y2010) - log(y2000))

imm7 <- 
  imm5 %>% 
  filter(year != 1990) %>% 
  arrange(nuts, cob_english, year) %>% 
  mutate(diff_val = log1p(value) - log1p(lag(value)),
         diff_val = diff_val * 100) %>% 
  filter(year == 2010) %>% 
  select(nuts, cob_english, diff_val) %>% 
  left_join(x=.,y=imm6,by=c('nuts','cob_english')) %>% 
  filter(!is.nan(diff_val), !is.nan(theta)) %>%
  mutate(in_sum = diff_val * theta) %>% 
  group_by(nuts) %>% 
  summarise(d_migr_hat = sum(in_sum))

pas_fr <- 
  crime1 %>% 
  filter(year %in% c(2002,2011)) %>% 
  mutate(
    nuts,
    d_hom  = homicide_rate      - lag(homicide_rate),
    d_vehc = log1p(vehicle_theft_rate) - log1p(lag(vehicle_theft_rate)),
    d_migr = migr - lag(migr) - (log(population) - log(lag(population))),
  ) %>% 
  filter(year != 2002) %>% 
  select(nuts, d_hom, d_vehc, d_migr, foreign, GDP3, empl)

dta2 <- 
  imm7 %>% 
  inner_join(pas_fr,by='nuts')


# Estimations -------------------------------------------------------------

ols_hom_cont <- lm(d_hom ~ d_migr + log(foreign) + log(GDP3) + log(empl), data = dta2)
ols_vehc_cont <- lm(d_vehc ~ d_migr + log(foreign) + log(GDP3) + log(empl), data = dta2)
ols_hom_hat_cont <- lm(d_hom ~ d_migr_hat + log(foreign) + log(GDP3) + log(empl), data = dta2)
ols_vehc_hat_cont <- lm(d_vehc ~ d_migr_hat + log(foreign) + log(GDP3) + log(empl), data = dta2)

fe_ols_hom <- feols(
  fml      = d_hom ~ d_MIGR + d_lpop + log(empl) + I(log(empl)^2) + log(GDP3) + I(log(GDP3)^2) + log(population) + I(log(population)^2) | nuts + country + year, 
  data     = dta, 
  cluster  = ~ nuts,
  panel.id = ~ nuts + year
)
fe_ols_vehc <- feols(
  fml      = d_vehc ~ d_MIGR + d_lpop + log(empl) + I(log(empl)^2) + log(GDP3) + I(log(GDP3)^2) + log(population) + I(log(population)^2) | nuts + country + year, 
  data     = dta, 
  cluster  = ~ nuts,
  panel.id = ~ nuts + year
)


# Export the Estimates ----------------------------------------------------

etable(fe_ols_hom, fe_ols_vehc,tex = T,fitstat = c('wald','wald.p'))
stargazer::stargazer(ols_hom_cont, ols_hom_hat_cont, ols_vehc_cont, ols_vehc_hat_cont)

