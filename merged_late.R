library(dplyr)
library(foreign)
library(AER)
library(yaml)

source("util.R")

config <- yaml.load_file("local_config.yaml")

merged.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round2.dta"), convert.factors=FALSE) %>% 
  filter(hhid != 92, !is.na(village)) %>% 
  inner_join(read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round3.dta"), convert.factors=FALSE) %>% 
               filter(average_food3 < 2500) %>% 
               select(hhid, average_calorie_perday2),
             by="hhid") %>% 
  mutate_each(funs(factor), hhid, village, upazila, lit) %>%
  mutate(average_calorie_perday2=average_calorie_perday2.x + average_calorie_perday2.y)
    
controls.r2 <-c("lit", "walls_good", "monga", "dhaka_remit", "dhaka_network", "exp_total_pc_r1", "subsistencer1", "num_adltmalesr1", "num_childrenr1", "avgQ13earned", "constrainedr1", "bankedr1")

merged.data %>% 
  ivreg(formula(sprintf("average_calorie_perday2 ~ migrant_new + upazila + %s | . - migrant_new + cash + credit + info", paste(controls.r2, collapse=" + "))), data=.) %>% 
  coeftest(vcov.clx(., cluster=merged.data %>% select(village)))