library(plyr)
library(dplyr)
library(magrittr)
library(foreach)
library(foreign)
library(plm)
library(lmtest)
library(yaml)

config <- yaml.load_file("local_config.yaml")

round1.controls.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round1_Controls_Table1.dta")) 
