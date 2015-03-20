library(plyr)
library(dplyr)
library(magrittr)
library(foreach)
library(foreign)
library(plm)
library(lmtest)
library(doParallel)
library(yaml)

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

# round2.data <- read.dta("~/Data/mobarak_seasonal_migration/Round2.dta")

block.bootstrap.factory <- function(cluster) {
  function(.data) {
    data_frame(original.cluster=.data[, cluster] %>% unique %>% sample(length(.), replace=TRUE)) %>%
      mutate(cluster.id=seq_len(n())) %>%
      plyr::rename(c("original.cluster"=cluster)) %>%
      inner_join(.data, by=cluster) 
  }
}

bootstrap.c <- function(original.data, regress.fun, bootstrap.fun, num.resample=1000) {
  foreach(seq_len(num.resample), .combine=cbind) %dopar% {
    bootstrap.fun(original.data) %>% 
      regress.fun %>% 
      extract("t.value") %>%
      abs
  } 
}

regress.fun.factory <- function(depvars, controls, cluster, coef) {
  function(.data) {
    data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do(plm(formula(sprintf("%s ~ %s", .$depvar, paste(c(coef, controls), collapse=" + "))), data=.data, model="pooling", index=cluster) %>%
           coeftest(., vcov=plm::vcovHC(., cluster="group")) %>% 
           extract(coef, , drop=FALSE) %>% 
           as.data.frame) %>%
      ungroup %>% 
      set_names(c("depvar", "est", "se", "t.value", "p.value"))
  }
}

round3.data <- read.dta("~/Data/mobarak_seasonal_migration/Round3.dta") %>% 
  filter(average_food3 < 2500)

depvars_T1_09 <- c("average_food2", "average_nonfood2", "average_exp2", "average_calorie_perday2", "average_protein_perday2")
appvars_AT4_09 <- c("average_edu_exp_kds2", "average_cloths_shoes2", "average_med_exp_f2", "average_med_exp_m2")
all.depvars <- c(depvars_T1_09, appvars_AT4_09)

controls.r3 <- c("litr1", "walls_good", "monga", "dhaka_remit", "dhaka_network", "exp_total_pcr1", "subsistencer1",
              "num_adltmalesr1", "num_childrenr1", "avgQ13earned", "constrainedr1", "bankedr1")

# all.treat <- c("cash", "credit", "info")
all.treat <- "incentivized"

test.results <- foreach(dep.var=all.depvars, .combine=rbind) %do% {
  plm(formula(sprintf("%s ~ %s + upazila", dep.var, paste(c(all.treat, controls.r3), collapse=" + "))), data=round3.data, model="pooling", index="village") %>%
    coeftest(vcov=plm::vcovHC(., cluster="group")) %>%
    extract(all.treat, , drop=FALSE) %>%
    as.data.frame %>%
    set_names(c("est", "se", "t.value", "p.value")) %>%
    mutate(depvar=dep.var,
           coef=all.treat)
} %>%
  arrange(desc(abs(t.value)))

critical.c <- bootstrap.c(round3.data, 
                          regress.fun.factory(depvars=test.results$depvar, c("upazila", controls.r3), all.treat, cluster="cluster.id"), 
                          block.bootstrap.factory("village"),
                          num.resample=500) %>% 
  set_names(paste(names(.), seq_len(ncol(.)), sep="."))

alpha <- 0.05

critical.c <- foreach(start=seq_len(nrow(critical.c))) %dopar% {
  critical.c[seq(start, nrow(critical.c)), ] %>% 
    summarise_each(funs(max)) %>% 
    quantile(probs=alpha/2)
} %>% 
  unlist

test.results %<>% cbind(critical.c)
