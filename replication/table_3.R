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
  registerDoParallel(cores=config$cores)
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

bootstrap.c <- function(original.data, regress.fun, bootstrap.fun, original.test.results, num.resample=1000) {
  foreach(seq_len(num.resample), .combine=cbind) %dopar% {
    bootstrap.fun(original.data) %>% 
      regress.fun %>%
      left_join(original.test.results %>% 
#                   select(depvar, t.value) %>% 
#                   rename(origin.t.value=t.value),
                  select(depvar, est) %>% 
                  rename(origin.est=est),
                .,
                by="depvar") %>%
      transmute(centered.stat=abs((est - origin.est)/se))
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

round3.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round3.dta")) %>% 
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

centered.stats <- bootstrap.c(round3.data, 
                          regress.fun.factory(depvars=test.results$depvar, c("upazila", controls.r3), all.treat, cluster="cluster.id"), 
                          block.bootstrap.factory("village"),
                          test.results,
                          num.resample=1000) %>% 
  set_names(paste(names(.), seq_len(ncol(.)), sep="."))

alpha <- c(0.10, 0.05)

max.w <- foreach(start=seq_len(nrow(centered.stats)), .combine=rbind) %dopar% {
  centered.stats[seq(start, nrow(centered.stats)), ] %>% 
    summarise_each(funs(max)) 
} 

critical.pts <- max.w %>% 
  as.matrix %>% 
  aaply(1, quantile, probs=1 - alpha) %>%
  set_colnames(paste0("critical.c.", 100 * alpha)) %>% 
  as.data.frame

(test.results %<>% bind_cols(critical.pts))
