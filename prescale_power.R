source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

get.resampler <- function(.data, branch="branch.office", num.villages.per.branch=1) {
  inner.resampler <- function(num.cash, num.loan=num.cash, num.control=num.cash, cluster.size=50) {
    incentive.resampler <- function(.incentive, num.resample, cluster.size) {
      .data %>% 
        filter(incentive %in% .incentive) %>% 
        select_(branch) %>% # Resample over branches
        unique %>% 
        sample_n(num.resample, replace=TRUE) %>%
        mutate(resample.branch.id=seq_len(n())) %>%
        inner_join(.data, by=branch) %>% 
        filter(incentive %in% .incentive) %>%
        group_by(resample.branch.id) %>%  # Resample villages per branch
        select(village) %>%
        unique %>%
        sample_n(num.villages.per.branch, replace=TRUE) %>%
        ungroup %>% 
        mutate(resample.village.id=seq_len(n())) %>%
        inner_join(.data, by="village") %>%
        group_by(resample.branch.id) %>% # Resample households per _branch_
        sample_n(cluster.size, replace=TRUE) %>%
        mutate_each(funs(sample), fake.migrant=migrant, fake.ngo.help.migrate=ngo.help.migrate) %>% # Fake dependent variables
        ungroup %>%
        mutate(resample.hh.id=seq_len(n())) %>%
        mutate_each(funs(paste(., paste(.incentive, collapse="-"), sep="-")), matches("resample\\..+\\.id"))
    }
    
    incentive.resampler("cash", num.cash, cluster.size) %>% 
      bind_rows(incentive.resampler("credit", num.loan, cluster.size)) %>% 
      bind_rows(incentive.resampler(c("control"), num.control, cluster.size)) %>% 
      inner_join(select(., resample.branch.id) %>% 
                   unique %>% 
                   mutate(fake.incentive=sample(c(rep("cash", num.cash), rep("credit", num.loan), rep("control", num.control))), # Adding a "fake" treatment dummies to evaluate rejection rates (Type I error)
                          fake.incentivized=1*(fake.incentive != "control"), 
                          fake.cash=1*(fake.incentive == "cash")), 
                 by="resample.branch.id") 
  }
  
  return(inner.resampler)
}

resample.from.data <- round2.data %>%
  filter(!is.na(branch.office), !is.na(village))

resampler <- get.resampler(resample.from.data)

original.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), control=NULL, out.intercept="(Intercept)", cluster=c("village", "hhid"), coef=c("incentivized", "cash"))
original.regress.fun.2 <- regress.fun.factory(depvars=c("migrant"), control=NULL, out.intercept="(Intercept)", cluster=c("village", "hhid"), coef=c("info", "cash", "credit"))
original.iv.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), coef="migrant_new", iv="incentivized",  out.intercept="(Intercept)", cluster=c("village", "hhid")) 
original.impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), coef=c("incentivized", "cash"), out.intercept="(Intercept)", cluster=c("village", "hhid")) 
original.impact.regress.fun.2 <- regress.fun.factory(depvars=depvars_T1_08, controls="upazila", coef=c("credit", "cash", "info"), cluster=c("village", "hhid")) 

regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), controls=c(controls.r2, "upazila"), cluster=c("resample.branch.id", "resample.hh.id"), coef=c("incentivized", "cash")) 
impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), , coef=c("incentivized", "cash"), cluster=c("resample.branch.id", "resample.hh.id")) 
iv.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), coef="migrant", iv="incentivized", cluster=c("resample.branch.id", "resample.hh.id")) 
fake.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), controls=c(controls.r2, "upazila", "incentivized", "cash"), cluster=c("resample.branch.id", "resample.hh.id"), coef=c("fake.incentivized", "fake.cash")) 

round2.data %>% filter(!is.na(village)) %>% original.regress.fun
round2.data %>% filter(!is.na(village)) %>% original.impact.regress.fun
round2.data %>% filter(!is.na(village)) %>% original.iv.regress.fun

resample.est.data <- foreach(rep.id=1:1000, .combine=bind_rows) %dopar% {
  resample.data <- resampler(15, 15, 30, 100) 
  
  regress.fun(resample.data) %>% 
    bind_rows(fake.regress.fun(resample.data)) %>%
#    bind_rows(iv.regress.fun(resample.data)) %>%
    bind_rows(impact.regress.fun(resample.data)) %>%
    mutate(rep.id=rep.id) 
}

sim.sum <- resample.est.data %>% 
  group_by(depvar, coef) %>% 
  summarize_each(funs(mean, sd), est, starts_with("se"))

sim.reject.rate <- resample.est.data %>% 
  mutate(rej.p.value=pnorm(abs(est/se.HC2), lower.tail=FALSE) * 2) %>%
  group_by(depvar, coef) %>% 
  summarize(reject.rate=mean(p.value <= 0.05))
