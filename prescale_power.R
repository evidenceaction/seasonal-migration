source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

get.cluster.resampler <- function(.data, cluster.level="village", num.village.per.cluster=1) {
  inner.resampler <- function(num.cash, num.loan=num.cash, num.control=num.cash) {
    resample.cluster.col.name <- sprintf("resample.%s.id", cluster.level)
    
    incentive.resampler <- function(.incentive, num.resample, cluster.size) {
      .data %>% 
        filter(incentive %in% .incentive) %>% 
        select_(cluster.level) %>% 
        unique %>% 
        sample_n(num.resample, replace=TRUE) %>% {
          .[[resample.cluster.col.name]] <- seq_len(nrow(.))
          return(.)
        } %>%
        inner_join(.data, by=cluster.level) %>% {
          vill.ids <- . %$% village %>% unique
          if (length(vill.ids) > num.village.per.cluster) {
            filter(., village %in% sample(vill.ids, num.village.per.cluster, replace=TRUE)) 
          } else {
            return(.)
          }
        } %>%
        group_by_(resample.cluster.col.name) %>% 
        mutate_each(funs(sample), fake.migrant=migrant, fake.ngo.help.migrate=ngo.help.migrate) %>% # Fake dependent variables
        ungroup %>% 
        mutate(resample.hh.id=seq_len(n())) %>%
        mutate_each(funs(paste(., paste(.incentive, collapse="-"), sep="-")), matches("resample\\..+\\.id"))
    }
    
    incentive.resampler("cash", num.cash) %>% 
      bind_rows(incentive.resampler("credit", num.loan)) %>% 
      bind_rows(incentive.resampler(c("control"), num.control)) %>% 
      inner_join(select_(., resample.cluster.col.name) %>% 
                   unique %>% 
                   mutate(fake.incentive=sample(c(rep("cash", num.cash), rep("credit", num.loan), rep("control", num.control))), # Adding a "fake" treatment dummies to evaluate rejection rates (Type I error)
                          fake.incentivized=1*(fake.incentive != "control"), 
                          fake.cash=1*(fake.incentive == "cash")), 
                 by=resample.cluster.col.name)
  }
  
  return(inner.resampler)
}

get.branch.office.resampler <- function(.data, branch="branch.office", num.villages.per.branch=1) {
  inner.resampler <- function(num.cash, num.loan=num.cash, num.control=num.cash, cluster.size=NULL) {
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
        { if (!is.null(cluster.size)) sample_n(., cluster.size, replace=TRUE) else . } %>%
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

get.stratified.branch.resampler <- function(.data, branch="branch.office") {
  function(num.resample, num.village.pairs.per.branch=1, total.num.villages=num.resample*num.village.pairs.per.branch*2, village.sample.size=NULL) {
    resample.villages <- function(.resampled.data, num.pairs=1, start.id=1) {
      .resampled.data %>% 
        inner_join(.data, by=branch) %>% 
        group_by(resample.branch.id) %>%  # Resample villages per branch
        select(village) %>%
        unique %>%
        sample_n(num.pairs * 2, replace=TRUE) %>%
        ungroup %>% 
        mutate(resample.village.id=seq_len(n()) + start.id - 1) %>%
        inner_join(.data, by="village") %>% 
        group_by(resample.village.id) %>% 
        { if (!is.null(village.sample.size)) sample_n(., village.sample.size, replace=TRUE) else . } %>%
        mutate(resample.hh.id=paste(resample.village.id, seq_len(n()), sep="-")) %>% 
        ungroup
    }
    
    handle.deficit.villages <- function(.resampled.data) {
      num.resample.villages <- n_distinct(.resampled.data$resample.village.id)
      village.deficit <- total.num.villages - num.resample.villages
        
      if (village.deficit > 0) {
        deficit.pairs <- min(village.deficit %/% 2, num.resample)
        
        .resampled.data %>% 
          select_("resample.branch.id", branch) %>% 
          distinct(resample.branch.id) %>% 
          sample_n(deficit.pairs, replace=FALSE) %>% 
          resample.villages(start.id=num.resample.villages + 1) %>% 
          bind_rows(.resampled.data) %>% 
          handle.deficit.villages # Recursive
      } else {
        return(.resampled.data)
      }
    }
  
    randomize.fake.branch.treat <- . %>% 
        select(resample.branch.id) %>% 
        unique %>% 
        mutate(fake.branch.treat=sample(rep(0:1, each=n_distinct(resample.branch.id) %/% 2))) 
    
    randomize.fake.village.treat <- . %>% 
        group_by(resample.branch.id) %>% # Resample households per _branch_
        select(resample.village.id) %>% 
        distinct(resample.village.id) %>% 
        mutate(fake.treat=sample(rep(0:1, each=n_distinct(resample.village.id) %/% 2))) %>% 
        ungroup 
    
    .data %>% 
      select_(branch) %>% # Resample over branches
      unique %>% 
      sample_n(num.resample, replace=TRUE) %>%
      mutate(resample.branch.id=seq_len(n())) %>%
      resample.villages(num.village.pairs.per.branch) %>%  # Get the min number of villages for each branch
      handle.deficit.villages %>% # Do we still need to find more villages in the selected branches?
      inner_join(randomize.fake.village.treat(.), by=c("resample.branch.id", "resample.village.id")) %>% 
      inner_join(randomize.fake.branch.treat(.), by="resample.branch.id") %>% 
      mutate(fake.branch.treat=fake.branch.treat * fake.treat) %>% 
      mutate(resample.hh.id=seq_len(n())) %>% 
      mutate_each(funs(factor), matches("^resample\\..+\\.id$"))
  }
}

resample.from.data <- round2.data %>%
  filter(!is.na(branch.office), !is.na(village)) %>% 
  group_by(branch.office) %>% 
  mutate(num.branch.vill=n_distinct(village)) %>% 
  ungroup %>% 
  filter(num.branch.vill > 1)

branch.office.resampler <- get.branch.office.resampler(resample.from.data)
village.resampler <- get.cluster.resampler(resample.from.data, "village")
union.resampler <- get.cluster.resampler(resample.from.data, "union")

original.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), control=NULL, out.intercept="(Intercept)", cluster=c("village", "hhid"), coef=c("incentivized", "cash"))
original.regress.fun.2 <- regress.fun.factory(depvars=c("migrant"), control=NULL, out.intercept="(Intercept)", cluster=c("village", "hhid"), coef=c("info", "cash", "credit"))
original.iv.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), coef="migrant_new", iv="incentivized",  out.intercept="(Intercept)", cluster=c("village", "hhid"), se.types="wild") 
original.impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, controls=c(controls.r2, "upazila"), coef=c("incentivized", "cash"), out.intercept="(Intercept)", cluster=c("village", "hhid")) 
original.impact.regress.fun.2 <- regress.fun.factory(depvars=depvars_T1_08, controls="upazila", coef=c("credit", "cash", "info"), cluster=c("village", "hhid")) 

regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.branch.id", "resample.hh.id")) 
village.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.village.id", "resample.hh.id")) 
union.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate"), coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.union.id", "resample.hh.id")) 
impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.branch.id", "resample.hh.id")) 
village.impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.village.id", "resample.hh.id")) 
union.impact.regress.fun <- regress.fun.factory(depvars=depvars_T1_08, coef=c("incentivized", "cash"), controls=c(controls.r2, "upazila"), cluster=c("resample.union.id", "resample.hh.id")) 
# iv.regress.fun <- regress.fun.factory(depvars=depvars_T1_08,  coef="migrant", controls=c(controls.r2, "upazila"), iv="incentivized", cluster=c("resample.branch.id", "resample.hh.id")) 
iv.regress.fun <- regress.fun.factory(depvars="average_calorie_perday2",  coef="migrant_new", controls=c(controls.r2, "upazila"), iv="incentivized", cluster=c("resample.branch.id", "resample.hh.id")) 
village.iv.regress.fun <- regress.fun.factory(depvars="average_calorie_perday2",  coef="migrant_new", controls=c(controls.r2, "upazila"), iv="incentivized", cluster=c("resample.village.id", "resample.hh.id")) 
fake.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate", "average_calorie_perday2"), coef=c("fake.incentivized", "fake.cash"), controls=c(controls.r2, "upazila", "incentivized", "cash"), cluster=c("resample.branch.id", "resample.hh.id")) 
village.fake.regress.fun <- regress.fun.factory(depvars=c("migrant", "ngo.help.migrate", "average_calorie_perday2"), coef=c("fake.incentivized", "fake.cash"), controls=c(controls.r2, "upazila", "incentivized", "cash"), cluster=c("resample.village.id", "resample.hh.id")) 

# round2.data %>% filter(!is.na(village)) %>% original.regress.fun
# round2.data %>% filter(!is.na(village), incentive != "info") %>% original.impact.regress.fun
# round2.data %>% filter(!is.na(village)) %>% original.iv.regress.fun

strat.regress.fun <- regress.fun.factory(depvars=c("average_calorie_perday2", "average_exp2", "migrant"), coef=c("fake.treat", "fake.branch.treat"), controls=c(controls.r2, "resample.branch.id", "cash", "credit", "info"), cluster=c("resample.branch.id", "resample.hh.id")) 
strat.incentive.takeup.regress.fun <- regress.fun.factory(depvars="ngo.help.migrate", coef=c("fake.treat", "fake.branch.treat"), controls=c(controls.r2, "cash", "credit", "info"), cluster=c("resample.branch.id", "resample.hh.id")) 

strat.resampler <- get.stratified.branch.resampler(resample.from.data)

# strat.resample.est.data <- foreach(rep.id=1:2, .combine=bind_rows) %do% {
strat.resample.est.data <- foreach(rep.id=1:10000, .combine=bind_rows) %dopar% {
  resample.data <- strat.resampler(30, num.village.pairs.per.branch=3) #, total.num.villages=180)

  tryCatch({ 
    resample.data %>% 
      # mutate(average_calorie_perday2=ifelse(fake.treat, average_calorie_perday2 + 80, average_calorie_perday2)) %>% 
      strat.regress.fun %>% 
      bind_rows(strat.incentive.takeup.regress.fun(resample.data)) %>% 
      mutate(rep.id=rep.id) 
  }, error=function(e) { 
    cat(sprintf("ERROR: %s\n", e))
    return(NULL) 
  })
}

# union.resample.est.data <- foreach(rep.id=1:1000, .combine=bind_rows) %dopar% {
#   # resample.data <- village.resampler(20, 20, 80) 
#   resample.data <- union.resampler(20, 20, 80) 
#   
# #   village.regress.fun(resample.data) %>% 
# #     bind_rows(village.fake.regress.fun(resample.data)) %>%
# #     village.impact.regress.fun(resample.data) %>%
# #       bind_rows(village.iv.regress.fun(resample.data)) %>% 
#   union.regress.fun(resample.data) %>%
#     bind_rows(union.impact.regress.fun(resample.data)) %>% 
#       mutate(rep.id=rep.id) 
# }

# resample.late.est.data <- foreach(rep.id=1:1000, .combine=bind_rows) %dopar% {
#   # resample.data <- resampler(15, 15, 60) 
#   resample.data <- resampler(37, 31, 32) 
#   
#   iv.regress.fun(resample.data) %>% 
#     bind_rows(fake.regress.fun(resample.data)) %>% 
#     mutate(rep.id=rep.id) 
# }

sim.sum <- . %>% 
  group_by(depvar, coef) %>% 
  summarize_each(funs(mean, sd), est, starts_with("se"))

sim.reject.rate <- . %>% 
  group_by(depvar, coef) %>% 
  summarize_each(funs(mean(. <= 0.05)), ends_with("p.value"))

bs.sim.reject.rate <- . %>% 
  group_by(depvar, coef) %>% 
  summarize(reject.rate=mean((pnorm(abs(est/sd(est)), lower.tail=FALSE) * 2) <= 0.05))
