source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

get.resampler <- function(.data, branch="RDRS.Office.Location.Name", num.villages.per.branch=1) {
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
        ungroup %>%
        mutate(resample.hh.id=seq_len(n())) %>%
        mutate_each(funs(paste(., paste(incentive, collapse="-"), sep="-")), matches("resample\\..+\\.id"))
    }
    
    incentive.resampler("cash", num.cash, cluster.size) %>% 
      bind_rows(incentive.resampler("credit", num.loan, cluster.size)) %>% 
      bind_rows(incentive.resampler(c("control", "info"), num.control, cluster.size)) %>% # Treating "info" treatment as control
      inner_join(select(., resample.branch.id) %>% unique %>% mutate(fake.incentivized=sample(c(rep.int(1, num.cash + num.loan), rep(0, num.control)))), by="resample.branch.id")
  }
  
  return(inner.resampler)
}

resample.from.data <- round2.data %>%
  filter(!is.na(RDRS.Office.Location.Name), !is.na(village))

resampler <- get.resampler(resample.from.data)

original.regress.fun <- regress.fun.factory(depvars="migrant", controls=NULL, cluster=c("village", "hhid"), coef=c("incentivized", "cash"))
# regress.fun <- regress.fun.factory(depvars="migrant", controls=NULL, cluster=c("resample.branch.id", "resample.hh.id"), coef=c("cash", "credit", "fake.incentivized")) #coef="incentivized")
# regress.fun <- regress.fun.factory(depvars="migrant", controls=NULL, cluster=c("resample.branch.id", "resample.hh.id"), coef=c("cash", "info", "fake.incentivized")) #coef="incentivized")
# regress.fun <- regress.fun.factory(depvars="migrant", controls=NULL, cluster=c("resample.branch.id", "resample.hh.id"), coef=c("cash")) 
regress.fun <- regress.fun.factory(depvars="migrant", controls=NULL, cluster=c("resample.branch.id", "resample.hh.id"), coef=c("incentivized", "cash")) 

round2.data %>% filter(!is.na(village)) %>% original.regress.fun

resample.est.data <- foreach(rep.id=1:1000, .combine=bind_rows, .verbose=TRUE) %dopar% {
  resampler(15, 15, 30, 20) %>% 
    regress.fun %>% 
    mutate(rep.id=rep.id) 
}
