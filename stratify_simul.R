source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

dgp <- function(num.office, num.vill=2, cluster.size=20, office.sd=0.7, vill.sd=0.2, hh.sd=1) {
  data.frame(office.id=rep(seq_len(num.office), each=num.vill*cluster.size),
             vill.id=rep(seq_len(num.vill*num.office), each=cluster.size),
             hh.id=rep(seq_len(num.office*num.vill*cluster.size)),
             office.shock=rep(rnorm(num.office, sd=office.sd), each=num.vill*cluster.size),
             vill.shock=rep(rnorm(num.vill*num.office, sd=vill.sd), each=cluster.size)) %>% 
    mutate_each(funs(factor), ends_with(".id")) %>% 
    mutate(hh.shock=rnorm(n(), sd=hh.sd),
           y=hh.shock + vill.shock + rnorm(n(), mean=office.shock, sd=0.1))
}

strata.treat <- function(.data, treat.effect=0.1, treat.effect.sd=0.05, cluster.lvl="vill.id", strata="office.id") {
  .data %>% 
    { if (!is.null(strata)) group_by_(., strata) else . } %>% { 
      # Stratified treatment
      inner_join(., 
                 select_(., cluster.lvl) %>% unique %>% mutate(treat=sample(c(rep.int(0, n()/2), rep.int(1, n()/2)))), 
                 by=c(strata, cluster.lvl))
    } %>%
    mutate(strata.treat.effect=rnorm(1, treat.effect, treat.effect.sd) %>% max(0),
           y=y + treat*strata.treat.effect) %>% 
    ungroup
}

vill.treat <- function(.data, treat.effect=0.1, treat.effect.sd=0.05, strata="office.id") {
  .data %>% 
    select(vill.id) %>% 
    distinct %>% 
    mutate(treat=sample(c(rep.int(0, n()/2), rep.int(1, n()/2)))) %>% 
    left_join(.data, ., by="vill.id") %>% 
    { if (!is.null(strata)) group_by_(., strata) else . } %>% 
    mutate(strata.treat.effect=rnorm(1, treat.effect, treat.effect.sd) %>% max(0),
           y=y + treat*strata.treat.effect) %>% 
    ungroup 
}

reg.fun <- . %>% plm(y ~ treat, data=., model="pooling", index="vill.id")
fe.reg.fun <- . %>% plm(y ~ treat + office.id, data=., model="pooling", index="vill.id")

strata.simulation <- function(num.office, num.vill, cluster.size) {
  foreach(rep.id=seq_len(1000), .combine=bind_rows) %dopar% {
    raw.data <- dgp(num.office=num.office, num.vill=num.vill, cluster.size=cluster.size)
    
    treated.strata.data <- raw.data %>% strata.treat(treat.effect=0, treat.effect.sd=0)
    treated.pure.data <- raw.data %>% vill.treat(treat.effect=0, treat.effect.sd=0)
    
    strata.reg.res <- fe.reg.fun(treated.strata.data) # strata dummies
    pure.reg.res <- reg.fun(treated.pure.data)
    
    data.frame(rep.id=rep.id, 
               strata.est=coefficients(strata.reg.res) %>% magrittr::extract("treat"),
               pure.est=coefficients(pure.reg.res) %>% magrittr::extract("treat"),
               strata.conv.se=vcov(strata.reg.res) %>% diag %>% sqrt %>% magrittr::extract("treat"),
               pure.conv.se=vcov(pure.reg.res) %>% diag %>% sqrt %>% magrittr::extract("treat"),
               strata.clust.se=plm::vcovHC(strata.reg.res, cluster="group", type="HC2") %>% diag %>% sqrt %>% magrittr::extract("treat"),
               pure.clust.se=plm::vcovHC(pure.reg.res, cluster="group", type="HC2") %>% diag %>% sqrt %>% magrittr::extract("treat"),
               strata.icc=icc(treated.strata.data, "y", "vill.id", "treat", stratify="office.id") %>% magrittr::extract("icc"),
               strata.2.icc=icc(treated.strata.data, "y", "vill.id", "treat", stratify="office.id", include.strata.var=TRUE) %>% magrittr::extract("icc"),
               pure.icc=icc(treated.pure.data, "y", "vill.id", "treat") %>% magrittr::extract("icc")) %>% 
      mutate(strata.clust.p.value=2 * pnorm(strata.est/strata.clust.se, lower.tail=FALSE),
             pure.clust.p.value=2 * pnorm(pure.est/pure.clust.se, lower.tail=FALSE))
  }
}

impute.icc <- function(clust.se, conv.se, cluster.size) ((clust.se/conv.se)^2 - 1)/(cluster.size - 1)

simul.num.strata <- 40
simul.num.cluster <- 4
simul.cluster.size <- 50

strata.simul.data <- strata.simulation(simul.num.strata, simul.num.cluster, simul.cluster.size)

strat.summ <- strata.simul.data %>% 
  mutate(strata.clust.reject=strata.clust.p.value <= 0.05,
         pure.clust.reject=pure.clust.p.value <= 0.05) %>%
  summarize_each(funs(sd, mean), ends_with("est"), ends_with("icc"), ends_with("se"), ends_with("reject")) %>% 
  mutate(impute.strata.icc=impute.icc(strata.clust.se_mean, strata.conv.se_mean, simul.cluster.size),
         impute.pure.icc=impute.icc(pure.clust.se_mean, pure.conv.se_mean, simul.cluster.size)) %>% 
  as.data.frame