library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(foreach)
library(foreign)
library(plm)
library(quantreg)
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
      mutate(bs.cluster.id=seq_len(n())) %>%
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

regress.fun.factory <- function(depvars, controls, cluster, coef, se.types=c("HC0", "HC1", "HC2")) {
  function(.data) {
    inner.reg <- function(depvar) {
      reg.res <- plm(formula(sprintf("%s ~ %s", depvar, paste(c(coef, controls), collapse=" + "))), data=.data, model="pooling", index=cluster) 
      
      reg.res$coefficients %>%
        magrittr::extract(coef) %>%
        as.data.frame %>%
        set_names("est") %>%
        bind_cols(maply(se.types, function(se.type) plm::vcovHC(reg.res, type=se.type) %>% diag %>% sqrt) %>% 
                    t %>% 
                    magrittr::extract(coef, , drop=FALSE) %>%
                    as.data.frame %>%
                    set_names(paste0("se.", se.types)) %>% 
                    mutate(se=apply(., 1, max))) %>%
        mutate(t.value=est/se,
               p.value=pnorm(abs(t.value), lower.tail=FALSE) * 2)
    }
    
    data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do(inner.reg(.$depvar)) %>%
#       do(plm(formula(sprintf("%s ~ %s", .$depvar, paste(c(coef, controls), collapse=" + "))), data=.data, model="pooling", index=cluster) %$%
#            coefficients %>%
#            extract(coef) %>%
#            as.data.frame %>%
#            set_names("est")) %>%
#            coeftest(vcov=plm::vcovHC(., type="HC0", cluster="group")) %>% 
#            extract(coef, , drop=FALSE) %>% 
           # as.data.frame) %>%
      ungroup #%>% 
      # set_names(c("depvar", "est", "se", "t.value", "p.value"))
  }
}

qr.fun.factory <- function(depvars, controls, cluster, coef, tau=0.5) {
  function(.data) tryCatch({
    data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do((function(.depvar) {
        sum.list <- rq(formula(sprintf("%s ~ %s", .depvar, paste(c(coef, controls), collapse=" + "))), tau=tau, data=.data) %>%
           summary(se="boot")
        
        if (length(tau) < 2) {
          sum.list %<>% list
        }
        
        ldply(sum.list, function(x) x$coefficients %>% magrittr::extract(coef, , drop=FALSE)) %>%
          mutate(tau=rep(tau, each=n()/length(tau)))
      })(.$depvar)) %>%
      ungroup %>%
      set_names(c("depvar", "est", "se", "t.value", "p.value", "tau")) %>% 
      arrange(depvar, tau)
  }, error=function(err) browser())
}

# Quantile Treatment Effect
iv.qr.fun.factory <- function(depvars, controls, cluster, endo.var, iv, tau=0.5, num.resample=1000) {
  resampler <- block.bootstrap.factory(cluster)
  
  self.fun <- function(.data, bootstrap.se=TRUE) tryCatch({
    est.data <- data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do((function(.depvar) {
        # Abadie kappa estimate
        probit.z.on.y.x.res <- glm(formula(sprintf("%s ~ (%s)*%s", iv, paste(c(.depvar, controls), collapse=" + "), endo.var)), data=.data, family=binomial(link="probit"))
        
        probit.z.on.x.res <- glm(formula(sprintf("%s ~ %s", iv, paste(controls, collapse=" + "))), data=.data, family=binomial(link="probit"))
        
        .data %<>% 
          merge(probit.z.on.y.x.res$fitted.values %>% as.data.frame %>% set_names("fitted.z.on.y.x"), by="row.names", all.x=TRUE) %>%
          select(-Row.names) %>%
          merge(probit.z.on.x.res$fitted.values %>% as.data.frame %>% set_names("fitted.z.on.x"), by="row.names", all.x=TRUE) %>%
          select(-Row.names) %>%
          filter(!is.na(fitted.z.on.y.x), !is.na(fitted.z.on.x)) %>%
          mutate(fitted.kappa=1 - ifelse(.[, endo.var] == 1, (1 - fitted.z.on.y.x)/(1 - fitted.z.on.x), fitted.z.on.y.x/fitted.z.on.x) %>% pmax(0) %>% pmin(1))
       
        # Weighted quantile regression 
        rq(formula(sprintf("%s ~ %s", .depvar, paste(c(endo.var, controls), collapse=" + "))), tau=tau, data=.data, weights=.data$fitted.kappa) %$%
          coefficients %>% 
          as.matrix %>%
          magrittr::extract(endo.var, , drop=FALSE) %>% 
          as.data.frame %>% 
          set_names(paste0("tau.", sub(".+\\.", "", names(.)))) %>% 
          gather(tau, est, starts_with("tau")) %>% 
          tidyr::extract(tau, into="tau", regex="(\\.\\d+)") %>%
          mutate(tau=as.numeric(tau))
      })(.$depvar)) %>%
      ungroup 
    
    if (bootstrap.se) {
      se <- foreach(seq_len(num.resample), .combine=cbind) %dopar% {
        # Recursive call to "self" function with resampled data
        resampler(.data) %>% self.fun(bootstrap.se=FALSE) %>% arrange(depvar, tau) %>% select(est)
      } %>% 
        as.matrix %>% 
        aaply(1, sd) # standard deviation of every row (per estimate)
     
      est.data %<>% 
        mutate(se=se,
               t.value=est/se,
               p.value=pnorm(abs(t.value), lower.tail=FALSE) * 2)
    }
    
    return(est.data)
  }, error=function(err) browser())
  
  return(self.fun)
}

iv.regress.fun.factory <- function(depvars, controls, cluster, endo.var, iv) {
  function(.data) tryCatch({
    data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do(ivreg(formula(sprintf("%s ~ %s | %s", 
                             .$depvar, 
                             paste(c(endo.var, controls), collapse=" + "), 
                             paste(c(iv, controls), collapse=" + "))), data=.data) %>%
#       do(plm(formula(sprintf("%s ~ %s | . - %s + %s", 
#                              .$depvar, 
#                              paste(c(endo.var, controls), collapse=" + "), 
#                              paste(endo.var, collapse=" + "), 
#                              paste(iv, collapse=" + "))), data=.data, model="pooling", index=cluster) %>%
#            coeftest(vcov=plm::vcovHC(., cluster="group")) %>% 
           # coeftest %>%
           coeftest(vcov=vcovHC(.)) %>%
           magrittr::extract(endo.var, , drop=FALSE) %>% 
           as.data.frame) %>%
      ungroup %>% 
      set_names(c("depvar", "est", "se", "t.value", "p.value"))
  }, error=function(err) browser())
}

# Relative likelihood of complier covariate (characteristic)
complier.covar.rl <- function(.data, cluster, covar, endo.var, iv, num.resample=1000) {
  resampler <- block.bootstrap.factory(cluster)
  
  self.fun <- function(bootstrap.se) {
    numerator.prob <- lm(formula(sprintf("%s ~ %s*%s", endo.var, iv, covar)), data=.data) %>% 
      coefficients %>% 
      magrittr::extract(c(iv, paste(endo.var, iv, sep=":"))) %>% 
      sum(na.rm=TRUE)
    
    first.stage.prob <- lm(formula(sprintf("%s ~ %s", endo.var, iv)), data=.data) %>%
      coefficients %>% 
      magrittr::extract(iv)
  
    ratio.est <- numerator.prob/first.stage.prob  
    
    if (bootstrap.se) {
      se <- foreach(seq_len(num.resample), .combine=c) %dopar% {
        # Recursive call to "self" function with resampled data
        resampler(.data) %>% self.fun(FALSE) 
      } %>% sd
      
      data.frame(est=ratio.est,
                 t.value=(est - 1)/se,
                 p.value=qnorm(abs(t.value), lower.tail=FALSE) * 2)
    } else {
      return(ratio.est)
    }
  }

  self.fun(TRUE)  
}

round3.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round3.dta")) %>% 
  tbl_df %>%
  filter(average_food3 < 2500) %>%
  mutate(cluster.id=village,
         migrant_r2=ifelse(is.na(migrant_r2), NA, ifelse(migrant_r2 == "Yes", 1, 0)))

depvars_T1_09 <- c("average_food2", "average_nonfood2", "average_exp2", "average_calorie_perday2", "average_protein_perday2")
appvars_AT4_09 <- c("average_edu_exp_kds2", "average_cloths_shoes2", "average_med_exp_f2", "average_med_exp_m2")
all.depvars <- c(depvars_T1_09, appvars_AT4_09)

controls.r3 <- c("litr1", "walls_good", "monga", "dhaka_remit", "dhaka_network", "exp_total_pcr1", "subsistencer1",
              "num_adltmalesr1", "num_childrenr1", "avgQ13earned", "constrainedr1", "bankedr1")

all.treat.1 <- "incentivized"
all.treat.2 <- c("cash", "credit", "info")
all.treat.3 <- c("cash", "credit")

ols.regress.fun <- regress.fun.factory(depvars=depvars_T1_09, 
                                       controls=c("upazila", controls.r3), 
                                       coef=all.treat.1, 
                                       cluster="cluster.id")

qr.fun <- qr.fun.factory(depvars=depvars_T1_09, 
                         controls=c("upazila", controls.r3), 
                         coef=all.treat.1, 
                         cluster="cluster.id",
                         tau=c(0.1, 0.25, 0.5, 0.75, 0.9))

iv.qr.fun <- iv.qr.fun.factory(depvars="average_calorie_perday2", #depvars_T1_09, 
                         controls=c("upazila", controls.r3), 
                         endo.var="migrant_r2",
                         iv=all.treat.1,
                         cluster="cluster.id",
                         tau=c(0.1, 0.25, 0.5, 0.75, 0.9))

first.stage.regress.fun <- regress.fun.factory(depvars="migrant_r2", 
                                               controls=c("upazila", controls.r3), 
                                               coef=all.treat.1, 
                                               cluster="cluster.id")

iv.regress.fun <- iv.regress.fun.factory(depvars=depvars_T1_09, 
                                         controls=c("upazila", controls.r3), 
                                         endo.var="migrant_r2",
                                         iv=all.treat.2, 
                                         cluster="cluster.id")

# test.results <- ols.regress.fun(round3.data) %>%
#   arrange(desc(abs(t.value)))

# qr.test.results <- qr.fun(round3.data) %>%
#   arrange(depvar, tau)

iv.qr.test.results <- iv.qr.fun(round3.data) 
        # 
# first.stage.test.results <- first.stage.regress.fun(round3.data)
# 
# iv.test.results <- iv.regress.fun(round3.data) %>%
#   arrange(desc(abs(t.value)))

centered.stats <- bootstrap.c(round3.data, 
                              ols.regress.fun,
                              block.bootstrap.factory("village"),
                              test.results,
                              num.resample=2000) %>% 
  set_names(paste(names(.), seq_len(ncol(.)), sep="."))

alpha <- c(0.15, 0.10, 0.05)

max.w <- foreach(start=seq_len(nrow(centered.stats)), .combine=rbind) %dopar% {
  centered.stats[seq(start, nrow(centered.stats)), ] %>% 
    summarise_each(funs(max)) 
} 

critical.pts <- max.w %>% 
  bind_cols(test.result) %>%
  as.matrix %>% 
  aaply(1, function(row, alpha=c(0.15, 0.10, 0.05)) {
    quantile(row, probs=1 - alpha) %>%
      as.data.frame %>%
      set_names(paste0("critical.c.", 100 * alpha)) %>% 
      mutate(stemp.p.value)
  })

(test.results %<>% bind_cols(critical.pts))
