library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(stringr)
library(foreach)
library(foreign)
library(plm)
library(quantreg)
library(AER)
library(lmtest)
library(car)
library(lme4)
library(doParallel)
library(yaml)

options(contrasts=c("contr.Treatment", getOption("contrasts")[2]))

block.bootstrap.factory <- function(cluster) {
  function(.data) {
    data_frame(original.cluster=as.data.frame(.data)[, cluster] %>% unique %>% sample(length(.), replace=TRUE)) %>%
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
                  select(depvar, est) %>% 
                  rename(origin.est=est),
                .,
                by="depvar") %>%
      transmute(centered.stat=abs((est - origin.est)/se))
  } 
}

regress.fun.factory <- function(depvars, controls, cluster, coef, se.types=c("HC0", "HC1", "HC2", "HC3")) {
  function(.data) {
    inner.reg <- function(depvar) {
      reg.res <- plm(formula(sprintf("%s ~ %s", depvar, paste(c(coef, controls), collapse=" + "))), data=.data, model="pooling", index=cluster) 
      
      # reg.res$coefficients %>%
      reg.res %>%
        coeftest %>%
        magrittr::extract(coef, , drop=FALSE) %>%
        as.data.frame %>%
        set_names(c("est", "se.conv", "t.value.conv", "p.value.conv")) %>%
        mutate(coef=str_extract(rownames(.), "((?<=\\[T\\.)[^\\]]+)|[^\\[\\]]+$")) %>%
        bind_cols(maply(se.types, function(se.type) plm::vcovHC(reg.res, type=se.type) %>% diag %>% sqrt) %>% 
                    t %>% 
                    magrittr::extract(coef, , drop=FALSE) %>%
                    as.data.frame %>%
                    set_names(paste0("se.", se.types))) %>%
        mutate(se=apply(select(., starts_with("se.")), 1, max),
               t.value=est/se,
               p.value=pnorm(abs(t.value), lower.tail=FALSE) * 2)
    }
    
    data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do(inner.reg(.$depvar)) %>%
      ungroup 
  }
}

qr.fun.factory <- function(depvars, controls, cluster, coef, tau=0.5, bootstrap.se=TRUE, num.resample=1000) {
  resampler <- block.bootstrap.factory(cluster)
  
  self.fun <- function(.data, .bootstrap.se=bootstrap.se) tryCatch({ est <- data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do((function(.depvar) tryCatch({
        rq(formula(sprintf("%s ~ %s", .depvar, paste(c(coef, controls), collapse=" + "))), tau=tau, data=.data) %>% 
          coefficients %>% 
          magrittr::extract(coef, , drop=FALSE) %>% 
          as.data.frame %>% 
          mutate(coef=rownames(.)) %>%
          gather(tau, est) %>%
          mutate(tau=sub("tau=\\s+", "", tau) %>% as.numeric)
      }, error=function(err) browser()))(.$depvar)) %>%
      ungroup %>%
      arrange(depvar, tau)
    
    if (.bootstrap.se) {
      se <- foreach(seq_len(num.resample), .combine=cbind) %dopar% {
        # Recursive call to "self" function with resampled data
        resampler(.data) %>% self.fun(FALSE) %>%  select(est)
      } %>% 
        as.matrix %>% 
        aaply(1, sd) # standard deviation of every row (per estimate)
      
      return(est %>% 
               mutate(se=se,
                      t.value=ifelse(!is.na(se), est/se, NA),
                      p.value=ifelse(!is.na(t.value), pnorm(abs(t.value), lower.tail=FALSE) * 2, NA)))
    } else {
      return(est)
    }
  }, error=function(err) browser())
  
  return(self.fun)
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

iv.regress.fun.factory <- function(depvars, controls, cluster, endo.var, iv, bootstrap.se=TRUE, num.resample=1000) {
  resampler <- block.bootstrap.factory(cluster)
  
  self.fun <- function(.data, .bootstrap.se=bootstrap.se) tryCatch({
    est <- data_frame(depvar=depvars) %>% 
      group_by(depvar) %>%
      do(ivreg(formula(sprintf("%s ~ %s | %s", 
                             .$depvar, 
                             paste(c(endo.var, controls), collapse=" + "), 
                             paste(c(iv, controls), collapse=" + "))), data=.data) %>%
           coefficients %>% 
#       do(plm(formula(sprintf("%s ~ %s | . - %s + %s", 
#                              .$depvar, 
#                              paste(c(endo.var, controls), collapse=" + "), 
#                              paste(endo.var, collapse=" + "), 
#                              paste(iv, collapse=" + "))), data=.data, model="pooling", index=cluster) %>%
#            coeftest(vcov=plm::vcovHC(., cluster="group")) %>% 
           # coeftest %>%
           # coeftest(vcov=vcovHC(.)) %>%
           magrittr::extract(endo.var, drop=FALSE) %>% 
           as.data.frame %>% 
           mutate(coef=endo.var)) %>%
      ungroup %>% 
      set_names(c("depvar", "est", "coef")) #, "se", "t.value", "p.value"))
    
    if (.bootstrap.se) {
      se <- foreach(seq_len(num.resample), .combine=cbind) %dopar% {
        # Recursive call to "self" function with resampled data
        resampler(.data) %>% self.fun(FALSE) %>%  select(est)
      } %>% 
        as.matrix %>% 
        aaply(1, sd) # standard deviation of every row (per estimate)
      
      est %<>% 
        mutate(se=se,
               t.value=ifelse(!is.na(se), est/se, NA),
               p.value=ifelse(!is.na(t.value), pnorm(abs(t.value), lower.tail=FALSE) * 2, NA))
    }
    
    return(est)
  }, error=function(err) browser())
  
  return(self.fun) 
}

# Relative likelihood of complier covariate (characteristic)
subpop.covar.rl <- function(.data, cluster, covar, endo.var, iv, bootstrap.se=TRUE, num.resample=1000, subpop=c("compliers", "always.takers", "takers")) {
  resampler <- block.bootstrap.factory(cluster)
  subpop <- match.arg(subpop)
  
  self.fun <- function(.data, bootstrap.se) {
    numerator.prob <- lm(formula(sprintf("%s ~ %s*%s", endo.var, iv, covar)), data=.data) %>% 
      coefficients %>% 
      (function(.coef) {
        magrittr::extract(.coef,
                          switch(subpop,
                                 compliers=c(iv, paste(covar, iv, sep=":")),
                                 always.takers=c("(Intercept)", covar),
                                 takers=c("(Intercept)", iv, paste(covar, iv, sep=":"), covar)))
      }) %>%
      sum(na.rm=TRUE)
    
    denom.prob <- lm(formula(sprintf("%s ~ %s", endo.var, iv)), data=.data) %>%
      coefficients %>% 
      (function(.coef) {
        magrittr::extract(.coef,
                          switch(subpop,
                                 compliers=iv,
                                 always.takers="(Intercept)",
                                 takers=c(iv, "(Intercept)")))
      }) %>%
      sum(na.rm=TRUE)
  
    ratio.est <- numerator.prob/denom.prob  
    
    if (bootstrap.se) {
      se <- foreach(seq_len(num.resample), .combine=c) %dopar% {
        # Recursive call to "self" function with resampled data
        resampler(.data) %>% self.fun(FALSE)
      } %>% sd
      
      data.frame(est=ratio.est) %>%
        mutate(t.value=(est - 1)/se,
               p.value=pnorm(abs(t.value), lower.tail=FALSE) * 2)
    } else {
      return(ratio.est)
    }
  }

  self.fun(.data, bootstrap.se)   
}

always.of.takers.rl <- function(.data, cluster, covar, endo.var, iv, estimator.type=c("diff", "ratio"), num.resample=1000) {
  estimator.type <- match.arg(estimator.type)
  
  if (estimator.type == "diff") {
    .data %>%
      plm(formula(sprintf("%s ~ %s*%s - %s", covar, endo.var, iv, iv)), data=., model="pooling", index=cluster) %>%
      (function(reg.res) {
        # lht.res <- reg.res %>% lht(paste(c(iv, paste(endo.var, iv, sep=":")), sep=" + "), test="F", vcov=plm::vcovHC(., cluster="group"))
        lht.res <- reg.res %>% lht(paste(c(paste(endo.var, iv, sep=":")), sep=" + "), test="F", vcov=plm::vcovHC(., cluster="group"))
        # data.frame(est=sum(reg.res$coefficients[c(iv, paste(endo.var, iv, sep=":"))], na.rm=TRUE),
        data.frame(est=sum(reg.res$coefficients[c(paste(endo.var, iv, sep=":"))], na.rm=TRUE),
                   f.value=lht.res$F[2],
                   p.value=lht.res[2, "Pr(>F)"])
      })
  } else {
    ratio.estimator <- . %>%
      lm(formula(sprintf("%s ~ %s*%s", covar, endo.var, iv)), data=.) %>%
      coefficients %>%
      (function(coef) sum(coef[c("(Intercept)", endo.var)], na.rm=TRUE)/sum(coef[c("(Intercept)", endo.var, iv, paste(endo.var, iv, sep=":"))], na.rm=TRUE))
    
    resampler <- block.bootstrap.factory(cluster)
    ratio.est <- ratio.estimator(.data)
    
    se <- foreach(seq_len(num.resample), .combine=c) %dopar% {
      # Recursive call to "self" function with resampled data
      resampler(.data) %>% ratio.estimator
    } %>% sd
        
    data.frame(est=ratio.est) %>%
      mutate(t.value=(est - 1)/se,
             p.value=pnorm(abs(t.value), lower.tail=FALSE) * 2)
  }
}

predict.rob <- function(x, vcov.=vcov(x), signif.level=0.05, newdata) {
  if (missing(newdata)) { 
    newdata <- x$model 
  }
  
  newdata.col.re <- newdata %>% names %>% paste(collapse="|") %>% sprintf("(^|:)(%s)($|\\[)", .)
  
  m.mat <- x %>% 
    terms %>% 
    (function(trms) {
      newdata.col.re %>% grepl(attr(trms, "term.labels")) %>% not %>% 
        (function(mask) { 
          if (any(mask)) drop.terms(trms, which(mask), keep.response=FALSE) else delete.response(trms) 
        }) 
    }) %>%
    model.matrix(data=newdata)
  
  used.terms <- function(names) {
    ut <- names %>% sub("\\[.+", "", .) %>% grepl(newdata.col.re, .) 
    ut[1] <- TRUE
    
    return(ut)
  }
  
  vcov. %<>% (function(vc) { vc %>% colnames %>% used.terms %>% vc[., .] })
  
  fit <- m.mat %*% x$coef[names(x$coef) %>% used.terms]
  se.fit <- (m.mat %*% vcov. %*% t(m.mat)) %>% diag %>% sqrt

  merr <- qnorm(signif.level/2, lower.tail=FALSE) * se.fit
  
  return(list(fit=fit,
              se.fit=se.fit,
              fit.max=fit + merr,
              fit.min=fit - merr))
}

icc <- function(.data, lhs, cluster, rhs=NULL) {
  formula.str <- sprintf("%s ~ (1|%s)", lhs, cluster)
  
  if (!is.null(rhs)) {
    formula.str %<>% paste(rhs, sep=" + ") 
  }
  
  .data %>%
    lmer(formula(formula.str), data=., REML=FALSE) %>%
    summary %>%
    (function(sum.obj) {
      within.var <- sum.obj$varcor %>% magrittr::extract(cluster) %>% unlist(use.names=FALSE)
      between.var <- sum.obj$sigma^2
      
      c(sqrt(within.var), sum.obj$sigma, within.var/(between.var + within.var)) %>%
        set_names(c("within.sd", "between.sd", "icc"))
    })
}