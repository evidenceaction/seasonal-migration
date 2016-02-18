source("util.R")

library(ggplot2)

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores)
}, error=function(err) {
  registerDoSEQ()
})

# DGP factory
logit.dgp.factory <- function(num.treatments = 2,
                              treatment.coefs = 0, # Logit coefficient
                              intercept = -0.69229, # Of logit index
                              branch.re.sd = 0.3041, # Random effects SDs
                              village.re.sd = 0.6323) {
  stopifnot(length(treatment.coefs) == 1) # not supported yet
  
  function(num.branches = 15,
           num.treat.villages.per.branch = 2,
           sample.size.per.village = 50) {
    
    branch.size <- num.treat.villages.per.branch * num.treatments * sample.size.per.village
    
    data.frame(branch.id = rep(seq_len(num.branches), each=branch.size),
               village.id = rep(seq_len(num.branches * num.treat.villages.per.branch * num.treatments), each=sample.size.per.village),
               hh.id = seq_len(num.branches * num.treat.villages.per.branch * num.treatments * sample.size.per.village)) %>% 
      mutate(village.id = paste(branch.id, village.id, sep="-"),
             hh.id = paste(village.id, hh.id, sep="-")) %>% 
      mutate_each(funs(factor), ends_with(".id")) %>% 
      mutate(branch.re = rep(rnorm(num.branches, sd = branch.re.sd), each=branch.size),
             village.re = rep(rnorm(num.branches * num.treat.villages.per.branch * num.treatments, sd = village.re.sd), each=sample.size.per.village)) %>% 
      group_by(branch.id) %>% 
      left_join(., select(., village.id) %>% 
                  distinct(village.id) %>% 
                  mutate(village.treated = sample(rep(0:1, num.treat.villages.per.branch))),
                by=c("branch.id", "village.id"))  %>% 
      left_join(., select(., village.id) %>% 
                  distinct(village.id) %>% 
                  mutate(fake.village.treated = sample(rep(0:1, num.treat.villages.per.branch))),
                by=c("branch.id", "village.id"))  %>% 
      group_by(village.id) %>% 
      mutate(village.model.index = intercept + village.treated*treatment.coefs[1] + branch.re + village.re,
             village.prob = plogis(village.model.index),
             y = rbinom(n(), 1, prob = village.prob),
             village.treatment.effect = treatment.coefs[1] * dlogis(village.model.index)) %>% # I don't remember what this is! 
      ungroup
  }
}

simul.mde <- function(dgp.fun = logit.dgp.factory(), cluster = "village.id", rep=1000, ...) {
  foreach(rep.id=seq_len(rep), .combine=bind_rows) %dopar% tryCatch({
    resampled.data <- dgp.fun(...) 
    
    reg.res.1 <- lm(y ~ village.treated, data=resampled.data)
    
    data.frame(rep.id = rep.id, 
               treat.effect.est = coefficients(reg.res.1) %>% magrittr::extract("village.treated"),
               treat.clust.se = tryCatch(vcov.clx(reg.res.1, cluster = resampled.data[, cluster]) %>% 
                                           { diag(sqrt(.)) } %>% 
                                           magrittr::extract("village.treated"), error = function(e) NA)) %>% 
               # treat.wild.bs.se = wild.bootstrap(y ~ village.treated, 
               #                                   .data = resampled.data, 
               #                                   est.callback = function(reg.res) reg.res$coefficients[1], 
               #                                   cluster=cluster, 
               #                                   bootstrap.rep=50) %>% unlist %>% sd(na.rm = TRUE),
               # treat.icc = icc(resampled.data, "y", cluster, "village.treated") %>% magrittr::extract("icc")) %>% 
      mutate(treat.1.reject = 2 * pnorm(abs(treat.effect.est/treat.clust.se), lower.tail=FALSE)) %>% 
      mutate_each(funs(1*(. <= 0.05)), ends_with("reject"))
  }, error = function(e) return(NULL)) %>%  
    summarize_each(funs(mean = mean(., na.rm=TRUE), sd = sd(., na.rm = TRUE)), ends_with("se"), treat.effect.est, ends_with("reject")) %>% 
    rename(bs.se_mean = treat.effect.est_sd) #%>% 
    # mutate_each(funs(. * 2.8), ends_with("se_mean")) 
}

simul.mde.data <- expand.grid(num.treat.villages.per.branch = 1:4,
                              sample.size.per.village = c(50, 100)) %>% 
  mutate(rep=500) %>% 
  by_row(~ do.call(simul.mde, .)) %>% 
  unnest %>% 
  bind_rows

 plot.simul.mde <- function(.data, y.var, y.lab) {
  .data %>% 
    ggplot(aes_string(x = "num.treat.villages.per.branch", y = y.var, color = "factor(sample.size.per.village)")) +
    geom_line() +
    geom_point() +
    scale_x_continuous("Number of Treatment Villages Per Branch") +
    scale_y_continuous(y.lab) +
    scale_color_discrete("Per Village Sample Size") +
    theme(legend.position = "top")
    # facet_grid(~ num.clust, labeller=function(var, val) sprintf("%d Clusters", val))
 }
 
 simul.mde.data %>% 
   mutate(mde = treat.clust.se_mean * 2.8) %>% 
   plot.simul.mde("mde", "Migration Min. Detectable Effect")
 
theme_set(theme_bw() + 
            theme(panel.border=element_rect(color=NA), 
                  axis.ticks=element_blank(), 
                  strip.background=element_rect(color=NA, size=2)))

simul.mde.data %>% 
  plot.simul.mde("treat.clust.se", "Minimum Detectable Effect") +
  ggtitle("Treatment vs Control Power Analysis") 

simul.mde.data %>% 
  plot.simul.mde("treat.reject", "Probability of Rejecting Size 5% Two-Sided Test") +
  # coord_cartesian(ylim=c(0, 0.10)) +
  ggtitle("Treatment vs Control Power Analysis") 

simul.mde.data %>% 
  plot.simul.mde("treat.variant.clust.se", "Minimum Detectable Effect") +
  ggtitle("Treatment Alternatives Power Analysis")

simul.mde.data %>% 
  plot.simul.mde("treat.variant.reject", "Probability of Rejecting Size 5% Two-Sided Test") +
  coord_cartesian(ylim=c(0, 0.10)) +
  ggtitle("Treatment Alternatives Power Analysis")