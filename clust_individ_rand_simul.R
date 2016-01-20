source("util.R")

library(ggplot2)

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores)
}, error=function(err) {
  registerDoSEQ()
})

dgp <- function(num.clust=100, cluster.size=20, icc=0.1) { # ICC = 0.095 from Ali's study
  data.frame(clust.id=rep(seq_len(num.clust), each=cluster.size),
             individ.id=rep(seq_len(num.clust*cluster.size)),
             clust.shock=rep(rnorm(num.clust, sd=sqrt(icc)), each=cluster.size)) %>% 
    mutate_each(funs(factor), ends_with(".id")) %>% 
    mutate(individ.shock=rnorm(n(), sd=sqrt(1 - icc)),
           y=clust.shock + individ.shock)
}

binary.dgp <- function(num.treat.cluster=50, 
                       num.control.cluster=50, 
                       cluster.size=20, 
                       treat.cluster.size=cluster.size, 
                       control.cluster.size=cluster.size, 
                       treat.effect=0, 
                       treat.variant.effect=0, 
                       baseline.prob=0.01, 
                       icc=0.1) { 
  generated.data <- data.frame(clust.id=c(rep(seq_len(num.treat.cluster), each=treat.cluster.size),
                                          num.treat.cluster + rep(seq_len(num.control.cluster), each=control.cluster.size)),
                               clust.prob=c(rep(rbeta(num.treat.cluster, 
                                                    shape1=baseline.prob * (1 - icc) / icc, 
                                                    shape2=(1 - baseline.prob) * (1 - icc) / icc), 
                                              each=treat.cluster.size),
                                            rep(rbeta(num.control.cluster, 
                                                    shape1=baseline.prob * (1 - icc) / icc, 
                                                    shape2=(1 - baseline.prob) * (1 - icc) / icc), 
                                              each=control.cluster.size))) %>% 
    mutate(individ.id=seq_len(n())) %>% 
    mutate_each(funs(factor), ends_with(".id")) 
  
    generated.data %>% 
      select(clust.id) %>% 
      distinct %>% 
      sample_frac(0.5) %>% 
      mutate(treat=1) %>% 
      left_join(generated.data, .) %>% 
      mutate(treat=ifelse(is.na(treat), 0, 1)) %>% 
      group_by(clust.id) %>% 
      mutate(treat.variant=sample(rep(0:1, n()/2))*treat,
             y=rbinom(n(), 1, prob=min(1, first(clust.prob) + treat.effect*treat + treat.variant.effect*treat.variant))) %>% 
      ungroup %>%  
      mutate(treat.variant.1=treat*(1 - treat.variant),
             treat.variant.2=treat*treat.variant)
}

simul.mde <- function(dgp.fun=binary.dgp, rep=1000, ...) {
  foreach(rep.id=seq_len(rep), .combine=bind_rows) %dopar% tryCatch({
    # treated.data <- bl.data %>% 
    #   select(clust.id) %>% 
    #   distinct %>% 
    #   sample_frac(0.5) %>% 
    #   mutate(treat=1) %>% 
    #   left_join(bl.data, .) %>% 
    #   mutate(treat=ifelse(is.na(treat), 0, 1),
    #          y=y + treat*treat.effect)
    # 
    # treated.data %<>% 
    #   group_by(clust.id) %>% 
    #   mutate(treat.variant=sample(rep(0:1, n()/2))*treat,
    #          y=y + treat.variant*treat.variant.effect) %>% 
    #   ungroup 
    
    resampled.data <- dgp.fun(...) 
    
    reg.res.1 <- lm(y ~ treat, data=resampled.data)
    reg.res.2 <- lm(y ~ treat + treat.variant.2, data=resampled.data)
    
    data.frame(rep.id=rep.id, 
               treat.1.est=coefficients(reg.res.1) %>% magrittr::extract("treat"),
               treat.2.est=coefficients(reg.res.2) %>% magrittr::extract("treat"),
               # treat.variant.1.est=coefficients(reg.res.2) %>% magrittr::extract("treat.variant.1"),
               treat.variant.2.est=coefficients(reg.res.2) %>% magrittr::extract("treat.variant.2"),
               treat.1.clust.se=tryCatch(vcov.clx(reg.res.1, cluster=resampled.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat"), error=function(e) NA),
               treat.2.clust.se=tryCatch(vcov.clx(reg.res.2, cluster=resampled.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat"), error=function(e) NA),
               # treat.variant.1.clust.se=tryCatch(vcov.clx(reg.res.2, cluster=resampled.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat.variant.1"), error=function(e) NA), 
               treat.variant.2.clust.se=tryCatch(vcov.clx(reg.res.2, cluster=resampled.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat.variant.2"), error=function(e) NA), 
               treat.icc=icc(resampled.data, "y", "clust.id", "treat") %>% magrittr::extract("icc")) %>% 
      mutate(treat.1.reject=2 * pnorm(abs(treat.1.est/treat.1.clust.se), lower.tail=FALSE),
             treat.2.reject=2 * pnorm(abs(treat.2.est/treat.2.clust.se), lower.tail=FALSE),
             # treat.variant.1.reject=2 * pnorm(abs(treat.1.variant.est/treat.variant.1.clust.se), lower.tail=FALSE),
             treat.variant.2.reject=2 * pnorm(abs(treat.variant.2.est/treat.variant.2.clust.se), lower.tail=FALSE)) %>% 
      mutate_each(funs(1*(. <= 0.05)), ends_with("reject"))
  }, error=function(e) return(NULL)) %>%  
    summarize_each(funs(mean(., na.rm=TRUE)), ends_with("se"), treat.icc, ends_with("reject")) %>% 
    mutate_each(funs(. * 2.8), ends_with("se")) 
}

simul.mde.data <- expand.grid(num.clust=seq(10, 50, by=5), cluster.size=seq(20, 60, by=20), treat.clust.size.incr=seq(0, 60, by=20)) %>% 
  mutate(mde=map3(num.clust, cluster.size, treat.clust.size.incr, ~ simul.mde(dgp.fun=function(...) binary.dgp(..., treat.effect=0.03), rep=1000, num.treat.cluster=.x, num.control.cluster=.x, treat.cluster.size=.y + .z, control.cluster.size=.y))) %>% 
  unnest %>% 
  mutate(treat.cluster.size=cluster.size + treat.clust.size.incr,
         num.clust=num.clust*2) %>% 
  arrange(num.clust, treat.cluster.size)  

 plot.simul.mde <- function(.data, y.var, y.lab) {
  .data %>% 
    ggplot(aes_string(x="cluster.size", y=y.var, color="factor(treat.clust.size.incr)")) +
    geom_line() +
    geom_point() +
    scale_x_continuous("Cluster Size") +
    scale_y_continuous(y.lab) +
    scale_color_discrete("Treatment Cluster Size Increase") +
    facet_grid(~ num.clust, labeller=function(var, val) sprintf("%d Clusters", val))
}
 
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