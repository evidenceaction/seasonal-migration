source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
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

simul.data <- foreach(rep.id=seq_len(1000), .combine=bind_rows) %dopar% {
  bl.data <- dgp() 
  
  treated.data <- bl.data %>% 
    select(clust.id) %>% 
    distinct %>% 
    sample_frac(0.5) %>% 
    mutate(treat=1) %>% 
    left_join(bl.data, .) %>% 
    mutate(treat=ifelse(is.na(treat), 0, 1))
  
  treated.data %<>% 
    group_by(clust.id) %>% 
    # left_join(select(., individ.id, treat) %>% sample_frac(0.5) %>% mutate(treat.variant=1*treat)) %>% 
    mutate(treat.variant=sample(rep(0:1, n()/2))*treat) %>% 
    ungroup #%>% 
    # mutate(treat.variant=ifelse(is.na(treat.variant), 0, 1))
  
  # reg.res <- plm(y ~ treat + treat.variant, data=treated.data, model="pooling", index="clust.id")
  reg.res <- lm(y ~ treat + treat.variant, data=treated.data)
  
  data.frame(rep.id=rep.id, 
             treat.est=coefficients(reg.res) %>% magrittr::extract("treat"),
             treat.variant.est=coefficients(reg.res) %>% magrittr::extract("treat.variant"),
             treat.variant.se=sandwich::vcovHC(reg.res, type="HC2") %>% diag %>% sqrt %>% magrittr::extract("treat.variant"),
             treat.clust.se=vcov.clx(reg.res, cluster=treated.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat"),
             treat.variant.clust.se=vcov.clx(reg.res, cluster=treated.data$clust.id) %>% diag %>% sqrt %>% magrittr::extract("treat.variant"),
             # treat.conv.se=vcov(reg.res) %>% diag %>% sqrt %>% magrittr::extract("treat"),
             # treat.clust.se=plm::vcovHC(reg.res, cluster="group", type="HC2") %>% diag %>% sqrt %>% magrittr::extract("treat"),
             # treat.variant.clust.se=plm::vcovHC(reg.res, cluster="group", type="HC2") %>% diag %>% sqrt %>% magrittr::extract("treat.variant"),
             treat.icc=icc(treated.data, "y", "clust.id", "treat") %>% magrittr::extract("icc")) %>% 
    mutate(treat.clust.p.value=2 * pnorm(abs(treat.est/treat.clust.se), lower.tail=FALSE),
           treat.variant.p.value=2 * pnorm(abs(treat.variant.est/treat.variant.se), lower.tail=FALSE),
           treat.variant.clust.p.value=2 * pnorm(abs(treat.variant.est/treat.variant.clust.se), lower.tail=FALSE))
}

simul.data %>% summarize_each(funs(sd), ends_with("est")) 
simul.data %>% summarize_each(funs(mean), ends_with("se")) %T>% { print(mutate_each(., funs(. * 2.8))) }
simul.data %>% summarize_each(funs(mean(. <= 0.05)), ends_with("p.value")) 
