
depvar.icc <- function(depvar, include.r3=TRUE, stratify="branch.office", treatments=c("credit", "cash", "info")) {
  rhs <- treatments 
  
  reg.variables <- c(depvar, "village", "incentive", stratify, rhs)
  
  round2.data[, reg.variables] %>% 
    mutate(round=2) %>% { 
      if (include.r3) {
        bind_rows(., round3.data[, reg.variables] %>% 
                    mutate(round=3))
      } else {
        return(.)
      }
    } %>%
    filter(incentive %in% treatments) %>% 
    mutate_each(funs(factor), round) %>%
    icc(depvar, cluster="village", rhs=rhs, stratify=stratify)
}

foreach(depvar=c(depvars_T1_08), .combine=bind_rows) %do% {
  depvar.icc(depvar) %>% t %>% as.data.frame %>% mutate(depvar=depvar)
} %>% 
  bind_rows(depvar.icc("migrant", include.r3=FALSE) %>% t %>% as.data.frame %>% mutate(depvar="migrant")) %>% 
  bind_rows(depvar.icc("ngo.help.migrate", stratify=c("branch.office"), include.r3=FALSE) %>% t %>% as.data.frame %>% mutate(depvar="ngo.help.migrate")) 
