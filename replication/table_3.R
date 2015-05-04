source("../util.R")

library(ggplot2)

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

r2.migrants.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/r2_migrants.dta"))
r2.mig.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round2_Migration_byHH.dta"))

round1.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round1_Controls_Table1.dta"))

# round2.sec16b1.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/internal/round2/hh_section16b1.dta"))

round2.sec16.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/internal/round2/hh_section16b.dta")) %>%
  transmute(hhid,
            mid,
            ngo.reason.migrate=q8 == "Help from NGO (credit, money,etc)",
            pay.migrate.1=q14_1,
            pay.migrate.2=q14_2,
            pay.migrate.3=q14_3,
            num.migrate=q3) %>% 
  mutate_each(funs(1*(. %in% c("Help received (Money) from NGO", "NGO credit for migrate"))), pay.migrate.1:pay.migrate.3) %>%
  filter(!is.na(mid)) %>%
  group_by(hhid) %>%
  summarize(ngo.help.migrate=1*(sum(pay.migrate.1, pay.migrate.2, pay.migrate.3, na.rm=TRUE) > 0),
            ngo.reason.migrate=1*any(ngo.reason.migrate, na.rm=TRUE),
            num.migrate=sum(num.migrate, na.rm=TRUE)) %>% 
  ungroup

vill.rdrs.dist.data <- read.csv(paste0(config$data_path, "/mobarak_seasonal_migration/internal/dist_vill_rdrs.csv")) %>%
  transmute(village.id=Village..Number,
            rdrs.dist=Distance) %>%
  mutate(rdrs.dist=sub("km", "", rdrs.dist, ignore.case=TRUE) %>% as.numeric,
         rdrs.dist.q=rdrs.dist %>% cut(quantile(., probs=seq(0, 1, 1/3))),
         rdrs.dist.f=rdrs.dist %>% cut(range(.) %>% multiply_by(c(-1, 1)) %>% sum %>% multiply_by(seq(0, 1, 0.25))),
         rdrs.dist.far=1*(rdrs.dist > 10))

round2.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round2.dta"), convert.factors=FALSE) %>% 
  tbl_df %>% 
  filter(hhid != 92, !is.na(village)) %>% 
  mutate(cluster.id=village,
         incentive2=ifelse(control + info > 0, "control", incentive) %>% factor(levels=c("control", "credit", "cash")),
         groupnature=factor(groupnature) %>% relevel(ref="individual")) %>% 
         # incentivized=1*(incentivized == "Incentivized")) %>%
  # mutate_each(funs(ifelse(is.na(.), NA, ifelse(. == "Yes", 1, 0))), migrant, migrant_new) %>%
  mutate_each(funs(factor), incentive) %>% 
  mutate(incentive=relevel(incentive, ref="control")) %>% 
  left_join(round2.sec16.data, by="hhid") %>% 
  mutate(ngo.help.migrate=1*(ngo.help.migrate & migrant)) %>% 
  left_join(round3.data %>% select(hhid, migrant_r2)) %>%
  left_join(vill.rdrs.dist.data, c("village"="village.id"))

round3.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/Round3.dta")) %>% 
  tbl_df %>%
  filter(average_food3 < 2500) %>%
  mutate(cluster.id=village,
         migrant_r2=ifelse(is.na(migrant_r2), NA, ifelse(migrant_r2 == "Yes", 1, 0))) %>% 
  merge(round1.data, by="hhid", all.x=TRUE, all.y=FALSE, suffixes=c("", ".y")) %>% 
  select(-ends_with(".y")) 

depvars_T1_08 <- c("average_food2", "average_nonfood2", "average_exp2",  "average_calorie_perday2") 
depvars_T1_09 <- c("average_food2", "average_nonfood2", "average_exp2", "average_calorie_perday2", "average_protein_perday2")
appvars_AT4_09 <- c("average_edu_exp_kds2", "average_cloths_shoes2", "average_med_exp_f2", "average_med_exp_m2")
all.depvars <- c(depvars_T1_09, appvars_AT4_09)

controls.r2 <-c("lit", "walls_good", "monga", "dhaka_remit", "dhaka_network", "exp_total_pc_r1", "subsistencer1",
                "num_adltmalesr1", "num_childrenr1", "avgQ13earned", "constrainedr1", "bankedr1")

controls.r3 <- c("litr1", "walls_good", "monga", "dhaka_remit", "dhaka_network", "exp_total_pcr1", "subsistencer1",
              "num_adltmalesr1", "num_childrenr1", "avgQ13earned", "constrainedr1", "bankedr1")

all.treat.1 <- "incentivized"
all.treat.2 <- c("cash", "credit", "info")
all.treat.3 <- c("cash", "credit")

ols.regress.fun.r2 <- regress.fun.factory(depvars=depvars_T1_08, 
                                       controls=c("upazila"), # controls.r2), 
                                       coef=all.treat.3, 
                                       cluster=c("cluster.id", "hhid"))

ols.regress.fun.r3 <- regress.fun.factory(depvars=depvars_T1_09, 
                                       controls=c("upazila", controls.r3), 
                                       coef=all.treat.1, 
                                       cluster="cluster.id")

ngo.help.ols.fun <- regress.fun.factory(depvars=depvars_T1_08, #"ngo.help.migrate", 
                                       controls=c("upazila", controls.r2), 
                                       coef=c("incentivized", "incentivized:cash"), #all.treat.3,
                                       cluster=c("cluster.id", "hhid"))

ngo.help.ols.fun.2 <- regress.fun.factory(depvars="ngo.help.migrate", 
                                       controls=c("upazila", controls.r2), 
                                       coef=c("cash"), #all.treat.3,
                                       cluster=c("cluster.id", "hhid"))

qr.fun.r2 <- qr.fun.factory(depvars=depvars_T1_08,
                         controls=c("upazila", controls.r2), 
                         coef=all.treat.3, 
                         cluster="cluster.id",
                         tau=c(0.1, 0.25, 0.5, 0.75, 0.9),
                         num.resample=10)

qr.fun.r3 <- qr.fun.factory(depvars=depvars_T1_09, 
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
                                         iv=all.treat.1, 
                                         cluster="cluster.id")

iv.regress.fun.r2 <- iv.regress.fun.factory(depvars=depvars_T1_08, 
                                         controls=c("upazila", controls.r2), 
                                         endo.var="migrant_r2",
                                         iv=all.treat.1, 
                                         cluster="cluster.id")

ngo.help.iv.fun <- iv.regress.fun.factory(depvars=depvars_T1_08, 
                                         controls=c("upazila", controls.r2), 
                                         endo.var="ngo.help.migrate",
                                         iv="cash", #all.treat.2, 
                                         cluster="cluster.id")

test.results.r2 <- ols.regress.fun.r2(round2.data) #%>%
  # arrange(desc(abs(t.value)))

# test.results <- ols.regress.fun(round3.data) %>%
#   arrange(desc(abs(t.value)))

ngo.help.results <- round2.data %>% filter(migrant == 1) %>% ngo.help.ols.fun

qr.test.results.r2 <- qr.fun.r2(round2.data) %>%
  arrange(depvar, tau)

# iv.qr.test.results <- iv.qr.fun(round3.data) 
        # 
# first.stage.test.results <- first.stage.regress.fun(round3.data)
# 
# iv.test.results <- iv.regress.fun(round3.data) %>%
#   arrange(desc(abs(t.value)))

ngo.help.iv.results <- round2.data %>% filter(incentivized == 1, migrant == 1) %>% ngo.help.iv.fun

incentive.pred.data <- data.frame(incentive2=levels(round2.data$incentive2) %>% factor(., levels=.))
dist.incentive.pred.data <- expand.grid(cash=0:1, rdrs.dist.far=0:1)

round2.data %>% 
  filter(!is.na(village)) %>% 
  plm(formula(sprintf("migrant ~ incentivized + incentivized:cash + upazila + %s", paste(controls.r2, collapse=" + "))), data=., model="pooling", index=c("village", "hhid")) %>% 
  coeftest(vcov=plm::vcovHC(., cluster="group"))

round2.data %>% 
  filter(!is.na(village), incentivized == 1, migrant == 1) %>% 
  plm(formula(sprintf("ngo.help.migrate ~ cash + upazila + %s", paste(controls.r2, collapse=" + "))), data=., model="pooling", index=c("village", "hhid")) %>% 
  coeftest(vcov=plm::vcovHC(., cluster="group"))

migrate.pred.data <- round2.data %>% 
  filter(!is.na(village)) %>% 
  plm(migrant ~ incentive2, data=., model="pooling", index=c("village", "hhid")) %>% 
  # coeftest(vcov=plm::vcovHC(., cluster="group")))
  predict.rob(signif.level=0.1, vcov=plm::vcovHC(., cluster="group"), newdata=incentive.pred.data) %>% 
  as.data.frame %>% 
  bind_cols(incentive.pred.data) %>% 
  mutate(depvar="migrant")

round2.data %>% 
  filter(!is.na(village)) %>% 
  plm(ngo.reason.migrate ~ incentive2, data=., model="pooling", index=c("village", "hhid")) %>% 
  coeftest(vcov=plm::vcovHC(., cluster="group"))

ngo.help.pred.data <- round2.data %>% 
  filter(!is.na(village)) %>% 
  plm(ngo.help.migrate ~ incentive2, data=., model="pooling", index=c("village", "hhid")) %>% 
  # coeftest(vcov=plm::vcovHC(., cluster="group"))
  predict.rob(signif.level=0.1, vcov=plm::vcovHC(., cluster="group"), newdata=incentive.pred.data) %>% 
  as.data.frame %>% 
  bind_cols(incentive.pred.data) %>% 
  mutate(depvar="ngo.help.migrate")

migrate.dist.pred.data <- round2.data %>% 
  filter(!is.na(village), incentivized == 1) %>% 
  plm(migrant ~ rdrs.dist.far*cash, data=., model="pooling", index=c("village", "hhid")) %>% 
  # coeftest(vcov=plm::vcovHC(., cluster="group")))
  predict.rob(signif.level=0.1, vcov=plm::vcovHC(., cluster="group"), newdata=dist.incentive.pred.data) %>% 
  as.data.frame %>% 
  bind_cols(dist.incentive.pred.data) %>% 
  mutate(depvar="migrant")

ngo.help.dist.pred.data <- round2.data %>% 
  filter(!is.na(village), incentivized == 1) %>% 
  plm(ngo.help.migrate ~ rdrs.dist.far*cash, data=., model="pooling", index=c("village", "hhid")) %>% 
  # coeftest(vcov=plm::vcovHC(., cluster="group"))
  predict.rob(signif.level=0.1, vcov=plm::vcovHC(., cluster="group"), newdata=dist.incentive.pred.data) %>% 
  as.data.frame %>% 
  bind_cols(dist.incentive.pred.data) %>% 
  mutate(depvar="ngo.help.migrate")

theme_set(theme_bw() + 
            theme(panel.border=element_rect(color=NA), 
                  axis.ticks=element_blank(), 
                  strip.background=element_rect(color=NA, size=2)))

migrate.pred.data %>% 
  bind_rows(ngo.help.pred.data) %>% 
  mutate_each(funs(ifelse(is.nan(.), 0, .)), fit.max, fit.min) %>% 
  ggplot(aes(x=incentive2, group=depvar)) + 
  geom_line(aes(y=fit, color=depvar)) + 
  geom_point(aes(y=fit, color=depvar)) +
  geom_ribbon(aes(ymax=fit.max, ymin=fit.min, fill=depvar), alpha=0.2) +
  scale_x_discrete("", labels=c("Control/info", "Credit", "Cash")) +
  scale_y_continuous("Proportion") +
  scale_fill_discrete("", labels=c("Migration", "NGO Assistance Take-up")) +
  scale_color_discrete("", labels=c("Migration", "NGO Assistance Take-up")) +
  theme(legend.position="top")
  
migrate.dist.pred.data %>% 
  bind_rows(ngo.help.dist.pred.data) %>% 
  mutate_each(funs(ifelse(is.nan(.), 0, .)), fit.max, fit.min) %>% 
  mutate_each(funs(factor), rdrs.dist.far, cash) %>% 
  mutate(depvar=factor(depvar, levels=c("migrant", "ngo.help.migrate"), labels=c("Migration", "NGO Assistance Take-up"))) %>%
  ggplot(aes(x=rdrs.dist.far, group=cash)) + 
  geom_line(aes(y=fit, color=cash)) + 
  geom_point(aes(y=fit, color=cash)) +
  geom_ribbon(aes(ymax=fit.max, ymin=fit.min, fill=cash), alpha=0.2) +
  scale_x_discrete("Distance from Nearest RDRS Office", labels=c("Near", "Far (>7 km)")) +
  scale_y_continuous("Proportion") +
  scale_fill_discrete("Incentive", labels=c("Credit", "Cash")) +
  scale_color_discrete("Incentive", labels=c("Credit", "Cash")) +
  facet_wrap(~ depvar) +
  theme(legend.position="top")

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

