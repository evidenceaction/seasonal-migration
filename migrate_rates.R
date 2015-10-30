source("util.R")

library(plm)
library(dplyr)
library(tidyr)
library(yaml)
library(ggplot2)
library(scales)

config <- yaml.load_file("local_config.yaml")

theme_set(theme_bw() + 
            theme(panel.border=element_rect(color=NA), 
                  axis.ticks=element_blank(), 
                  strip.background=element_rect(color=NA, size=2)))

dodge_obj <- position_dodge(0.9)

# read.csv(paste0(config$data_path, "/mobarak_seasonal_migration/internal/Migration_Rates.csv")) %>% 
#   mutate(id=rownames(.)) %>% 
#   gather(key, value, -id) %>% 
#   tidyr::extract(key, c("question", "year"), "([^\\.]+)\\.(\\d+)") %>% 
  # spread(question, value) %>% 
migrate.rates <- read.csv(paste0(config$data_path, "/mobarak_seasonal_migration/internal/Migration_Rates_2.csv")) %>% 
  gather(key, value, -hh.id, -village.id, -incentivized.08) %>% 
  tidyr::extract(key, c("question", "year"), "([^\\.]+)\\.(\\d+)") %>% 
  spread(question, value) %>% 
  mutate_each(funs(factor), year, hh.id, village.id, incentivized.08) %>% 
  filter(!is.na(village.id)) %>% 
  plm(migrant ~ year*incentivized.08, data=., model="pooling", index="village.id") %>% 
  { bind_cols(predict.rob(., signif.level=0.1, vcov=plm::vcovHC(., type="HC2", cluster="group")) %>% as.data.frame,
              .$model) %>% 
      select(-migrant) 
  } %>% 
  distinct %>% 
  arrange(year, desc(incentivized.08))

migrate.rates %>%
  ggplot(aes(x=year, y=fit, group=incentivized.08)) +
  geom_bar(aes(fill=incentivized.08), stat="identity", alpha=0.5, position=dodge_obj, width=0.9) +
  geom_errorbar(aes(ymax=fit.max, ymin=fit.min), position=dodge_obj, width=0.2) +
  geom_text(aes(y=fit + 0.015, label=sprintf("%.0f%%", round(fit*100))), size=3, hjust=-0.15, position=dodge_obj) +
  scale_x_discrete("", labels=c("2008", "2009", "2011")) +
  scale_y_continuous("Migration Rate", labels=percent) +
  scale_fill_discrete("", labels=c("Non-Incentivized", "Incentivized in 2008")) +
  theme(legend.position="bottom")
