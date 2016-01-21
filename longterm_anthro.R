anthro.data <- read.dta("~/Downloads/anthropometrics (1).dta") %>%
  mutate(incentive=factor(incentive)) %>%
  rename(sex=s1_1_q3_,
         age=s1_1_q5_1_,
         weight=s1_1_6_,
         height=s1_1_7_,
         muac=s1_1_8_)