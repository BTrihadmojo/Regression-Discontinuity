##Bambang Trihadmojo 2YP##
##Nov 9, 2021, update##

setwd("")

library(jtools)
library(broom)
library(rdd)
library(rddtools)
library(estimatr)
library(tidyverse)
library(stargazer)
library(psych)
library(apaTables)
library(tableone)
library(Hmisc)
library(corrplot)
library(huxtable)
library(modelsummary)


#Importing variable
d <- read.csv("2YPTrihadmojoData.csv", header = T)

#creating logGNI
d$logGNI <- log(d$GNI)

#checking missing value
is.na(d)

#cleaning up missing value
d <- na.omit(d)

#fuzzy or sharp (Polio) --> fuzzy
ggplot(data = d,
       mapping = aes(x = logGNI,
                     y = Pol3GAVIa,
                     color = Pol3GAVIa)) + 
  geom_vline(xintercept = log(1580), linetype = "dashed", color = "black") +
  geom_point(alpha = 1, position = position_jitter(width = 0, height = 0.3)) +
  labs(x = "Log GNI per Capita", y = "Received GAVI Aid") + 
  ggtitle("GAVI’s 2019 Polio Vaccine Beneficiary Countries by Income") + 
  theme(legend.position = "none")

#identifying crossovers
d %>% 
  group_by(Pol3GAVIa, logGNI <= 7.36518) %>% 
  summarize(count = n()) %>% 
  group_by(Pol3GAVIa) %>% 
  mutate(prop =  count / sum(count)) 

#fuzzy or sharp (PENT) --> fuzzy
ggplot(data = d,
       mapping = aes(x = logGNI,
                     y = PENTGAVIa,
                     color = PENTGAVIa)) + 
  geom_vline(xintercept = log(1580), linetype = "dashed", color = "black") +
  geom_point(alpha = 1, position = position_jitter(width = 0, height = 0.3)) +
  labs(x = "Log GNI per Capita", y = "Received GAVI Aid") +
  ggtitle("GAVI’s 2019 Pentavalent Vaccine Beneficiary Countries by Income") +
  theme(legend.position = "none")

#identifying crossovers
d %>% 
  group_by(PENTGAVIa, logGNI <= 7.36518) %>% 
  summarize(count = n()) %>% 
  group_by(PENTGAVIa) %>% 
  mutate(prop =  count / sum(count)) 

# constructing rdd_data 
polio_rdd <- rdd_data(d$Polio, d$logGNI, cutpoint = log(1580))
PENT_rdd <- rdd_data(d$PENT, d$logGNI, cutpoint = log(1580))

#determining bandwidth for logGNI
polio_bw <- rdd_bw_ik(polio_rdd)
PENT_bw <- rdd_bw_ik(PENT_rdd)

#Fuzzy estimate (parametric lm model)
polio_est <- RDestimate(Polio ~ logGNI + Pol3GAVIa, data = d, 
                        cutpoint = log(1580), bw = 1.4, model = T)

summary(polio_est$model$iv[[2]])
summary(polio_est$model$iv[[1]])
summary(polio_est$model$iv[[3]])

plot(polio_est, gran = log(1580), which = 1)
abline(h = 0, v = 7.36518, col = "red")
title(xlab = "Log GNI per Capita", 
      ylab = "Polio Vaccine Coverage (%)",
      main = "Regression Discontinuity Plot of Polio Vaccine Coverage by GAVI aid")


PENT_est <- RDestimate(PENT ~ logGNI + PENTGAVIa, data = d, 
                       cutpoint = log(1580), bw = 1.54, model = F)

plot(PENT_est, gran = log(1580), which = 1)
abline(h = 0, v = 7.36518, col = "red")
title(xlab = "Log GNI per Capita", 
      ylab = "Pentavalent Vaccine Coverage (%)",
      main = "Regression Discontinuity Plot of Pentavalent Vaccine Coverage by GAVI aid")

summary(PENT_est$model$iv[[2]])
summary(PENT_est$model$iv[[1]])
summary(PENT_est$model$iv[[3]])

#OLS Polio
lmp <- lm(Polio ~ as.factor(Pol3GAVI), data = d)

lmp2 <- lm(Polio ~ as.factor(Pol3GAVI) + logGNI, data = d)

lmp3 <- lm(Polio ~ as.factor(Pol3GAVI) + logGNI + VA + PS + RL + RQ +
             GE + CC, data = d)

lmp4 <- lm(Polio ~ as.factor(Pol3GAVI) + logGNI + VA + PS + RL + RQ +
             GE + CC + Importance + 
             Safety + Effectiveness, data = d)

stargazer(lmp, lmp2, lmp3, lmp4, type = "text", 
         title = "Summary of multiple regression analysis for variables predicting Polio vaccine coverage",
         dep.var.labels = "Polio vaccine coverage", align = T)

#OLS Pentavalent
lmd <- lm(PENT ~ as.factor(PENTGAVI), data = d)

lmd2 <- lm(PENT ~ as.factor(PENTGAVI) + logGNI, data = d)

lmd3 <- lm(PENT ~ as.factor(PENTGAVI) + logGNI + VA + PS + RL + RQ +
             GE + CC , data = d)

lmd4 <- lm(PENT ~ as.factor(PENTGAVI) + logGNI + VA + PS + RL + RQ +
             GE + CC + Importance + 
             Safety + Effectiveness, data = d)

stargazer(lmd, lmd2, lmd3, lmd4, type = "text", 
          title = "Summary of multiple regression analysis for variables predicting Pentavalent vaccine coverage",
          dep.var.labels = "Pentavalent vaccine coverage", align = T)

##################################The End######################################

