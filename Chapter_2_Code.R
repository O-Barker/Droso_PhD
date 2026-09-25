#################################### All scripts for Chapter 2 

########################## Mating duration and mating latency
######
library(ggplot2)
library(scales)
library(gridExtra)
library(tidyr)

######

## Survival packages 
library(survival)
library(survminer) # for Kaplan-Meier plots
library(broom) # for tidy output 
library(ggplot2) # (loaded by survminer but I typed it again out of habit)
library(coxme)

## And then my usual friends because coxme will be a bad fit

library(dplyr) # For filtering
library(glmmTMB) # For modelling
library(lme4)
library(DHARMa) # To check model fits

## Post-hoc tests
library(emmeans) # For Tukey tests


HS_master <- read.csv("Heat_shock_master.csv")

head(HS_master)

### Time-to-event like survival

KM_fit <- survfit(Surv(Mating_latency, Mated_binary) ~ Sex + Rearing_diet + Heat_treatment, data = HS_master)
print(KM_fit)

Plot_lat <- ggsurvplot(
  KM_fit, 
  data = HS_master, 
  fun = "event", 
  conf.int = TRUE, 
  legend.title = "Heat treatment",
  xlim = c(0, 120), 
  break.x.by = 10, 
  xlab = "Observation window (minutes)",
  ylim = c(0.00, 1.00), 
  ylab = "Probability of mating",
  facet.by = c("Sex", "Rearing_diet"),
  palette = c("#409bb0", "#e19e6a"),   # only 2 colors now, one per Heat_treatment level
  ggtheme = theme_bw()
)
Plot_lat


### Now the same but mating duration
duration <- filter(HS_master, Mated_binary=="1")

KM_fit <- survfit(Surv(Mating_duration, Mated_binary) ~ Sex + Rearing_diet + Heat_treatment, data = duration)
print(KM_fit)

Plot_dur <- ggsurvplot(
  KM_fit, 
  data = HS_master, 
  fun = "event", 
  conf.int = TRUE, 
  legend.title = "Heat treatment",
  xlim = c(0, 60), 
  break.x.by = 10, 
  xlab = "Time from mating onset (minutes)",
  ylim = c(0.00, 1.00), 
  ylab = "Probability that mating will end",
  facet.by = c("Sex", "Rearing_diet"),
  palette = c("#409bb0", "#e19e6a"),   # only 2 colours now, one per Heat_treatment level
  ggtheme = theme_bw()
)
Plot_dur

#####################################################
#####################################################

## coxme


####### for latency

Cox_fit <- coxme(
  Surv(Mating_latency, Mated_binary) ~ Sex * Heat_treatment * Rearing_diet + (1 | Vial_ID_rand),
  data = HS_master
)

summary(Cox_fit)
cox.zph(Cox_fit)


## 3-way n.s.

Cox_fit1 <- coxme(
  Surv(Mating_latency, Mated_binary) ~ Sex * Heat_treatment + Rearing_diet + (1 | Vial_ID_rand),
  data = HS_master
)

summary(Cox_fit1)
cox.zph(Cox_fit1)
plot(cox.zph(Cox_fit1))
## Looks good, proportional hazards assumption holds!


## Approximate residuals with coxph (no random effect)
Cox_fit_diag <- coxph(
  Surv(Mating_latency, Mated_binary) ~ Sex * Heat_treatment + Rearing_diet,
  data = HS_master
)

# Martingale residuals - check functional form / outliers
mart_res <- residuals(Cox_fit_diag, type = "martingale")
plot(mart_res)

# Deviance residuals - check for outliers
dev_res <- residuals(Cox_fit_diag, type = "deviance")
plot(dev_res)

# dfbeta - check influential observations
dfbeta_res <- residuals(Cox_fit_diag, type = "dfbeta") ### some random points but not many
plot(dfbeta_res)


###### For duration

Cox_fit <- coxme(
  Surv(Mating_duration, Mated_binary) ~ Sex + Heat_treatment + Rearing_diet + (1 | Vial_ID_rand),
  data = HS_master
)

summary(Cox_fit)
cox.zph(Cox_fit)
plot(cox.zph(Cox_fit))


## Approximate residuals with coxph (no random effect)
Cox_fit_diag <- coxph(
  Surv(Mating_latency, Mated_binary) ~ Sex + Heat_treatment + Rearing_diet,
  data = HS_master
)

# Martingale residuals - check functional form / outliers
mart_res <- residuals(Cox_fit_diag, type = "martingale")
plot(mart_res)

# Deviance residuals - check for outliers
dev_res <- residuals(Cox_fit_diag, type = "deviance")
plot(dev_res)

# dfbeta - check influential observations
dfbeta_res <- residuals(Cox_fit_diag, type = "dfbeta") ### some random points but not many
plot(dfbeta_res)

#################################
################################

### Now post-hoc tests

### latency first

# Heat_treatment effect within each Sex
emmeans(Cox_fit1, pairwise ~ Heat_treatment | Sex, adjust = "tukey")

# Sex effect within each Heat_treatment
emmeans(Cox_fit1, pairwise ~ Sex | Heat_treatment, adjust = "tukey")

##################################

## now duration 

#  Sex
emmeans(Cox_fit, pairwise ~ Sex, adjust = "tukey")

# Heat_treatment
emmeans(Cox_fit, pairwise ~ Heat_treatment, adjust = "tukey")

# Diet
emmeans(Cox_fit, pairwise ~ Rearing_diet, adjust = "tukey")


################################################################# Now week 1 offspring
#################################################################

######
library(ggplot2)
library(scales)
library(gridExtra)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
######

HS_master <- read.csv("Heat_shock_master.csv")
head(HS_master)

## Remove individuals which didn't make it to day 7

HS_master <- filter(HS_master, Dead_escaped_sterile!="Died_first_week") ## 264/ 280 left

## Make block a factor

is.factor(HS_master$Block)

HS_master$Block <- as.factor(HS_master$Block)
is.factor(HS_master$Block)

################################################## started with 4-way interaction and simplified 

M11 <- glmmTMB(Week1_offspring ~ Rearing_diet*Heat_treatment*Weight_mg + Heat_treatment*Weight_mg*Sex + Rearing_diet*Sex + Heat_treatment*Rearing_diet*Sex +
                 Block + (1|Vial_ID_rand),
               family = "nbinom1", dispformula = ~1, ziformula = ~1, 
               data = HS_master)
s11 <- simulateResiduals(M11)
plot(s11) ## nice!
summary(M11)
car::Anova(M11)

################################################## Now emmeans and emtrends 

### Because I'm baffled they aren't interacting: 
emmdhs <- emmeans(M11, ~ Rearing_diet * Sex | Heat_treatment)
pairs(emmdhs, adjust = "tukey")
plot(emmdhs)

### Shows expected pattern from graph without sig interaction in model 
### Other view
emmdhs <- emmeans(M11, ~ Rearing_diet * Heat_treatment |Sex )
pairs(emmdhs, adjust = "tukey")
plot(emmdhs)



## Diet, heat, weight

et_diet_heat <- emtrends(M11, ~ Rearing_diet * Heat_treatment, var = "Weight_mg")
et_diet_heat
pairs(et_diet_heat, adjust = "tukey")
plot(et_diet_heat)

## Heat, weight, sex

et_heat_sex <- emtrends(M11, ~ Heat_treatment * Sex, var = "Weight_mg")
et_heat_sex
pairs(et_heat_sex, adjust = "tukey")
plot(et_heat_sex)

##################################################

### Diet and sex only 
## 3 Heat treatment x Sex 
emmdh <- emmeans(M11, ~ Rearing_diet * Sex)
pairs(emmdh, adjust = "tukey")
plot(emmdh)

#################################################

### Heat and sex only 

emmhs <- emmeans(M11, ~ Heat_treatment * Sex)
pairs(emmhs, adjust = "tukey")
plot(emmhs)


################################################

### Weight and sex only emtrends


et_sex_weight <- emtrends(M11, ~ Sex, var = "Weight_mg")
et_sex_weight
pairs(et_sex_weight, adjust = "tukey")
plot(et_sex_weight)

### N.S - like ANOVA

###############################################

### Diet and heat only 

emmdh <- emmeans(M11, ~ Heat_treatment * Rearing_diet)
pairs(emmdh, adjust = "tukey")
plot(emmdh)




#############################################

### Heat only 

emmh <- emmeans(M11, ~ Heat_treatment)
pairs(emmh, adjust = "tukey")
plot(emmh)


###########################################

### Diet only 

emmd <- emmeans(M11, ~ Rearing_diet)
pairs(emmd, adjust = "tukey")
plot(emmd)


######################################################## Now week 2 offspring
########################################################

HS_master <- read.csv("Heat_shock_master.csv")
head(HS_master)

## Remove individuals which didn't make it to day 14 (N = 38)

HS_master <- filter(HS_master, Dead_escaped_sterile!="Died_first_week") ## 264/ 280 left
HS_master <- filter(HS_master, Dead_escaped_sterile!="Died_second_week") ## 242/ 280 left

## Make block a factor

is.factor(HS_master$Block)

HS_master$Block <- as.factor(HS_master$Block)
is.factor(HS_master$Block)

################################################## started with 4-way interaction and simplified
hist(HS_master$Week2_offspring)

########################################## M1 was final model! 
########################################## Removing n.s. 3-way led to worse fit so keeping it in

M1 <- glmmTMB(Week2_offspring ~ Week1_offspring*Rearing_diet*Heat_treatment + Weight_mg*Rearing_diet + Week1_offspring*Weight_mg + Sex*Weight_mg + Sex*Heat_treatment*Rearing_diet + Block + (1| Vial_ID_rand), ziformula=~1, dispformula=~1, family = "compois", data = HS_master)
s1 <- simulateResiduals(M1)
plot(s1) 
testDispersion(s1)
testZeroInflation(s1)
summary(M1)
car::Anova(M1)

### Compois with zi and disp achieved best fit  after trying many variations

######################################################################################### 

## Also tried:

### scaled first week offs and weight to see if this improves fit 

HS_master$Week1_offspring_z <- as.numeric(scale(HS_master$Week1_offspring))
HS_master$Weight_mg_z <- as.numeric(scale(HS_master$Weight_mg))

## Trying nbinom1 again in the first instance

M2 <- glmmTMB(Week2_offspring ~ Week1_offspring_z*Rearing_diet*Heat_treatment 
              + Sex*Rearing_diet + Weight_mg_z*Rearing_diet 
              + Week1_offspring_z*Weight_mg_z + Block + (1|Vial_ID_rand),
              ziformula = ~1, dispformula = ~1, family = "nbinom1", data = HS_master)
s2 <- simulateResiduals(M2)
plot(s2) 
testDispersion(s2)
testZeroInflation(s2)
summary(M2)
car::Anova(M2)

## Slightly worse than compois but not terrible, trying compois + scaled to compare 

M3 <- glmmTMB(Week2_offspring ~ Week1_offspring_z*Rearing_diet*Heat_treatment*Sex 
              + Sex*Rearing_diet + Weight_mg_z*Rearing_diet 
              + Week1_offspring_z*Weight_mg_z + Block + (1|Vial_ID_rand),
              ziformula = ~1, dispformula = ~1, family = "compois", data = HS_master)
s3 <- simulateResiduals(M3)
plot(s3) 
testDispersion(s3)
testZeroInflation(s3)
summary(M3)
car::Anova(M3)

## No 3-way this time! 

##AICs: 

# M3 2130.6
# M2 2168
# M1 2130.6 

## M1 and M3 identical residuals and AIC. Using M1 as final model

##################################################### Post-hoc 
#####################################################

## Diet x weight

et_weight <- emtrends(M1, ~ Rearing_diet, var = "Weight_mg")
et_weight
pairs(et_weight, adjust = "tukey")
plot(et_weight)

## 3-way N.S. but it improved fit so I want to see what's going on with it

et_week1 <- emtrends(M1, ~ Rearing_diet * Heat_treatment | Sex, var = "Week1_offspring")
et_week1
pairs(et_week1, adjust = "tukey")
plot(et_week1)

## not much!

###################################### diet x week 1 only

et_week1 <- emtrends(M1, ~ Rearing_diet, var = "Week1_offspring")
et_week1
pairs(et_week1, adjust = "tukey")
plot(et_week1)

###################################### sex only 
emmd <- emmeans(M1, ~ Sex)
pairs(emmd, adjust = "tukey")
plot(emmd)


###################################### sex x diet (marginal)
emmds <- emmeans(M1, ~ Rearing_diet * Sex)
pairs(emmds, adjust = "tukey")
plot(emmds)
