########## Analyses of adult diet choice, oviposition choice, larval choice and survival
########## Compiled into single script as ESM due to issues with repository server 


############################## Diet choice analysis: negative binomial GLMM

library(glmmTMB)
library(tidyr)
library(DHARMa)
library(hablar)
library(dplyr)

Diet <- read.csv("Diet_choice.csv") ## In raw data file
head(Diet)

hist(Diet$Flies)
### zeroinflated

var(Diet$Flies)/mean(Diet$Flies) ## Poisson unsuitable

##################### Model selection using dHARMA and ANOVA
####################

m0 <- glmmTMB(Flies ~ Sex*Age*Flies_choice + (1|Time) + (1|Plate) + (1|Obs), zi=~1, family = 'nbinom1', data = Diet)
s0 <- simulateResiduals(m0)
plot(s0)
testZeroInflation(s0) ##zeroinflated
testDispersion(s0) ## overdispersed

# 

m1 <- glmmTMB(Flies ~ Sex*Age*Flies_choice + (1|Time) + (1|Plate) + (1|Obs), zi=~1, family = 'truncated_nbinom1', data = Diet)
s1 <- simulateResiduals(m1)
plot(s1)
testZeroInflation(s1) ## 
testDispersion(s1) ## 

anova(m0, m1)

## m0 is better, stick with negative binomial 

m0 <- glmmTMB(Flies ~ Sex*Age*Flies_choice + (1|Time) + (1|Plate) + (1|Obs), zi=~1, family = 'nbinom1', data = Diet)
s0 <- simulateResiduals(m0)
plot(s0)
testZeroInflation(s0) 
testDispersion(s0) 

m1 <- glmmTMB(Flies ~ Sex*Age*Flies_choice + (1|Time) + (1|Plate), zi=~1, family = 'nbinom1', data = Diet)
s1 <- simulateResiduals(m1)
plot(s1)
testOutliers(s1, type = 'bootstrap') 
testZeroInflation(s1)
testDispersion(s1)  

anova(m0, m1)
## Removing OLRE reduced AIC, but keeping RE for observation intervals

m2 <- glmmTMB(Flies ~ Sex*Age + Age*Flies_choice + Sex*Flies_choice + (1|Time) + (1|Plate), zi=~1, family = 'nbinom1', data = Diet)
s2 <- simulateResiduals(m2)
plot(s2)
 

anova(m0, m1, m2)
## m1 is still better 

########################### Final model is M1

summary(m1) 


################################################# Plot

Plot1 <- ggplot(Diet, aes(
  x = factor(Flies_choice, level = M.order),
  y = Flies,
  colour = factor(Age, level = Age.order)
)) +
  geom_errorbar(
    data = emm_df,
    aes(
      y    = response,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      group = factor(Age, level = Age.order)
    ),
    width = 0.4,
    linewidth = 1.2,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    data = emm_df,
    aes(
      y     = response,
      group = factor(Age, level = Age.order)
    ),
    size = 3,
    position = position_dodge(width = 0.75)
  ) +
  ylab("Flies per patch per recording") +
  xlab("P:C content of diet media") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = pretty_breaks(),
    oob = scales::squish
  ) +
  scale_colour_manual(values = c("#78dede", "#51A4Ad", "#4588e0")) +
  geom_jitter(outlier.shape = NA, 
              position = position_jitterdodge(
                jitter.width = 0.2,
                jitter.height = 0.4,
                dodge.width = 0.75
              ), alpha = 0.15,
              size = 1.5
  ) +
  facet_wrap(~Sex) + labs(color = "Age")
Plot1

################################################################################# End of diet choice script
#################################################################################
#################################################################################


################################# Start of oviposition choice analysis script: betabinomial GLMM

library(glmmTMB)
library(tidyr)
library(DHARMa)
library(hablar)
library(dplyr)

Ovi <- read.csv("Age_eggs.csv")
head(Ovi)

Ovi <- Ovi %>% mutate(Block=factor(Block))
Ovi <- Ovi %>% mutate(Age=factor(Age))
Ovi <- Ovi %>% mutate(Replicate=factor(Replicate))
Ovi <- Ovi %>% mutate(Obs=factor(Obs))
Ovi <- Ovi %>% mutate(Diet=factor(Diet))



# Assuming `total_eggs` is the total number of trials (e.g., total eggs laid)
Ovi$successes <- Ovi$Eggs
Ovi$failures <- Ovi$Total_eggs - Ovi$successes

# Create a two-column matrix for the response variable
Ovi$response_matrix <- cbind(Ovi$successes, Ovi$failures)

# Fit the model with the beta-binomial family
m1 <- glmmTMB(response_matrix ~ Age*Diet + (1|Replicate) + (1|Obs), 
              family = betabinomial, 
              data = Ovi)
s1 <- simulateResiduals(m1)
plot(s1) 
testDispersion(s1) ## Looks good

summary(m1)
car::Anova(m1)


m0 <- glmmTMB(response_matrix ~ Age*Diet + (1|Replicate) + (1|Obs), 
              family = betabinomial, 
              data = Ovi)
s0 <- simulateResiduals(m0)
plot(s0) ## 
testDispersion(s0) ## 

anova(m0, m1)

## both the same 

#write.csv(Ovi, "C:...")


########### emmeans

Emm_ovi <- emmeans(m1, "Age", "Diet", type = "response")
pairs(Emm_ovi)
emm_df <- as.data.frame(Emm_ovi)

###############

#### Orders 
Age.order <- c('5 days', '15 days', '25 days')

M.order <- c('1.4', '1.1', '4.1')



####### Raw egg counts

Plot <- ggplot(Ovi, aes(
  x = factor(Diet, level = M.order),
  y = Eggs,
  colour = factor(Age, level = Age.order)
)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Raw counts of eggs laid per diet patch by groups of 10 females") +
  xlab("P:C content of diet media") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 400),
    breaks = pretty_breaks(),
    oob = scales::squish
  ) +
  scale_colour_manual(values = c("#7D94B3", "#51a4ad", "#4588e0")) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.10,# Horizontal jitter
      jitter.height = 0.4, # vertical jitter
      dodge.width = 0.75 # Dodge by Age
    ), alpha = 0.2,
    size = 1.5 # Adjust point size for clarity
  ) + labs(color = "Age")

Plot


### 


#### Now errbr

Plot <- ggplot(Ovi, aes(
  x = factor(Diet, level = M.order),
  y = Eggs,
  colour = factor(Age, level = Age.order)
)) +  geom_errorbar(
  data = emm_df,
  aes(
    y    = prob,
    ymin = asymp.LCL,
    ymax = asymp.UCL,
    group = factor(Age, level = Age.order)
  ),
  width = 0.2,
  linewidth = 1.2,
  position = position_dodge(width = 0.75)
) +
  geom_point(
    data = emm_df,
    aes(
      y     = prob,
      group = factor(Age, level = Age.order)
    ),
    size = 3,
    position = position_dodge(width = 0.75)
  )  +
  ylab("Probability of oviposition site choice by a female") +
  xlab("P:C content of diet media") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = pretty_breaks(),
    oob = scales::squish
  ) +
  scale_colour_manual(values = c("#7D94B3", "#51a4ad", "#4588e0")) +
  labs(color = "Age")

Plot

### Model used binary 0 - 1 for egg laid vs not laid so can't use jitter to show data points 

################################################################################# End of oviposition choice script
#################################################################################
#################################################################################


################################# Start of larval diet choice analysis script: Poisson GLMM




library(glmmTMB)
library(tidyr)
library(DHARMa)
library(hablar)
library(dplyr)


Larv <- read.csv("Larval_choice.csv")
head(Larv)


Larv <- Larv %>% mutate(Parental_age = factor(Parental_age)) 
Larv <- Larv %>% mutate(Diet_pcb = factor(Diet_pcb))
Larv <- Larv %>% mutate(Replicate = factor(Replicate))
Larv <- Larv %>% mutate(Obs = factor(Obs))

hist(Larv$Larvae)



m0 <- glmmTMB(Larvae ~ Diet_pcb*Parental_age + (1|Replicate) + (1|Obs), data = Larv, family = 'poisson' )

s0 <- simulateResiduals(m0)
plot(s0) 


testOverdispersion(s0) 
testZeroInflation(s0)  

## no need to fix what isn't broken! 

car::Anova(m0)

summary(m0)

##################### Plot

## Make relevant things a factor 

Larv <- Larv %>% mutate(Parental_age = factor(Parental_age))
Larv <- Larv %>% mutate(Diet_pcb = factor(Diet_pcb))

## Orders

Age.order <- c('5 days', '15 days', '25 days')

M.order <- c('1.4', '1.1', '4.1')


## Best model was just an easy poisson

m0 <- glmmTMB(Larvae ~ Diet_pcb*Parental_age + (1|Replicate) + (1|Obs), data = Larv, family = 'poisson' )

emms <- emmeans(m0, "Parental_age", "Diet_pcb", type = "response")
emm_df <- as.data.frame(emms)


Plot <- ggplot(Larv, aes(
  x = factor(Diet_pcb, level = M.order),
  y = Larvae,
  colour = factor(Parental_age, level = Age.order)
)) +
  geom_errorbar(
    data = emm_df,
    aes(
      y    = rate,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      group = factor(Parental_age, level = Age.order)
    ),
    width = 0.2,
    linewidth = 1.2,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    data = emm_df,
    aes(
      y     = rate,
      group = factor(Parental_age, level = Age.order)
    ),
    size = 3,
    position = position_dodge(width = 0.75)
  ) +
  ylab("L2 larvae recovered from each patch after 48hrs") +
  xlab("P:C content of diet media") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = pretty_breaks(),
    oob = scales::squish
  ) +
  scale_colour_manual(values = c("#7D94B3", "#51a4ad", "#4588e0")) +
  geom_jitter(outlier.shape = NA, 
              position = position_jitterdodge(
                jitter.width = 0.2,
                jitter.height = 0.4,
                dodge.width = 0.75
              ), alpha = 0.2,
              size = 1.5
  ) +
  labs(color = "Parental age")
Plot

################################################################################# End of larval diet choice script
#################################################################################
#################################################################################


################################# Start of survival and development analyses: event history 

## Survival packages 
library(survival)
library(survminer) # for Kaplan-Meier plots
library(broom) # for tidy output 
library(ggplot2) 
library(coxme)

## And then my usual friends because coxme was a poor fit

library(dplyr) # For filtering
library(glmmTMB) # For modelling
library(lme4)
library(DHARMa) # To check model fits

## Post-hoc tests
library(emmeans) # For Tukey tests

#################### Load data

Data <- read.csv("Reformat_long_cox.csv")
head(Data)

################### Kaplan-Meier plot

## Kaplan-Meier
KM_fit <- survfit(Surv(Stop, Event) ~ Parent_age + Medium, data = Data)
print(KM_fit)

## Basic KM plot, fun = "event" becomes time to event rather than the default, which is time to death (which in this case would show probability of NOT eclosing rather than probability of eclosing)
Plot1 <- ggsurvplot(KM_fit, data = Data, pval = TRUE, fun = "event")$plot
Plot1

## initial plot to see what data looks like

Plot4 <- ggsurvplot(KM_fit, data = Data, fun = "event", conf.int = TRUE, legend.title = "Rearing diet", xlim = c(8,14), break.x.by = 1, xlab = "Age of offspring (days)", ylim = c(0.00, 1.00), ylab = c("Probability of eclosion"), facet.by = "Parent_age", palette = c("#339966", "#2E9FDF", "#9966CC"), ggtheme = theme_minimal()) 

Plot4

################### Cox proportional hazards model

survdiff(Surv(Stop, Event) ~ Parent_age_num + Medium, data=Data)
### Checks for difference in survival curves between treatments

Cox1 <- coxph(Surv(Stop, Event) ~ Parent_age_num*Medium, data=Data)
### Cox proportional hazards model (not mixed so no random effects)

summary(Cox1)
car::Anova(Cox1)

### Test model fit
cox.zph(Cox1)

### Global P is very significant, which indicates violation of the proportional hazards assumption - the hazard ratios for medium and/ or parental age do not remain constant over time

### The below is a mixed cox proportional hazards model - tried out of interest to see if it was any different, but still violated the proportional hazards assumption. It also crashed R nearly every time I tried to run it. 

cox2 <- coxme(Surv(Stop, Event) ~ Parent_age_num*Medium + (1|Vial) + (1|Individual_ID), data=Data)

summary(cox2)

cox.zph(cox2) ## Global P still very significant so an event history analysis is needed

# Need to reformat the data again into a different kind of long format, so every timepoint has a 0 or a 1 for eclosion
Data_long <- Data %>%
  rowwise() %>%
  do({
    row <- .
    data.frame(
      Individual_ID = row$Individual_ID,
      Vial = row$Vial,
      Parent_age = row$Parent_age,
      Medium = row$Medium,
      time = 1:row$Stop,
      event = c(rep(0, row$Stop - 1), row$Event)
    )
  }) %>%
  ungroup()

## Saving for quick loading next time in case that's useful 
write.csv(Data_long, "Data_long.csv", row.names = FALSE)

## Making everything relevant a factor
Data_long$Medium <- factor(Data_long$Medium)
Data_long$Parent_age <- factor(Data_long$Parent_age)
Data_long$Individual_ID <- factor(Data_long$Individual_ID)
Data_long$Vial <- factor(Data_long$Vial)

########################################################### Modelling probability of eclosing 


## Simplified data (easier to see mess-ups)
event_summary <- Data_long %>%
  group_by(Individual_ID, Vial, Parent_age, Medium) %>%
  summarise(event_occurred = max(event), .groups = "drop")

## Fitting binomial glmm
glm_event <- glmmTMB(event_occurred ~ Parent_age * Medium + (1 | Vial),
                     family = binomial, data = event_summary)
summary(glm_event)

car::Anova(glm_event)

## Checking model fit

glm_event_resid <- simulateResiduals(glm_event)
plot(glm_event_resid) 

### Post-hoc Tukey tests

## Medium by parent age
emm_eclosion_prob <- emmeans(glm_event, "Parent_age", "Medium", data = event_summary)
pairs(emm_eclosion_prob, adjust = "tukey")
plot(emm_eclosion_prob)

## Parent age by medium
emm_eclosion_prob2 <- emmeans(glm_event, "Medium", "Parent_age",  data = event_summary)
pairs(emm_eclosion_prob2, adjust = "tukey")
plot(emm_eclosion_prob2)
############################################################ Now modelling time to eclosion 

event_times <- Data %>%
  filter(Event == 1)  # only those who eclosed are relevant for this one

hist(event_times$Stop) # very oddly dispersed: lots of faff was required for this step to find the right fit to fix the dispersion

event_times <- subset(event_times, Medium != "P") ## needed to drop P entirely because there were so few observations in the old-parent P group

#### Additional trial and error: 

#glmm_time <- glmer(Stop ~ Parent_age * Medium + (1|Individual_ID) + (1|Vial), family = poisson, data = event_times)
#summary(glmm_time)

## Singular fit for both vial and individual
## Removing these so it's now a glm not a glmm
#glm_time <- glm(Stop ~ Parent_age * Medium, family = poisson, data = event_times)
#summary(glm_time)

#glm_sim <- simulateResiduals(glm_time)
#plot(glm_sim)
#testDispersion(glm_sim)

## This fit was still underdispersed: finally, tried conway-maxwell poisson family 
## Added a random effect (vial) back in as it's a glmm

glmm_time <- glmmTMB(Stop ~ Parent_age * Medium + (1|Vial), family = compois, data = event_times)
summary(glmm_time)

car::Anova(glmm_time)

glmm_sim <- simulateResiduals(glmm_time)
plot(glmm_sim)
testDispersion(glmm_sim)

#### Post-hoc Tukey tests using emmeans

## Medium by parent age
emm_eclosion_time <- emmeans(glmm_time, "Parent_age", "Medium", data = event_times)
pairs(emm_eclosion_time, adjust = "tukey")
plot(emm_eclosion_time)

## Parent age by medium
emm_eclosion_time2 <- emmeans(glmm_time, "Medium", "Parent_age",  data = event_times)
pairs(emm_eclosion_time2, adjust = "tukey")
plot(emm_eclosion_time2)


##################################### final KM plot

head(Data)

Data <- Data %>% mutate(
  Medium = recode(Medium,
                  "B" = "1:1",
                  "P" = "4:1",
                  "C" = "1:4"),
  Parent_age = recode(Parent_age,
                      "Young" = "5 days old",
                      "Midlife" = "15 days old" ) )

Data <- Data %>%
  mutate(Medium = fct_relevel(Medium, "1:4", "1:1", "4:1"))

Data <- Data %>%
  mutate(Parent_age = fct_relevel(Parent_age, "5 days old", "15 days old"))

fit <- survfit(Surv(Stop, Event) ~ Parent_age + Medium, data = Data)


Plot4 <- ggsurvplot(fit, data = Data, fun = "event", conf.int = TRUE, legend.title = "Larval diet", xlim = c(8,14), break.x.by = 1, xlab = "Age of offspring (days)", ylim = c(0.00, 1.00), ylab = c("Probability of eclosion"), facet.by = "Parent_age", palette = c("thistle3", "wheat4", "rosybrown3"), ggtheme = theme_bw()) 
Plot4

############################################# Wing / thorax proportionality

## New analysis for proportionality of wing/ thorax for ch3 dataset 08/07/2025

library(glmmTMB)
library(dplyr)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(scales)

setwd("C:/Users/barke/OneDrive - University of East Anglia/Experiments/2024/Age effects B2/No_choice_devo/Wing_thorax_stuff")

Wingth <- read.csv("Individuals.csv") 
head(Wingth)

############## plot first 

Age.order <- c('Young', 'Midlife')

M.order <- c('C', 'B', 'P')

## 
Plot <- ggplot(Wingth, aes(x=factor(Medium, level = M.order), y = Wing_over_thorax, colour = factor(Parental_age, level = Age.order))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha=0.3, position = position_jitterdodge()) + ylab("Wing length / thorax length (mm)") + xlab("P:C content of diet media") + theme_minimal() + scale_y_continuous(limits = c(1, 4), breaks = pretty_breaks(), oob = scales::squish) + scale_colour_manual(values = c( "#5e76b3","#9e5eb3")) + facet_wrap(~Sex) + theme_bw()
Plot

####### Now model 

hist(Wingth$Wing_over_thorax) ## Probably Gaussian

m0 <- glmmTMB(Wing_over_thorax ~ Sex*Medium*Parental_age + (1|Vial_no_corrected) + (1|Individual_ID), family = "gaussian", data = Wingth)

s0 <- simulateResiduals(m0)
plot(s0)  

testDispersion(s0) ## Not over/ underdispersed

summary(m0)
car::Anova(m0) ## keep sex* parental age and medium* parental age but drop 4-way 


m1 <- glmmTMB(Wing_over_thorax ~ Sex*Medium + Sex*Parental_age + Parental_age*Medium + (1|Vial_no_corrected), family = "gaussian", data = Wingth)

s1 <- simulateResiduals(m1)
plot(s1) ## 

testDispersion(s1) ## Not over/ underdispersed

summary(m1)
car::Anova(m1)


## Drop sex*medium

m2 <- glmmTMB(Wing_over_thorax ~ Sex*Parental_age + Parental_age*Medium + (1|Vial_no_corrected), family = "gaussian", data = Wingth)

s2 <- simulateResiduals(m2)
plot(s2) ## Slightly wonky  - probably fine? 

testDispersion(s2) ## Not over/ underdispersed

summary(m2)
car::Anova(m2)

AIC(m0, m1, m2)

## M2 wins and can't be further simplified

#### Now do tukeys in emmeans

### Need to model separately for each sex to look at the sex by parental age effect 

Fs <- filter(Wingth, Sex=="F")
Ms <- filter(Wingth, Sex=="M")

## For Fs
m2_F <- glmmTMB(Wing_over_thorax ~ Parental_age*Medium + (1|Vial_no_corrected), family = "gaussian", data = Fs)

emm_thorax_Fs <- emmeans(m2_F, "Medium", "Parental_age", data=Fs)
pairs(emm_thorax_Fs)

emm_thorax_Fs <- emmeans(m2_F,  "Parental_age", "Medium", data=Fs)
pairs(emm_thorax_Fs)

## For Ms

m2_M <- glmmTMB(Wing_over_thorax ~ Parental_age*Medium + (1|Vial_no_corrected), family = "gaussian", data = Ms)

emm_thorax_Ms <- emmeans(m2_M, "Medium", "Parental_age", data=Ms)
pairs(emm_thorax_Ms)

emm_thorax_Ms <- emmeans(m2_M,  "Parental_age", "Medium", data=Ms)
pairs(emm_thorax_Ms)

Plot

