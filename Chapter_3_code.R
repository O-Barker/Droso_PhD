####################################### All scripts for Chapter 3 combined
######### This script gives a brief overview of final models and model selection without the trial and error which led to them - these were the first analyses of my PhD, and original scripts involved a lot of data exploration and trying new things. For readability, rather than recapitulating the learning process which led to these final models, comments summarise what was tried. 

########################## Sucrose-casein diet choice in adult males and females

library(glmmTMB)
library(tidyr)
library(DHARMa)
library(hablar)
library(dplyr)

Diet <- read.csv("Diet_choice.csv") ## "Adult_diet_choice_synthetic" in data file
head(Diet)

hist(Diet$Flies)
### Veeery zeroinflated

var(Diet$Flies)/mean(Diet$Flies) ## Var > mean so Poisson unsuitable

## Made everything a factor - repeated this step for all data 

Diet <- Diet %>% mutate(Concentration=(factor(Concentration)))
Diet <- Diet %>% mutate(Sex=(factor(Sex)))
Diet <- Diet %>% mutate(Patch=(factor(Patch)))
Diet <- Diet %>% mutate(Plate=(factor(Plate)))
Diet <- Diet %>% mutate(Observation_interval=(factor(Observation_interval)))

## Tried nbinom1, nbinom2, truncated versions of both, and compois. Nbinom1 had the lowest AIC

M1 <- glmmTMB(No_flies ~ Patch*Concentration*Treatment*Sex + Block + (1|Plate) + (1|Observation_interval) + (1|obs), data = Diet, family = "nbinom1", ziformula=~1)

S1 <- simulateResiduals(M1)
plot(S1) 

summary(M1)
car::Anova(M1)

############################### Sucrose-casein Larval diet choice

Larv <- read.csv("Larval_choice.csv") # "Larval_diet_choice_synthetic" in data file
head(Larv)

### Began with 3-way interaction, best fit/ lowest AIC achieved with the following: 

M2 <- glmmTMB(No_larvae ~ Patch*Concentration + Patch*Block + Block*Concentration + (1|Replicate), zi=~1, family = "nbinom2", control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")), data = Larv1)
s2 <- simulateResiduals(M2)
plot(s2)
testDispersion(s2) 

car::Anova(M2)
summary(M2)

############################# Sucrose-casein oviposition choice

Ovi <- read.csv("Oviposition_choice_synthetic.csv") # Same in data file

hist(Ovi$Eggs)
## Very zeroinflated, negative binomial for ovi too 
M0 <- glmmTMB(Eggs ~ Patch_PCB*Concentration*Block + (1|Replicate), zi=~1, family = "nbinom1", data = Ovi)
s0 <- simulateResiduals(M0)
plot(s0) ### no problems
testDispersion(s0) ## 
testZeroInflation(s0) ## 

car::Anova(M0)
summary(M0)

############################# Female diet choice for agar content, with dyes

Diet <- read.csv("Female_Diet_Choice_SYA_Dyes.csv") # Same in data file

## Tried Poisson with and without ziformula, nbinom, truncated poisson and compois, Poisson had lowest AIC  


M1 <- glmmTMB(Chosen ~ Agar_content + Colour + (1|Plate) + (1|Obs), family = "poisson", data = Diet)
s1 <- simulateResiduals(M1)
plot(s1)
testZeroInflation(s1) 
testDispersion(s1) 

# bbmle::AICctab(M1, M2, M3, M4)


########################### Larval diet choice for agar content, with dyes

## Poisson was a good fit, tried with and without the interaction term and the observation-level random effect, and the model including both had the lowest AIC.

Larv2 <- read.csv("Larval_diet_choice_SYA_dyes.csv") # Same in data file

M3 <- glmmTMB(Larvae_choice ~ Agar_content*Colour + (1|Plate) + (1|Obs), family = "poisson", data = Larv2)
s3 <- simulateResiduals(M3)
plot(s3) 

# AICtab(M1, M2, M3)
## M3 and M2 identical 

summary(M3)
car::Anova(M3)


######################### Oviposition choice for agar content, without the use of dyes

Ovi2 <- read.csv("Oviposition_choice_SYA_no_dyes.csv")
head(Ovi2)

## The best fit was a simple poisson with nested random effects for plate and observation. Overdispersion was fixed very neatly by including the observation-level random effect. 

m0 <- glmmTMB(Eggs_total ~ Diet + (1|Plate) + (1|Obs),  data = Ovi, family = 'poisson')

s0 <- simulateResiduals(m0)
plot(s0) 
testDispersion(s0) 

########################### Tukey tests

## For all models, pairwise comparisons resembled the following:

### e.g. for diet choice for agar content with the use of dyes in females:
emm_diet <- emmeans(M1, "Agar_content", "Colour",  data=Diet)
pairs(emm_diet, adjust="tukey")

plot(emm_diet, comparisons = TRUE)
