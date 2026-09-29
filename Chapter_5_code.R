###################################### Scripts for chapter 5 

###################################### Mating latency in females

######
## usual friends
library(ggplot2)
library(scales)
library(gridExtra)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
######
## Survival packages 
library(survival)
library(survminer) # for Kaplan-Meier plots
library(broom) # for tidy output 
library(ggplot2) # (loaded by survminer but I typed it again out of habit)
library(coxme) # bet it'll be a bad fit IT'S ALWAYS A BAD FIT
######

Massays <- read.csv("Massays_offs.csv") # "Mating_assays_offspring_counts.csv" in data file
head(Massays)

Massays_female <- filter(Massays, Sex=="Female")

head(Massays_female)

#########################

KM_fit <- survfit(Surv(Mating_latency, Mated) ~ Diet_group + Parent_age + Individual_age, 
                  data = Massays_female)
Massays_female$Individual_age <- factor(Massays_female$Individual_age, levels = c("Young", "Old"))

Plot_lat <- ggsurvplot_facet(
  KM_fit,
  data = Massays_female,
  fun = "event",
  conf.int = TRUE,
  xlim = c(0, 90),
  break.x.by = 10,
  xlab = "Observation window (minutes)",
  ylim = c(0.00, 1.00),
  ylab = "Probability of mating",
  facet.by = c("Diet_group", "Individual_age"),
  legend.title = "Parent age",
  palette = c("#2d6647","#61b888"),
  ggtheme = theme_bw()
)
Plot_lat


#################

Cox_fit <- coxme(
  Surv(Mating_latency, Mated) ~ Diet_group * Parent_age * Individual_age  + (1 | Vial_ID),
  data = Massays_female
)

summary(Cox_fit)
cox.zph(Cox_fit)


#### mating latency in females DOES NOT meet the proportional hazards assumption, so moving on to event history analysis 

### make data long

Massays_female <- Massays_female %>%
  mutate(Mating_latency = ifelse(is.na(Mating_latency) & Mated == 0, 
                                 90, Mating_latency))

Massays_long_female <- Massays_female %>%
  rowwise() %>%
  do({
    row <- .
    data.frame(
      Vial_ID = row$Vial_ID,
      Parent_age = row$Parent_age,
      Sex = row$Sex,
      Diet_group = row$Diet_group,
      Individual_age = row$Individual_age,
      time = 1:row$Mating_latency,
      event = c(rep(0, row$Mating_latency - 1), row$Mated)
    )
  }) %>%
  ungroup()

####

glmm_time <- glmmTMB(event ~ time + Parent_age*Individual_age*Diet_group + (1|Vial_ID), family = binomial, data = Massays_long_female)
summary(glmm_time)

car::Anova(glmm_time, type = "III")
s0 <- simulateResiduals(glmm_time)
plot(s0)
testDispersion(s0)
### looks good

#### now Tukeys 


#### Mother's age * individual age 
em1 <- emmeans(glmm_time, ~ Individual_age * Parent_age, data = Massays_long_female)
pairs(em1, adjust = "tukey")
plot(em1)

### Mother's age * age * diet group
em1 <- emmeans(glmm_time, ~ Individual_age * Parent_age | Diet_group, data = Massays_long_female)
pairs(em1, adjust = "tukey")
plot(em1)


######################################### Mating latency in males

Massays_male <- filter(Massays, Sex=="Male")

head(Massays_male)
##############

KM_fit <- survfit(Surv(Mating_latency, Mated) ~ Diet_group + Parent_age + Individual_age, 
                  data = Massays_male)
Massays_male$Individual_age <- factor(Massays_male$Individual_age, levels = c("Young", "Old"))

Plot_lat <- ggsurvplot_facet(
  KM_fit,
  data = Massays_male,
  fun = "event",
  conf.int = TRUE,
  xlim = c(0, 90),
  break.x.by = 10,
  xlab = "Observation window (minutes)",
  ylim = c(0.00, 1.00),
  ylab = "Probability of mating",
  facet.by = c("Diet_group", "Individual_age"),
  legend.title = "Parent age",
  palette = c("#a8682b","#D4AC3E"),
  ggtheme = theme_bw()
)
Plot_lat

## colours before: #4c57ad","#5edbcf", "#d968d3", "#623670
#################

Cox_fit <- coxme(
  Surv(Mating_latency, Mated) ~ Diet_group * Parent_age * Individual_age  + (1 | Vial_ID),
  data = Massays_male
)

summary(Cox_fit)
cox.zph(Cox_fit)

#### Mating latency in males only meets proportional hazards!!! 

car::Anova(Cox_fit, type = "III")

#### Tukeys

## Individual age

em1 <- emmeans(Cox_fit, "Individual_age", data = Massays_male)
pairs(em1, adjust = "tukey")
plot(em1)

## Parent age * individual age
em1 <- emmeans(Cox_fit, ~ Individual_age * Parent_age, data = Massays_male)
pairs(em1, adjust = "tukey")
plot(em1)

## Diet group * parent age * individual age

em1 <- emmeans(Cox_fit, ~ Individual_age * Parent_age | Diet_group, data = Massays_male)
pairs(em1, adjust = "tukey")
plot(em1)

## other view 

em1 <- emmeans(Cox_fit, ~ Diet_group * Individual_age  | Parent_age , data = Massays_male)
pairs(em1, adjust = "tukey")
plot(em1)

################################################################

############################## Mating duration in females

head(Massays_female)
##############

KM_fit <- survfit(Surv(Mating_duration, Mated) ~ Diet_group + Parent_age + Individual_age, 
                  data = Massays_female)
Massays_female$Individual_age <- factor(Massays_female$Individual_age, levels = c("Young", "Old"))

Plot_dur <- ggsurvplot_facet(
  KM_fit,
  data = Massays_female,
  fun = "event",
  conf.int = TRUE,
  xlim = c(0, 90),
  break.x.by = 10,
  xlab = "Observation window (minutes)",
  ylim = c(0.00, 1.00),
  ylab = "Probability of uncoupling from mating onset",
  facet.by = c("Diet_group", "Individual_age"),
  legend.title = "Parent age",
  palette = c("#2d6647","#61b888"),
  ggtheme = theme_bw()
)
Plot_dur

####################

Cox_fit <- coxme(
  Surv(Mating_duration, Mated) ~ Diet_group * Parent_age * Individual_age  + (1 | Vial_ID),
  data = Massays_female
)

summary(Cox_fit)
cox.zph(Cox_fit)

### Does not violate proportional hazards!!! 

### significant 3-way so can't simplify further

car::Anova(Cox_fit, type = "III")


### Tukey
em1 <- emmeans(Cox_fit, ~ Individual_age * Parent_age | Diet_group, data = Massays_female)
pairs(em1, adjust = "tukey")
plot(em1)


em2 <- emmeans(Cox_fit, ~ Individual_age * Parent_age, data = Massays_female)
pairs(em2, adjust = "tukey")
plot(em2)

################################################# Mating duration in males

head(Massays_male)
##############

KM_fit <- survfit(Surv(Mating_duration, Mated) ~ Diet_group + Parent_age + Individual_age, 
                  data = Massays_male)
Massays_male$Individual_age <- factor(Massays_male$Individual_age, levels = c("Young", "Old"))

Plot_dur <- ggsurvplot_facet(
  KM_fit,
  data = Massays_male,
  fun = "event",
  conf.int = TRUE,
  xlim = c(0, 90),
  break.x.by = 10,
  xlab = "Observation window (minutes)",
  ylim = c(0.00, 1.00),
  ylab = "Probability of uncoupling from mating onset",
  facet.by = c("Diet_group", "Individual_age"),
  legend.title = "Parent age",
  palette = c("#a8682b","#D4AC3E"),
  ggtheme = theme_bw()
)
Plot_dur

####################

Cox_fit <- coxme(
  Surv(Mating_duration, Mated) ~ Diet_group * Parent_age * Individual_age  + (1 | Vial_ID),
  data = Massays_male
)

summary(Cox_fit)
cox.zph(Cox_fit)

#### violates proportional hazards assumption


#### Make data long then do binomial 

### exclude all unmated so no NAs

Massays_male <- filter(Massays_male, Mated=="1")

### two cases where time of uncoupling was missed - need to get rid 

Massays_male <- Massays_male %>%
  filter(!(Vial_ID %in% c("#337", "#859")))

### need "time of uncoupling" column instead of time to mating 

Massays_long_male <- Massays_male %>%
  rowwise() %>%
  do({
    row <- .
    data.frame(
      Vial_ID = row$Vial_ID,
      Parent_age = row$Parent_age,
      Sex = row$Sex,
      Diet_group = row$Diet_group,
      Individual_age = row$Individual_age,
      time = 1:row$Mating_duration,
      event = c(rep(0, row$Mating_duration - 1), 1)
    )
  }) %>%
  ungroup()

head(Massays_long_male)

####

glmm_time <- glmmTMB(event ~ time + Parent_age*Individual_age*Diet_group + (1|Vial_ID), family = binomial, data = Massays_long_male)
summary(glmm_time)

car::Anova(glmm_time) 
s0 <- simulateResiduals(glmm_time)
plot(s0)
testDispersion(s0)

### Can simplify this a lot, no sig interactions  

glmm_time <- glmmTMB(event ~ time + Parent_age + Individual_age + Diet_group + (1|Vial_ID), family = binomial, data = Massays_long_male)
summary(glmm_time)

car::Anova(glmm_time) 
s0 <- simulateResiduals(glmm_time)
plot(s0)
testDispersion(s0)


## Tukeys for individual age and diet group 

em1 <- emmeans(glmm_time, "Individual_age", data = Massays_long_male)
pairs(em1, adjust = "tukey")
plot(em1)


em1 <- emmeans(glmm_time, "Diet_group", data = Massays_long_male)
pairs(em1, adjust = "tukey")
plot(em1)


########################################################## 

##############################  Offspring number 

Data <- read.csv("Mating_assays_offspring_counts.csv") # same as in data file

head(Data)

Data <- Data %>% mutate(Parent_age = factor(Parent_age)) 
Data <- Data %>% mutate(Diet_group = factor(Diet_group))
Data <- Data %>% mutate(Individual_age = factor(Individual_age))

###########################################

Males <- filter(Data, Sex=="Male")
Females <- filter(Data, Sex=="Female")

########################### For plots

Age.order <- c('Young', 'Old')
M.order <- c('DR-DR', 'FF-DR', 'DR-FF','FF-FF')

########################################## Final models

############################### Final model for males - some deviation from predicted values, model selection showed that truncated compois was the best family (lowest AIC) and trial and error led to the ziformula and dispformula shown. Each model iteration was checked by visual inspection of residuals and comparison of AIC 

male_offs <- glmmTMB(
  Offspring ~ Diet_group * Individual_age * Parent_age, 
  ziformula = ~ Individual_age*Parent_age + Diet_group, 
  dispformula = ~ Individual_age*Parent_age + Diet_group,
  family = truncated_compois, 
  data = Males
)

m_res <- simulateResiduals(male_offs)
plot(m_res)
testDispersion(m_res)
testZeroInflation(m_res)
summary(male_offs)
car::Anova(male_offs, type = "III")

############################ Final model for females - fit with no problems after similar effort

female_offs <- glmmTMB(
  Offspring ~ Diet_group * Individual_age * Parent_age, 
  ziformula = ~ Diet_group + Individual_age + Parent_age, 
  dispformula = ~ Diet_group + Individual_age,
  family = nbinom2, 
  data = Females
)

f_res <- simulateResiduals(female_offs)
plot(f_res)
testDispersion(f_res)
testZeroInflation(f_res)
summary(female_offs)
car::Anova(female_offs, type = "III")

############################################################################## Emmeans Tukey tests
##############################################################################


############################### Males
###############################

emms_m <- emmeans(male_offs, ~ Parent_age * Individual_age | Diet_group,
                  type = "response")
pairs(emms_m)


############################## Females
##############################

emms_f <- emmeans(female_offs, ~ Parent_age * Individual_age | Diet_group,
                  type = "response")
pairs(emms_f)

### And just the 2-way because 3 way was marginally n.s.

emms_f <- emmeans(female_offs, ~ Parent_age * Individual_age,
                  type = "response")
pairs(emms_f)




############ Plots with emmeans 


########################################## Males
##########################################

emms_m <- emmeans(male_offs, ~ Individual_age + Diet_group | Parent_age,
                  type = "response")

emm_df_m <- as.data.frame(emms_m)

class(emm_df_m)

Plot <- ggplot(Males, aes(
  x      = Diet_group,
  y      = Offspring,
  colour = Parent_age
)) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width  = 0.2,
      jitter.height = 0.4,
      dodge.width   = 0.75
    ),
    alpha = 0.4,
    size  = 1.5
  ) +
  geom_errorbar(
    data = emm_df_m,
    aes(
      x     = Diet_group,
      y     = response,
      ymin  = asymp.LCL,
      ymax  = asymp.UCL,
      group = Parent_age
    ),
    width     = 0.2,
    linewidth = 1.2,
    position  = position_dodge(width = 0.75)
  ) +
  geom_point(
    data = emm_df_m,
    aes(
      x     = Diet_group,
      y     = response,
      group = Parent_age
    ),
    size     = 3,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~Individual_age) +
  scale_y_continuous(
    limits = c(0, 250),
    breaks = scales::pretty_breaks(),
    oob    = scales::squish
  ) +
  scale_colour_manual(values = c("#d4ac3e", "#a8682b"), labels = c("Young" = "Young mother (5d)", "Old" = "Old mother (35d)")) +
  labs(
    x     = "Mother's diet - Son's diet",
    y     = "Offspring (over 4 days with two 5d unmated females)",
    color = "Parent_age"
  ) +
  theme_bw()
Plot

##################### ordering

# Fix the main data
Males$Parent_age     <- factor(Males$Parent_age, levels = c("Young", "Old"))
Males$Individual_age <- factor(Males$Individual_age, levels = c("Young", "Old"))

# Fix the summary data
emm_df_m$Parent_age     <- factor(emm_df_m$Parent_age, levels = c("Young", "Old"))
emm_df_m$Individual_age <- factor(emm_df_m$Individual_age, levels = c("Young", "Old"))


Plot <- ggplot(Males, aes(
  x      = Diet_group,
  y      = Offspring,
  colour = Parent_age   # This now controls the young/old order on the left/right dodge
)) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width  = 0.2,
      jitter.height = 0.4,
      dodge.width   = 0.75
    ),
    alpha = 0.2,
    size  = 1.5
  ) +
  geom_errorbar(
    data = emm_df_m,
    aes(
      x      = Diet_group,
      y      = response,
      ymin   = asymp.LCL,
      ymax   = asymp.UCL,
      colour = Parent_age,  # Explicitly map colour here so it matches the new factor levels
      group  = Parent_age
    ),
    width     = 0.2,
    linewidth = 1.2,
    position  = position_dodge(width = 0.75)
  ) +
  geom_point(
    data = emm_df_m,
    aes(
      x      = Diet_group,
      y      = response,
      colour = Parent_age,  # Explicitly map colour here so it matches the new factor levels
      group  = Parent_age
    ),
    size     = 3,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~Individual_age) +  # Facets will now display 'Young' left and 'Old' right
  scale_y_continuous(
    limits = c(0, 250),
    breaks = scales::pretty_breaks(),
    oob    = scales::squish
  ) +
  scale_colour_manual(
    values = c("Young" = "#d4ac3e", "Old" = "#a8682b"), # Named values prevent color swapping mismatches
    labels = c("Young" = "Young mother (5d)", "Old" = "Old mother (35d)")
  ) +
  labs(
    x     = "Mother's diet - Son's diet",
    y     = "Offspring (over 4 days with two 5d unmated females)",
    color = "Parent age"
  ) +
  theme_bw()

Plot


################################################### Females
###################################################

emms_f <- emmeans(female_offs, ~ Individual_age + Diet_group | Parent_age,
                  type = "response")

emm_df_f <- as.data.frame(emms_f)


# Fix the main data
Females$Parent_age     <- factor(Females$Parent_age, levels = c("Young", "Old"))
Females$Individual_age <- factor(Females$Individual_age, levels = c("Young", "Old"))

# Fix the summary data
emm_df_f$Parent_age     <- factor(emm_df_f$Parent_age, levels = c("Young", "Old"))
emm_df_f$Individual_age <- factor(emm_df_f$Individual_age, levels = c("Young", "Old"))


Plot <- ggplot(Females, aes(
  x      = Diet_group,
  y      = Offspring,
  colour = Parent_age   # This now controls the young/old order on the left/right dodge
)) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width  = 0.2,
      jitter.height = 0.4,
      dodge.width   = 0.75
    ),
    alpha = 0.2,
    size  = 1.5
  ) +
  geom_errorbar(
    data = emm_df_f,
    aes(
      x      = Diet_group,
      y      = response,
      ymin   = asymp.LCL,
      ymax   = asymp.UCL,
      colour = Parent_age,  # Explicitly map colour here so it matches the new factor levels
      group  = Parent_age
    ),
    width     = 0.2,
    linewidth = 1.2,
    position  = position_dodge(width = 0.75)
  ) +
  geom_point(
    data = emm_df_f,
    aes(
      x      = Diet_group,
      y      = response,
      colour = Parent_age,  # Explicitly map colour here so it matches the new factor levels
      group  = Parent_age
    ),
    size     = 3,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~Individual_age) +  # Facets will now display 'Young' left and 'Old' right
  scale_y_continuous(
    limits = c(0, 150),
    breaks = scales::pretty_breaks(),
    oob    = scales::squish
  ) +
  scale_colour_manual(
    values = c("Young" = "#61b888", "Old" = "#2d6647"), # Named values prevent color swapping mismatches
    labels = c("Young" = "Young mother (5d)", "Old" = "Old mother (35d)")
  ) +
  labs(
    x     = "Mother's diet - Daughter's diet",
    y     = "Offspring (over 4 days with two 5d unmated males)",
    color = "Parent age"
  ) +
  theme_bw()

Plot

############################################################# Lifespan - P0 (mothers)

################## 


library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(BaSTA)
library(tidyverse)
library(tidylog)
library(coxme)

LS <- read.csv("p0_LS.csv")
head(LS)

###################################### KM plot

KM_fit <- survfit(Surv(Death_age, Died) ~ Diet_group, data = LS)

Plot_lat <- ggsurvplot(
  KM_fit,
  data = LS,
  conf.int = FALSE,
  xlim = c(0, 100),
  break.x.by = 10,
  xlab = "Time (days)",
  ylim = c(0.00, 1.00),
  ylab = "Survival",
  legend.title = "Diet",
  palette = c("#6CACC4","#1770A3"),
  ggtheme = theme_bw()
)
Plot_lat
###################################### Coxph out of interest

Cox_fit <- coxme(
  Surv(Death_age, Died) ~ Diet_group  + (1 | Vial_uncorrected),
  data = LS)

summary(Cox_fit)
cox.zph(Cox_fit) ### high global p as usual


########################################### Lifespan analysis with survreg - tried BaSTA but model repeatedly failed to converge with no informative error messages to diagnose the problem. For P0 my interest is only in whether DR significantly extended lifespan so decided that simpler weibull analysis using survreg was sufficient for P0. 

Weibullmodel<-survreg(Surv(Death_age,Died)~Diet_group,data=LS,dist='weibull')
summary(Weibullmodel)

library(emmeans)
em <- emmeans(Weibullmodel, ~ Diet_group, data = LS)
pairs(em, adjust = "tukey")
plot(em)

WeibullDiag(Surv(Death_age,Died)~Diet_group,data=LS) 


################################################### F1 Lifespan analysis with BaSTA

parent_dat <- read_csv("F1_LS_Basta_census.csv")

census_dat <- data.frame(ID = 1:nrow(parent_dat),
                         Birth.Date = parent_dat$Birth.Date,
                         Min.Birth.Date = parent_dat$Birth.Date,
                         Max.Birth.Date = parent_dat$Birth.Date,
                         Entry.Date = parent_dat$Birth.Date,
                         Depart.Date = parent_dat$Depart.Date,
                         Depart.Type = "D",
                         treatment = paste(parent_dat$Sex, parent_dat$Diet.Group, parent_dat$Mother.Age))

### multibasta

bas_1 <- multibasta(census_dat, dataType = "census",
                    nsim = 4, parallel = TRUE, ncpus = 10,
                    niter = 150000, burnin = 15001, thinning = 150,
                    formulaMort = ~ treatment -1)

summary(bas_1, digits = 3)

### Weibull simple is the best model

bas_2 <- basta(census_dat, dataType = "census",
               model = "WE", shape = "simple",
               nsim = 4, parallel = TRUE, ncpus = 10,
               niter = 150000, burnin = 15001, thinning = 150,
               formulaMort = ~ treatment -1)

summary(bas_2, digits = 3)


plot(bas_2, type = "fancy", plot.type = "demorates")

####################################
####################################

## Now aiming for individual plots for each sex and diet group, with parental ages in different colours, like the KM plots

#### Filtering for different groups 

### Mixing sexes might be fine? test with DRDR first

DRDR <- filter(parent_dat, Diet.Group =="DR-DR")

## reformat as census
census_DRDR <- data.frame(ID = 1:nrow(DRDR),
                          Birth.Date = DRDR$Birth.Date,
                          Min.Birth.Date = DRDR$Birth.Date,
                          Max.Birth.Date = DRDR$Birth.Date,
                          Entry.Date = DRDR$Birth.Date,
                          Depart.Date = DRDR$Depart.Date,
                          Depart.Type = "D",
                          treatment = paste(DRDR$Sex, DRDR$Mother.Age))


## Weibull simple on just DRDR

bas_DRDR <- basta(census_DRDR, dataType = "census",
                  model = "WE", shape = "simple",
                  nsim = 4, parallel = TRUE, ncpus = 10,
                  niter = 150000, burnin = 15001, thinning = 150,
                  formulaMort = ~ treatment -1)

summary(bas_DRDR, digits = 3)

plot(bas_DRDR, type = "fancy", densities = TRUE, col = c("#2d6647","#61b888", "#a8682b", "#d4ac3e")) 


################################ 

FFFF <- filter(parent_dat, Diet.Group =="FF-FF")
## reformat as census
census_FFFF <- data.frame(ID = 1:nrow(FFFF),
                          Birth.Date = FFFF$Birth.Date,
                          Min.Birth.Date = FFFF$Birth.Date,
                          Max.Birth.Date = FFFF$Birth.Date,
                          Entry.Date = FFFF$Birth.Date,
                          Depart.Date = FFFF$Depart.Date,
                          Depart.Type = "D",
                          treatment = paste(FFFF$Sex, FFFF$Mother.Age))


## Weibull simple on just FFFF

bas_FFFF <- basta(census_FFFF, dataType = "census",
                  model = "WE", shape = "simple",
                  nsim = 4, parallel = TRUE, ncpus = 10,
                  niter = 150000, burnin = 15001, thinning = 150,
                  formulaMort = ~ treatment -1)

summary(bas_FFFF, digits = 3)

plot(bas_FFFF, type = "fancy", densities = TRUE, col = c("#2d6647","#61b888", "#a8682b", "#d4ac3e")) 

####################################

DRFF <- filter(parent_dat, Diet.Group =="DR-FF")
## reformat as census
census_DRFF <- data.frame(ID = 1:nrow(DRFF),
                          Birth.Date = DRFF$Birth.Date,
                          Min.Birth.Date = DRFF$Birth.Date,
                          Max.Birth.Date = DRFF$Birth.Date,
                          Entry.Date = DRFF$Birth.Date,
                          Depart.Date = DRFF$Depart.Date,
                          Depart.Type = "D",
                          treatment = paste(DRFF$Sex, DRFF$Mother.Age))


## Weibull simple on just DRFF

bas_DRFF <- basta(census_DRFF, dataType = "census",
                  model = "WE", shape = "simple",
                  nsim = 4, parallel = TRUE, ncpus = 10,
                  niter = 150000, burnin = 15001, thinning = 150,
                  formulaMort = ~ treatment -1)

summary(bas_DRFF, digits = 3)

plot(bas_DRFF, type = "fancy", densities = TRUE, col = c("#2d6647","#61b888", "#a8682b", "#d4ac3e")) 


#######################################
FFDR <- filter(parent_dat, Diet.Group =="FF-DR")
## reformat as census
census_FFDR <- data.frame(ID = 1:nrow(FFDR),
                          Birth.Date = FFDR$Birth.Date,
                          Min.Birth.Date = FFDR$Birth.Date,
                          Max.Birth.Date = FFDR$Birth.Date,
                          Entry.Date = FFDR$Birth.Date,
                          Depart.Date = FFDR$Depart.Date,
                          Depart.Type = "D",
                          treatment = paste(FFDR$Sex, FFDR$Mother.Age))


## Weibull simple on just DRFF

bas_FFDR <- basta(census_FFDR, dataType = "census",
                  model = "WE", shape = "simple",
                  nsim = 4, parallel = TRUE, ncpus = 10,
                  niter = 150000, burnin = 15001, thinning = 150,
                  formulaMort = ~ treatment -1)

summary(bas_FFDR, digits = 3)

plot(bas_FFDR, type = "fancy", densities = TRUE, col = c("#2d6647","#61b888", "#a8682b", "#d4ac3e")) 

