library(ggplot2)
library(plyr)
library(dplyr)
library(tibble)
library(nvmix)
library(readxl)
library(pracma)
library(fit.models)
library(sqldf)
library(writexl)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(readr)
library(data.table)
library(dygraphs)
library(xts)
library(lubridate)
library(plotly)
library(hrbrthemes)
library(forecast)
library(bnlearn)
library(R.utils)
library(Rgraphviz)
library(gRain)
library(RBGL)
library(graph)
library(RCurl)
library(survival)
library(survminer)

################################################################
########################## WEEK 1 ##############################
### Import the dataset ###
g <- read.csv("C:/Users/liry9/Desktop/Tesi/Survival Analysis in R for Public Health_Imperial College/Week 1/simulated_HF.csv", header=TRUE, sep=',')
#g <- na.omit(g)
dim(g)
head(g)

#### Conversion of variables into categorical ###
gender <- as.factor(g[,"gender"]) # R calls categorical variables factors
fu_time <- g[,"fu_time"] # continuous variable (numeric) 
death <- g[,"death"] # binary variable (numeric) 

### Run an overall Kaplan-Meier plot ###
km_fit <- survfit(Surv(fu_time, death) ~ 1)

plot(km_fit)

## "Survfit" fits a simple survival model that depends only on
# gender in terms of predictors. in this case there aren't any
# predictors, so the model just has the intercept, denoted by
# "1". The two arguments used by "Surv" are the follow-up time
# for each patient and whether they died.
# In our data, death=1 for people who had died by the end of
# the follow-up period, and death=0 for those still alive at
# that time.
## Technically, those people still alive are censored,
# because we don't know when they'll die (of course everyone
# does at some point). The survfit() function produces the
# Kaplan-Meier estimates of the probability of survival over
# time that are used by "plot" to produce the Kaplan-Meier
# curve above.
# These estimates can be seen by typing:
summary(km_fit, times = c(1:7,30,60,90*(1:10)))

# The "time" argument gives us control over what time
# periods we want to see. The above code asks for output
# every day for the first week, then at 30, 60 and 90 days,
# and then every 90 days thereafter.

# Whereas all but about 1% make it past the first day, at
# 900 days after a first emergency admission for heart
# failure, the probability of surviving is just 38%.

## Let's extend this by splitting the curve by gender: 
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

plot(km_gender_fit)

### To compare survival by gender, we can run a logrank test
# or Mantel-Haenszel test.
# There are many ways to do this because of different
# versions for different scenarios, e.g. particular types
# of censored data, the standard one is:
survdiff(Surv(fu_time, death) ~ gender, rho=0)

## The function returns a list of components, including:
# - n: the number of subjects in each group.
# - obs: the weighted observed number of events in each group.
# - exp: the weighted expected number of events in each group.
# - chisq: the chisquare statistic for a test of equality.
# - strata: optionally, the number of subjects contained in each stratum.

# Just to be sure you've got that, I'd like you to compare
# survival by broad age group: those aged 65 and above versus
# those aged under 65. To do this, you'll need to dichotomise
# age, for instance using the "ifelse" command, which was
# covered in the first course in this specialisation.
## My solution:
#g <- g %>% mutate(category=cut(age, breaks=c(-Inf, 65, Inf), labels=c(0,1)))
# 0 = under 65, 1 = over 65 y.o.
#colnames(g)[32] <- "age2"
#age2 <- as.factor(g[,"category"])
#summary(age2)
### Run the Kaplan-Meier plot ###
#km_age_fit <- survfit(Surv(fu_time, death) ~ age2) 
#plot(km_age_fit)

# Log-rank test #
#survdiff(Surv(fu_time, death) ~ age2, rho=0)

## Codice alternativo
age_65plus <- ifelse(g[,"age"]>=65,1,0) # dichotomise age
table(age_65plus, exclude = NULL) # inspect the numbers - always a good idea
age_65plus
table(g$age,age_65plus, exclude = NULL) #
survdiff(Surv(fu_time, death) ~ age_65plus, rho=0)

###############################################################
######################### WEEK 2 ##############################
#### COX METHOD ####
### Age as the only predictor
## we can read the follow-up time variable and the outcome
# variable, namely death
cox_age <- coxph(Surv(fu_time, death) ~ age, data=g)
summary(cox_age)

## the coefficients are the hazard ratios. 1,06 means that for
# each increase of 1 year in age, the hazard of death at any given
# time point, since the beginning of the study goes up by 6%
# Every time you get a year older, your hazard goes up by 6%.
# se(coef) is the standard error = tiny and it is typical when we
# have only one predictor and many patients.
# z value = similar to t value and converted to p value: in our case
# it is very tiny hence very significant.

## Then we find the confidence levels and some technical features:
# - concordance =  is rather like the discrimination or c-statistic
# in logistic regression it's the fraction of pairs of patients
# in a sample where the patient with the highest survivor time
# has the higher probability of survival as predicted by the
# model -> high values are good and what interest us.
# - r-squared = the proportion of variation in the outcome
# that is explained by the model. Again high values are good.
# - likelihood ratio test = close to 0 -> null hyp not accepted
# - wold test = = close to 0 -> null hyp not accepted
# - score or log-rank test = = close to 0 -> null hyp not accepted
## All this means that the model with age in it is significantly better than
# w/o it, namely one with nothing in it.

## we run the same code as before for the ethnicity and we
# will get one coefficient for ethnicgroup. This is because
# unless you tell it otherwise, R will assume that all your
# variables are continuous. Ethnicity is very obviously not
# a continuous variable, but R doesn't know that unless you
# tell it!
ethnicgroup <- factor(g[,"ethnicgroup"]) # can also use "as.factor" rather than "factor"

cox_ethnic <- coxph(Surv(fu_time, death) ~ ethnicgroup)
summary(cox_ethnic)

###################################################################
############################ WEEK 3 ###############################
summary(g$age)

# For prior DNAs, i.e. prior missed outpatient appointments,
# the "summary" command isn't very useful.
# More informative is the "table" command:  
summary(g$prior_dnas)
table(g$prior_dnas)
t_prior_dnas <- table(g$prior_dnas, exclude=NULL) 
addmargins(t_prior_dnas) # adds the total (a "sum" column) 
round(100*prop.table(t_prior_dnas),digits=1)
# So nearly three in four patients had missed no appointments,
# but nearly three percent had missed five or more, with a
# maximum of ten.

# COPD as a comorbidity was pretty straightforward too:
table(g$copd, exclude=NULL)
# 24% of patients had COPD, with no missing values - though
# remember what I said in the earlier video about missing
# values masquerading as regular values. It's actually likely
# that some patients have COPD but haven't been recorded as
# having it. Such underrecording of comorbidities is common
# with administrative data for various reasons.

# Also, note that both categories have lots of patients in.
table(g$gender, exclude=NULL)
t <- table(gender, exclude=NULL)
addmargins(t) # adds the total (a "sum" column)
round(100*prop.table(t),digits=1) #no missing values.

table(g$ethnicgroup, exclude=NULL)

# You've already seen age before, so you can get the same
# histogram / kernel density plot as before. The main things
# to note here are that there's a pretty wide spread of values
# and none is missing. 

##### How to run multiple Cox model in R #####
cox_multiple1 <- coxph(Surv(fu_time, death) ~ g$age + gender + g$copd + g$prior_dnas + ethnicgroup)

summary(cox_multiple1)

### Cox model including age, gender, copd, quintile, and ethnic group
quintile <- factor(g[,"quintile"])
cox_multiple2 <- coxph(Surv(fu_time, death) ~ g$age + gender + g$copd + quintile + ethnicgroup)

summary(cox_multiple2)

### Option 1: Change the quintile reference
quintile <- relevel(quintile, ref = 1) # make the reference category quintile=1

### Option 2 and 3: To combine quintile values, you can do this:
## best start again with the original data set, not from the
# existing object called "quintile" 
quintile_5groups <- g[,"quintile"] 

## # This picks the individuals with quintile=0 (note the
# double equals sign) and sets them to 5
quintile_5groups[quintile_5groups==0] <- 5 

## # lastly, tell R that this is a categorical variable
# and not a continuous one
quintile_5groups <- factor(quintile_5groups)

###################################################################
############################# WEEK 4 #############################
##### Checking the proportionality assumption #####
## The default arguments are the following ##
#cox.zph(fit, transform="km", global=TRUE) #

fit <- coxph(Surv(fu_time, death) ~ g$gender) # fit the desired model

temp <- cox.zph(fit)# apply the cox.zph function to the desired model

print(temp) # display the results

plot(temp) # plot the curves

### Include the interaction term and test whether it is
## statistically significant.

# "tt" is the time-transform function 
fit2 <- coxph(Surv(fu_time, death) ~ g$gender + tt(g$gender)) 
summary(fit2)

## Codice del corso ##
# fit <- coxph(Surv(fu_time, death) ~ gender)  
#temp <- cox.zph(fit)  
#print(temp)  
#plot(temp)

##### Using the other types of residuals in Cox regression #####
res.cox <- coxph(Surv(fu_time, death) ~ g$age) 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 


### Martingale residual is used to test the assumption
# that the countinuous variables we assumed to have a linear relation
# with the outcome actually do have a linear relation. 
ggcoxfunctional(Surv(fu_time, death) ~ g$age + log(g$age) + sqrt(g$age)) 

######################################################################
##### Results of the exercise on model selection and backwards elimination #####
# make the other covariates 

ihd <- factor(g[,'ihd']) 

valvular <- factor(g[,'valvular_disease']) 

pvd <- factor(g[,'pvd']) 

stroke <- factor(g[,'stroke']) 

copd<- factor(g[,'copd'])

pneumonia <- factor(g[,'pneumonia']) 

ht <- factor(g[,'hypertension'])

renal <- factor(g[,'renal_disease']) 

ca <- factor(g[,'cancer']) 

mets <- factor(g[,'metastatic_cancer']) 

mental_health <- factor(g[,'mental_health']) 

los <- g[,'los']

prior_dna <- g[,'prior_dnas']

# generate cognitive impairment variable (senility and dementia combined)

cog_imp <- as.factor(ifelse(g$dementia == 1 | g$senile == 1, 1, 0))

# run the full model 

cox <- coxph(Surv(fu_time, death) ~ g$age + gender + ethnicgroup + ihd + 
               
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               
               ca + mets + mental_health + cog_imp + los + prior_dna) 

summary(cox)

## Backwards elimination ##
cox2 <- coxph(Surv(fu_time, death) ~ g$age + gender + valvular + pneumonia + mets + cog_imp) 

summary(cox2)

table(cog_imp)
t <- table(cog_imp,death) 
t
round(100*prop.table(t,1),digits=1) 

### Testing the proportionality assumption on the remaining variables ###
fit <- coxph(Surv(fu_time, death) ~ g$age + gender + valvular + pneumonia + 
               
               mets + cog_imp) # test them all in the same model 

temp <- cox.zph(fit)  

print(temp) 
