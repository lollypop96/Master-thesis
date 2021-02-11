# Install and load relevant packages
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
library(splines)
library(lattice)
library(JM)
library(cmprsk)

# Load dataset
breast <- read.csv("C:/Users/liry9/Desktop/Tesi/Dataset/breast_cancer_2010_2015.txt", header=TRUE, sep = "\t", dec = ".")
View(breast)
sum(is.na(breast)) ## 0 nan values

## We Drop the "Behavior code ICD.o.3" variable since it does not provide any useful info
drops <- c("Behavior.code.ICD.O.3")
breast <- breast[ , !(names(breast) %in% drops)]
rm(drops)

## We rename the variables in order for them to be more manageable
colnames(breast)[1] <- "patient_id"
colnames(breast)[2] <- "ethnicity"
colnames(breast)[3] <- "year_birth"
colnames(breast)[4] <- "year_diagnosis"
colnames(breast)[5] <- "state_county"
colnames(breast)[6] <- "grade"
colnames(breast)[7] <- "laterality"
colnames(breast)[8] <- "hist_recode" ## histology recode
colnames(breast)[9] <- "sequence_number" 
colnames(breast)[10] <- "number_of_malignant_tumors"
colnames(breast)[11] <- "number_of_benign_tumors"
colnames(breast)[12] <- "seer_other_cause_of_death"
colnames(breast)[13] <- "survival_months"
colnames(breast)[14] <- "cause_of_death"
colnames(breast)[15] <- "surgery_of_primary_site"

### Create a new column for the Cause of death variable where:
# 0 = alive
# 1 = death
## Then we categorize the 
breast$status <- ifelse(breast$cause_of_death=="Alive", 0, 1)

## We modify the "surgery of primary site" variable
#breast$surgery_of_primary_site2 <- with(breast, ifelse(surgery_primary_site== "0", "None", ifelse(surgery_primary_site>="10"|surgery_primary_site<="19"),"external","internal"))

breast$surgery_of_primary_site2 <- 'None' # default
breast$surgery_of_primary_site2 <- with(breast, replace(surgery_of_primary_site2, surgery_of_primary_site >= 10 &  surgery_of_primary_site <= 19, 'Tumor destruction'))
breast$surgery_of_primary_site2 <- with(breast, replace(surgery_of_primary_site2, surgery_of_primary_site >= 20 & surgery_of_primary_site <= 80, 'Resection'))
breast$surgery_of_primary_site2 <- with(breast, replace(surgery_of_primary_site2, surgery_of_primary_site == 90, 'Surgery but no info'))
breast$surgery_of_primary_site2 <- with(breast, replace(surgery_of_primary_site2, surgery_of_primary_site == 98, 'Resection'))
breast$surgery_of_primary_site2 <- with(breast, replace(surgery_of_primary_site2, surgery_of_primary_site >= 99, 'Ill-defined'))

breast$survival_months <- as.numeric(breast$survival_months)
na_values <- breast[rowSums(is.na(breast)) > 0,]
breast <- breast[complete.cases(breast), ]

## Define the variables
# Define variables 
survival_time <- as.numeric(breast[,"survival_months"])
#survival_time_days <- as.numeric(breast[,"survival_days"])
death <- breast[,"status"] 
breast$age <- with(breast, (year_diagnosis-year_birth))
age <- breast[,"age"] 
ethnicgroup <- factor(breast[,"ethnicity"]) 
state_county <- factor(breast[,"state_county"]) 
grade <- factor(breast[,'grade']) 
laterality <- factor(breast[,'laterality']) 
hist_recode <- factor(breast[,'hist_recode']) 
seq_number <- factor(breast[,'sequence_number']) 
num_malignant_tumors <- breast[,'number_of_malignant_tumors'] 
num_benign_tumors <- breast[,'number_of_benign_tumors']
seer_cause_of_death <- factor(breast[,'seer_other_cause_of_death']) 
cause_of_death <- factor(breast[,'cause_of_death']) 
surgery_primary_site2 <- factor(breast[,'surgery_of_primary_site'])

### We have the survival time expressed in moths. Let's assume
# that each month has 30 days (in absence of any other information).
# We will therefore add a column called "survival_days" which is the
# result of the product of the values in the "survival_month" column
# times 30.
breast$survival_days <- breast$survival_months*30
#######################################################################
##################### Plotting a Kaplan-Meier curve ################### 
#######################################################################
### 1.a Generate the survival curve, for the survival_months 
km_fit_months <- survfit(Surv(survival_time, death) ~ 1, data=breast)
plot(km_fit_months, ylim=c(0.72,1))

### 1.b Generate the survival curve, for the survival_months 
km_fit_days <- survfit(Surv(survival_days, death) ~ 1, data=breast)
plot(km_fit_days, ylim=c(0.65,1))

### 1.c Since the plot is not exactly what we would expect we try with
# the following: we select a small sample of our dataset, in particular
# 10 rows and repeat the exact procedure.
### Mini Example ###
popol_example <- breast[sample(nrow(breast),10),]

## KM curve
fit_example <- survfit(Surv(survival_days, status) ~ 1, data=popol_example)
plot(fit_example)
text(fit_example$time, 0, format(fit_example$n.risk), cex = 0.7 )
summary(fit_example)

# If we do not want the confidence interval we can write:
plot(fit_example, conf.int=F)
text(fit_example$time, 0, format(fit_example$n.risk), cex = 0.7 )
summary(fit_example)
#####################

# The basic output of function survfit() gives the number of subjects,
# the number of events, the median survival time and its 95% confidence
# interval.
km_fit_days
### 2. Let's output the probability of surviving after a certain amount
# of days
summary(km_fit_days, times = c(500, 1000, 1500, 2000))

### 3. The quantile() method computes instead the corresponding follow-up
# times at which the survival probability takes a specific value. 
quantile(km_fit_days, probs = 1 - c(0.5, 0.6))
# Na values -> which is coherent with our graph, whose lower value
# is at circa 0.75

quantile(km_fit_days, probs = 1 - c(0.7, 0.75))

## In the probs argument of quantile() we have to specify
# one minus our target survival probabilities; this is because the
# function works under the cumulative distribution function (CDF)
# convention, and the CDF equals one minus survival probability.

### The plot() method produces the figure of the estimated curve;
# by default, the 95% confidence interval is included when we have
# only one curve:
plot(km_fit_days, xlab = "Time to Death (days)", ylab = "Survival Probability", 
     main = "Kaplan-Meier Estimate of S(t) for the Breast cancer dataset", ylim=c(0.7,1))

## Interpretation of the output:
# The survival columns represents the survival probability (the y axis of the previous plot)
# For example, let's deduce from this that the cumulative probability of survival at the end
# of 900 days is 78.7%.
# Since at the end of the time period considered there are still 6 patients alive, then we
# can deduct that the survival probability is not equal to zero.

#######################################################################################################
### Breslow Estimator ###
br_fit <- survfit(Surv(survival_days, status) ~ 1, data = breast, 
                            type = "fleming-harrington")
br_fit

plot(br_fit, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     ylim=c(0.7,1),
     main = "Breslow Estimate of S(t) for the SEER Breast Cancer Data")

## We try the same procedure for our small subset ##
br_fit_example <- survfit(Surv(survival_days, status) ~ 1, data = popol_example, 
                  type = "fleming-harrington")
br_fit_example

plot(br_fit_example, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     main = "Breslow Estimate of S(t) for the Example dataset")

#############################################################################
######################### Statistical tests #################################
#############################################################################
#### Log-Rank test ####
logrank_breast_ethnicity <- survdiff(Surv(survival_days, status == 1) ~ ethnicity, data = breast)
logrank_breast_ethnicity

# There is an important difference between the caucasian ethnicity and
# the black and the other ethnicities.

#### Peto & Peto Gehan-Wilcoxon test ####
peto_peto_breast <- survdiff(Surv(survival_days, status == 1) ~ ethnicity, data = breast, rho = 1)
peto_peto_breast

# The conclusion remains the same with also the same p-value.

#############################################################################
############# Checking for Proportional-Hazard assumption ###################
km_fit_breast_pha <- survfit(Surv(survival_days, status == 1) ~ ethnicity, data = breast)
plot(km_fit_breast_pha, fun = "cumhaz")
