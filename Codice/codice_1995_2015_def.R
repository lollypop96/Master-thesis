## We split our dataset into three subsample: train, validation and
# test sample
#write.csv(breast2,"C:/Users/liry9/Desktop/breast_py.csv", row.names = FALSE)
set.seed(127)

## We first drop the Surgery of primary site variable
#breast2 <- breast2[ , -which(names(breast2) %in% c("surgery_of_primary_site_labeled"))]

### Drop unnecessary columns
breast2 <- breast2[ , -which(names(breast2) %in% c('patient_id',
                                                   'type_follow_up_expected',
                                                   'sequence_number',
                                                   'primary_international',
                                                   'breast_ajcc',
                                                   'breast_ajcc_T',
                                                   'breast_ajcc_N',
                                                   'breast_ajcc_M',
                                                   'seer_hist_stage',
                                                   'cod',
                                                   'state_county',
                                                   'marital_status',
                                                   'hist_recode',
                                                   'year_birth',
                                                   'year_diagnosis',
                                                   'month_diagnosis',
                                                   'cause_of_death'))]

breast2 <- breast2[!(breast2$laterality == "Bilateral, single primary"),]
breast2 <- breast2[!(breast2$laterality == "Only one side - side unspecified"),]
breast2 <- breast2[!(breast2$primary_site_labeled == "C50.9-Breast, NOS"),]

spec = c(train = .7, testing = .15, validate = .15)

g_20y = sample(cut(
  seq(nrow(breast2)), 
  nrow(breast2)*cumsum(c(0,spec)),
  labels = names(spec)
))

res_20y = split(breast2, g_20y)

training_20y <- as.data.frame(res_20y$train)
test_20y <- as.data.frame(res_20y$testing)
valid_20y <- as.data.frame(res_20y$validate)

#######################################################################
##################### Plotting a Kaplan-Meier curve ################### 
#######################################################################
### 1.a Generate the survival curve, for the survival_months
km_fit_months <- survfit(Surv(survival_months, status) ~ 1,
                         data=training_20y)
plot(km_fit_months, ylim=c(0,1))
plot(km_fit_months, ylim=c(0,1),
     mark.time = FALSE,
     xlab = "Months", 
     ylab = "Overall survival probability")
summary(km_fit_months)
km_fit_months

survfit(Surv(survival_months, status) ~ 0,
                         data=training_20y)

ggsurvplot(km_fit_months,
           conf.int=TRUE,
           pval=TRUE,
           #risk.table=TRUE, 
           #legend.labs=c("Grade II","Grade III","Grade IV","Grade I"),
           #legend.title="Tumor Grade:",  
           palette=c("orchid"), 
           title="Kaplan-Meier Curve", 
           risk.table.height=.4,
           ylim=c(0.0,1))

ggsurvplot(fit = survfit(Surv(survival_months, status) ~ 1,
                data=training_20y), 
  palette=c("orchid"),
  title="Kaplan-Meier fit for Breast cancer",
  xlab = "Months", 
  ylab = "Overall survival probability")


# Mean and median survival time
print(km_fit_months, print.rmean=TRUE)

## Computing of confidence intervals based on 
# the complementary log-log transformation
km_fit_months_log <- survfit(Surv(survival_months, status) ~ 0,
                             conf.type="log-log",
                             data=training_20y)
plot(km_fit_months_log, ylim=c(0,1),
     mark.time = TRUE,
     xlab = "Months", 
     ylab = "Overall survival probability")
summary(km_fit_months_log)
km_fit_months_log

# Mean and median survival time
print(km_fit_months_log, print.rmean=TRUE)

### Kaplan-Meier in months for non censored patients
km_fit_days_non_censored <- survfit(Surv(survival_days, status) ~ 1,
                                conf.type="log-log",
                                data=training_20y)
plot(km_fit_days_non_censored, ylim=c(0,1),
     mark.time = TRUE,
     xlab = "Days", 
     ylab = "Overall survival probability")
summary(km_fit_days_non_censored)
km_fit_days_non_censored

# Mean and median survival time
print(km_fit_days_non_censored, print.rmean=TRUE)

### Kaplan Meier in days for censored patients
km_fit_days_censored <- survfit(Surv(survival_days, status) ~ 0,
                                  conf.type="log-log",
                                  data=training_20y)
plot(km_fit_days_censored, ylim=c(0,1),
     mark.time = TRUE,
     xlab = "Days", 
     ylab = "Overall survival probability")
summary(km_fit_days_censored)
km_fit_days_censored

# Mean and median survival time
print(km_fit_days_censored, print.rmean=TRUE)

### Nelson-Altschuler = alternative survival estimator
km_fit_months_na <- survfit(Surv(survival_months, status) ~ 1,
                            conf.type="log-log",
                            type="fh",
                            data=training_20y)
summary(km_fit_months_na)

# Mean and median survival time
print(km_fit_months_na, print.rmean=TRUE)

### Mean and median survival time (1st method)
# Median survival time
status <- training_20y$status
survival_months <- training_20y$survival_months
status_followup <- 1 - status
survfit(Surv(survival_months, status_followup) ~ 1)

median(survival_months) # 83 months

delta.followup <- 1-training_20y$status
survfit(Surv(survival_months, delta.followup) ~ 1, data=training_20y)

## Mean and median survival time (2nd method)
print(km_fit_months, print.rmean=TRUE)

summary(km_fit_months, times=seq(min(training_20y$survival_months), max(training_20y$survival_months), 60))

## Median survival time
#lung %>% 
#  filter(status == 2) %>% 
#  summarize(median_surv = median(time))

training_20y %>% 
  filter(training_20y$status == 1) %>% 
  summarize(median_surv = median(training_20y$survival_days))

# https://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Advanced-statistics/SurvivalAnalysis.html
km_fit_prodlim <- prodlim(Hist(survival_days,status)~1, data=training_20y)
plot(km_fit_prodlim, ylim=c(0,1),
     mark.time = TRUE,
     xlab = "Days", 
     ylab = "Overall survival probability")
summary(km_fit_prodlim)

par(mar=c(7,7,5,5),                          # margin of figure
    mgp=c(4,1,0))                            # move axis label away from figure
plot(km_fit_prodlim,
     xlab="Years",                           # label for x-axis
     ylab="Absolute risk of Death",          # label for y-axis
     type="cuminc",                          # increasing risks = 1-survival instead of decreasing survival
     axis1.at=seq(0,2900,365.25),            # time grid for x-axis
     axis1.labels=0:7,                       # time labels for x-axis
     axis2.las=2,                            # rotate labels of y-axis
     atrisk.dist=1,                          # adjust numbers below the figure
     atrisk.labels="Number of \npatients: ") # labels for numbers below figure

publish(km_fit_prodlim, times=seq(0,2900,365.25), org=TRUE) # tavola attuariale

# Tavola attuariale ho deciso io una divisione degli intervalli di tempo
# il k-m gli intervalli se li fa da solo -> crea tanti intervalli quanti sono gli eventi
# La tavola attuariale fa degli assunti un po' diversi rispetto alla km.
# Sono due approcci diversi

### 1.b Generate the survival curve, for the survival_days 
km_fit_days <- survfit(Surv(survival_days, status) ~ 1, data=training_20y)
plot(km_fit_days,
     ylim=c(0.3,1),
     xlab = "Days",
     ylab = "Overall survival probability",
     main="K-M Breast cancer (2010-2015)")
summary(km_fit_days)
km_fit_days

## Mean survival time
print(km_fit_days, print.rmean=TRUE)

## Median survival time
training_20y %>% 
  filter(status == 1) %>% 
  summarize(median_surv = median(survival_days))

### 1.c:Estimating x-year survival
# 1-year survival probability:
summary(survfit(Surv(survival_days, status) ~ 1,
                data = training_20y),
        times = 365.25)

# 5-year survival probability -> 12*5 = 60 months
summary(survfit(Surv(survival_months, status) ~ 1,
                data = training_20y),
        times = 60)

# 10-year survival probability -> 12*10 = 1200 months
summary(survfit(Surv(survival_months, status) ~ 1, data = training_20y), times = 120)

### 1.d: Fit Kaplan Meyer for tumor grade
km_fit_grade <- prodlim(Hist(survival_days, status)~grade,
                        data=training_20y)
par(mar=c(7,7,5,5), mgp=c(3,1,0))

plot(km_fit_grade,
     atrisk.labels=paste(c("Well differentiated; Grade I","Moderately differentiated; Grade II","Poorly differentiated; Grade III","Undifferentiated; anaplastic; Grade IV"),": "),
     atrisk.title="",
     xlab="Years",  # label for x-axis
     axis1.at=seq(0,2900,365.25), # time grid for x-axis
     axis1.labels=0:7, # time labels for x-axis
     legend.x="bottomleft", # positition of legend
     legend.cex=0.6, # font size of legend
     legend.title="Tumor Grade\n", # 
     logrank=TRUE) # show log-rank p-value

### 1.e: Fit the Kaplan-Meyer model for the estrogen receptor
km_fit_estrogen <- prodlim(Hist(survival_days, status)~estrogen_receptor, data=training_20y)
par(mar=c(7,7,5,5), mgp=c(3,1,0))
plot(km_fit_estrogen,
     atrisk.labels=paste("ER: ",c("Borderline","Negative","Positive"),": "),
     atrisk.title="",
     xlab="Years",  # label for x-axis
     axis1.at=seq(0,2900,365.25), # time grid for x-axis
     axis1.labels=0:7, # time labels for x-axis
     legend.x="bottomleft", # positition of legend
     legend.cex=0.8, # font size of legend
     legend.title="Tumor Grade\n", # 
     logrank=TRUE) # show log-rank p-value

### 1.e: Smooth survival plot - quantile of survival
sm.options(
  list(
    xlab = "Age (years)",
    ylab = "Time to death (months)")
)

#sm.survival(
#  x = prova3$age,
#  y = prova3$survival_months,
#  status = prova3$status,
#  h = (1/4)*sd(prova3$age) / nrow(prova3)^(-1/4)
#)

### 1.g Fit km for the follwing variables:
# Tumor grade
fit_km_grade <- survfit(Surv(survival_days, status) ~ grade,
                        data=training_20y)
ggsurvplot(fit_km_grade, ylim=c(0,1))

ggsurvplot(fit_km_grade, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Grade II","Grade III","Grade IV","Grade I"), legend.title="Tumor Grade:",  
           palette=c("dodgerblue2", "orchid2", "lightgreen", "lightblue4"), 
           title="Kaplan-Meier Curve: Tumor Grade", 
           risk.table.height=.4,
           ylim=c(0.35,1))

## Log-rank test:
survdiff(Surv(survival_days, status) ~ grade, data=training_20y)

# Ethnicity
fit_km_etn <- survfit(Surv(survival_days, status) ~ ethnicity, data=training_20y)
ggsurvplot(fit_km_etn, ylim=c(0,1))

ggsurvplot(fit_km_etn, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Black","Other","White"), legend.title="Ethnicity:",  
           palette=c("dodgerblue2", "orchid2", "lightgreen"), 
           title="Kaplan-Meier Curve: Ethnicity", 
           risk.table.height=.4,
           ylim=c(0,1), 
           ggtheme = theme_bw())

## Log-rank test:
survdiff(Surv(survival_days, status) ~ ethnicity, data=training_20y)

# Surgery of primary_international site
#fit_km_surg <- survfit(Surv(survival_days, status) ~ surgery_of_primary_site_labeled, data=training_20y)
#ggsurvplot(fit_km_surg, ylim=c(0,1), legend = "right")

#ggsurvplot(fit_km_surg,
#           pval=TRUE,
#           risk.table=TRUE, 
#           #legend = "right",
#           title="Kaplan-Meier Curve: Primary by international rules",
#           risk.table.height=.2,
#           ylim=c(0,1), 
#           ggtheme = theme_bw())

## Log-rank test:
#survdiff(Surv(survival_days, status) ~ surgery_of_primary_site_labeled, data=training_20y)

# Estrogen receptor
fit_km_er <- survfit(Surv(survival_days, status) ~ estrogen_receptor, data=training_20y)
ggsurvplot(fit_km_er, ylim=c(0,1), legend = "right")
ggsurvplot(fit_km_er, conf.int=TRUE, pval=TRUE, risk.table=TRUE,
           title="Kaplan-Meier Curve: Estrogen Receptor",
           legend.labs=c("Borderline","Negative","Positive"),
           legend.title="ER result:",  
           palette=c("dodgerblue2", "orchid2", "lightgreen"),
           risk.table.height=.3,
           ylim=c(0.3,1), 
           ggtheme = theme_bw())

## Log-rank test:
survdiff(Surv(survival_days, status) ~ estrogen_receptor, data=training_20y)

# Progesteron receptor
fit_km_pr <- survfit(Surv(survival_days, status) ~ progesteron_receptor, data=training_20y)
ggsurvplot(fit_km_pr, ylim=c(0,1), legend = "right")
ggsurvplot(fit_km_pr, conf.int=TRUE, pval=TRUE, risk.table=TRUE,
           title="Kaplan-Meier Curve: Progesteron Receptor",
           legend.labs=c("Borderline","Negative","Positive"),
           legend.title="PR result:",  
           palette=c("dodgerblue2", "orchid2", "lightgreen"),
           risk.table.height=.3,
           ylim=c(0.3,1), 
           ggtheme = theme_bw())

## Log-rank test:
survdiff(Surv(survival_days, status) ~ progesteron_receptor, data=training_20y)

## TNM classification
# T
fit_km_T <- survfit(Surv(survival_days, status) ~ breast_ajcc_T, data=training_20y)
ggsurvplot(fit_km_T,
           legend="left",
           ylim=c(0.2,1),
           ggtheme = theme_bw())

# N
fit_km_N <- survfit(Surv(survival_days, status) ~ breast_ajcc_N, data=training_20y)
ggsurvplot(fit_km_N,
           legend="bottom",
           ylim=c(0.4,1),
           ggtheme = theme_bw())

# M
fit_km_M <- survfit(Surv(survival_days, status) ~ breast_ajcc_M, data=training_20y)
ggsurvplot(fit_km_M,
           legend="bottom",
           ylim=c(0.15,1), 
           ggtheme = theme_bw())

### 1.f Since the plot is not exactly what we would expect we try with
# the following: we select a small sample of our dataset, in particular
# 10 rows and repeat the exact procedure.
### Mini Example ###
popol_example <- training_20y[sample(nrow(training_20y),100),]

## KM curve
fit_example <- survfit(Surv(survival_months, status) ~ 1, data=popol_example)
plot(fit_example, ylim=c(0,1),
     xlab = "Days", 
     ylab = "Overall survival probability")
#text(fit_example$time, 0, format(fit_example$n.risk), cex = 0.7 )
summary(fit_example)

# If we do not want the confidence interval we can write:
plot(fit_example, conf.int=F, ylim=c(0,1))
#text(fit_example$time, 0, format(fit_example$n.risk), cex = 0.7 )
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
# is at circa 0.8

quantile(km_fit_days, probs = 1 - c(0.8, 0.9))

# https://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Advanced-statistics/SurvivalAnalysis.html

## In the probs argument of quantile() we have to specify
# one minus our target survival probabilities; this is because the
# function works under the cumulative distribution function (CDF)
# convention, and the CDF equals one minus survival probability.

### The plot() method produces the figure of the estimated curve;
# by default, the 95% confidence interval is included when we have
# only one curve:
plot(km_fit_days, xlab = "Time to Death (days)", ylab = "Survival Probability", 
     main = "Kaplan-Meier Estimate of S(t) for the training_20y cancer dataset",
     ylim=c(0,1))

### 4.a) Quantiles for the age variable
training_20y$fage <- with(training_20y, cut(age, quantile(age), include = TRUE))

## don't use `attach()`; use the `data` argument of model fitting routine
fit_quantile <- survfit(Surv(survival_days, status) ~ training_20y$fage, data = training_20y, type="kaplan-meier")
summary(fit_quantile)
fit_quantile

fit <- survfit(Surv(survival_days, status) ~ 1, data=training_20y)
quantile(fit)

## Comment and interprepation:
# The median survival time is defined to be the time at which the
# survival curve crosses 50% survival. If the curve doesn't cross
# 50% (because survival is greater than 50% at the last time point),
# then median survival is simply undefined.
# More precisely, it is greater than the last time point on your
# survival curve. The only way around this, is if you fit some kind
# of model so are willing to extrapolate beyond the time span you
# have data for.

### 4.b) Quantiles for the size variable
training_20y$fsize <- with(training_20y, cut(tumor_size, quantile(tumor_size), include = TRUE))

## don't use `attach()`; use the `data` argument of model fitting routine
fit_quantile_size <- survfit(Surv(survival_days, status) ~ training_20y$fsize, data = training_20y, type="kaplan-meier")
summary(fit_quantile_size)
fit_quantile_size

#######################################################################################################
####### Breslow Estimator ########
br_fit <- survfit(Surv(survival_days, status) ~ 1, data = training_20y, 
                  type = "fleming-harrington")
br_fit

plot(br_fit, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     ylim=c(0,1),
     main = "Breslow Estimate of S(t) for the SEER training_20y Cancer Data")

## We try the same procedure for our small subset ##
br_fit_example <- survfit(Surv(survival_months, status) ~ 1, data = popol_example, 
                          type = "fleming-harrington")
br_fit_example

plot(br_fit_example, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     main = "Breslow Estimate of S(t) for the Example dataset",
     ylim=c(0,1))

#############################################################################
######################### Statistical tests #################################
#############################################################################
##### In order to verify the equality or proportionality whether or not the survival curves of different
# categories of the varibale taken into account we can use three tests:
# - the log-rank test,
# - the Breslow test,
# - the Peto & Peto Gehan Wilcoxon test
## The first 2 test are distributed as a Chi Square witn (m-1) degrees of
# freedom where m in the number of different group categories

#### Log-Rank test ####
## The log-rank test is the most powerful test when the proportional hazards
# assumption is satisfied. To check this assumption we can plot the cumulative
# hazard functions for the 2 groups. When PH is satisfied the two curves will
# be proportional to each other

logrank_breast_ethnicity <- survdiff(Surv(survival_days, status == 1) ~ ethnicity,
                                     data = training_20y)
logrank_breast_ethnicity

# There is an important difference between the caucasian ethnicity and
# the black and the other ethnicities.

#### Peto & Peto Gehan-Wilcoxon test ####
peto_peto_breast <- survdiff(Surv(survival_days, status == 1) ~ ethnicity,
                             data = training_20y,
                             rho = 1)
peto_peto_breast

# The conclusion remains the same with also the same p-value.

#############################################################################
############# Checking for Proportional-Hazard assumption ###################
km_fit_breast_pha <- survfit(Surv(survival_days, status == 1) ~ ethnicity,
                             data = training_20y)
plot(km_fit_breast_pha, fun = "cumhaz")
####### dopo cox c'è un altro test: probabilmente applicandolo è poco robusto
# quindi applicandolo a numeri classici di solito non rifiutiamo il test di
# PH -> cox.zph residui di Shoenfeld 

#############################################################################
#################### Accelerated Failure Time Models ####################
#########################################################################
###### Model Fitting ######
### Accelerated failure models are a parametric model that provides an
# alternative to the commonly used proportional hazards models.
# The AFT model assumes that the effect of a covariate is to accelerate
# or decelerate the life course of a disease by some constant.
## Unlike the PH models the regression parameter estimates from the AFT
# models are robust to omitted covariates.

### We fit an AFT model assuming the Weibull distribution
# for the breast dataset. We control for primary_international surgery, age at diagnosis. The code is:
fit_weibull <- survreg(Surv(survival_months, status) ~ ethnicity +
                         grade +
                         estrogen_receptor + 
                         progesteron_receptor +
                         ns(age, 3) +
                         regional_nodes_positive +
                         regional_nodes_examined +
                         ns(tumor_size,3),
                         data = training_20y)
summary(fit_weibull)
coef(fit_weibull)
vcov(fit_weibull)
fit_weibull$scale

fit_km_2 <- survfit(Surv(survival_months, status) ~ ethnicity +
                         grade +
                         estrogen_receptor + 
                         progesteron_receptor +
                         ns(age, 3) +
                         regional_nodes_positive +
                         regional_nodes_examined +
                         ns(tumor_size,3),
                       data = training_20y)

### Fit the same model but with the exponential distribution,
# the code is:
fit_exp <- survreg(Surv(survival_months, status) ~ ethnicity +
                     grade +
                     estrogen_receptor + 
                     progesteron_receptor +
                     ns(age, 3) +
                     regional_nodes_positive +
                     regional_nodes_examined +
                     ns(tumor_size,3),
                     data = training_20y,
                     dist = 'exponential')
summary(fit_exp)
coef(fit_exp)
vcov(fit_exp)
fit_exp$scale

### To fit the same model but with the log-normal distribution,
# the code is:
fit_lnorm <- survreg(Surv(survival_months, status) ~ ethnicity +
                       grade +
                       estrogen_receptor + 
                       progesteron_receptor +
                       ns(age, 3) +
                       regional_nodes_positive +
                       regional_nodes_examined +
                       ns(tumor_size,3),
                       data = training_20y,
                       dist = "lognormal")
summary(fit_lnorm)
coef(fit_lnorm)
vcov(fit_lnorm)
fit_lnorm$scale

### To fit the same model but with the log-logistic distribution,
# the code is:
fit_llogis <- survreg(Surv(survival_months, status) ~ ethnicity +
                        grade +
                        estrogen_receptor + 
                        progesteron_receptor +
                        ns(age, 3) +
                        regional_nodes_positive +
                        regional_nodes_examined +
                        ns(tumor_size,3),
                        data = training_20y,
                        dist = "loglogistic")
summary(fit_llogis)
coef(fit_llogis)
vcov(fit_llogis)
fit_llogis$scale

### Log Likelihood value is a measure of goodness of fit for any model.
# The higher the value, the better is the model. We should remember
# that Log Likelihood can lie between -Inf to +Inf -> therefore, based on our results
# we will opt for the Weibull model.
## We look for confirmation by using the anova function:
# this allows us to compute analysis of variance (or deviance) tables
# for one or more fitted model objects.

anova(fit_weibull, fit_exp, fit_lnorm, fit_llogis)

### AIC and BIC
# http://www.sthda.com/english/articles/38-regression-model-validation/158-regression-model-accuracy-metrics-r-square-aic-bic-cp-and-more/
# We all know that the main criteria for assessing the quality of the fit
# are the R^2 and the RMSE, however they have the flow that by increasing
# the size of the sample we cannot worsen the model so that we would
# go normally with the larger model.
## Akaike Information Criterion (AIC):
# AIC = -2logL(Beta_estimated)+2p
## where p is the penalty, which is a function of the size of the model
# and in particular p is the number of Beta parameters in the model.
# The smaller the aic, the better.

## The Bayesian Information Criterion (BIC)
# It is similar to the AIC, but it has a larger penalty.
# It quantifies the trade off between a model which fits well and
# the number of model parameters.
# As per the AIC we will prefer a small BIC.

glance(fit_exp) %>%
  dplyr::select(AIC, BIC)

glance(fit_weibull) %>%
  dplyr::select(AIC, BIC)

glance(fit_lnorm) %>%
  dplyr::select(AIC, BIC)

glance(fit_llogis) %>%
  dplyr::select(AIC, BIC)

## Based on both the AIC and the BIC we would choose the log logistic model
## Effect Plots ######
### Effect plots come into use when we consider statistical models
# with complex interactions and non-linear terms.

## As an example, we fit a model to the dataset that contains
# the grade variable, the linear effect of age,
# the quadratic effect of age. We assume that the survival times follow
# the Weibull distribution.
# The code to fit the model is:

fit_ef1 <- survreg(Surv(survival_months, status) ~ (grade)*(age + I(age^2)), 
                   data = training_20y)

summary(fit_ef1)

#### expand grid (surgery type con le età) -> dovrei trovare una prediction
# a tutti i tempi per quei surgery type e per quelle età

### Create the dataset with expand.grid
grade2 <- factor(training_20y[,'grade'])
breast_grid <- with(training_20y, expand.grid(
  age = seq(18, 80, length.out = 10),
  grade = levels(grade2)
))

breast_grid

#require(rms)
#cfit <- cph(Surv(survival_months, status) ~ surgery_of_primary_site_labeled, data = training_20y, surv=TRUE)
#survest(cfit, newdata=expand.grid(x=levels(training_20y$surgery_of_primary_site_labeled)) , 
#        times=seq(30,75,by=10)
#)$surv

## - Function seq() creates a sequence of numbers, in our example from 30
# to 101 years old, with the length of the sequence equal to 10.
## - Function levels() extracts from the training_20y database the possible
# values of the surgery of primary_international site variables.

head(breast_grid, 10) # to print the first ten rows of this dataset.

## Note: In the definition of the data frame breast_grid, all variables that
# appear in the model need to have their values specified. Also, you
# need to respect the class of these variables, e.g., in the example
# above and in the original data frame pbc2.id, variable surgery of primary_international sity is a
# factor. This is why we define its levels in ND.

### Now the predictions (in the log survival times scale) and the corresponding
# standard errors.

breast_prediction <- predict(fit_ef1,
                             breast_grid,
                             se.fit = TRUE,
                             type = "lp")
summary(breast_prediction)
head(breast_prediction$fit)
length(breast_prediction$fit)
#breast_grid$pred <- breast_prediction$fit# predictions
breast_grid$pred <- breast_prediction[[1]]
breast_grid$se <- breast_prediction[[2]] # standard errors
breast_grid$lo <- breast_grid$pred - 1.96 * breast_grid$se # lower limit of a 95% confidence interval for the log survival times.
breast_grid$up <- breast_grid$pred + 1.96 * breast_grid$se # upper limit of a 95% confidence interval for the log survival times.

head(breast_grid, 50)

### Finally, we create the plot using function xyplot() from package lattice:
xyplot(pred + lo + up ~ age | (grade),
       data = breast_grid,
       type = "l", 
       col = "black", lwd = 2, xlab = "Age", 
       ylab = "Log Survival Time")

anova(fit_ef1, fit_weibull, fit_exp, fit_lnorm, fit_llogis)

### Residuals ###
### To confirm that our Weibull model is correct we could also use Residuals
# Next, we construct the residuals based on their definition

## First we fit the model of interest, which is the fit_weibull one
fit_weibull <- survreg(Surv(survival_months, status) ~ ethnicity +
                         grade +
                         estrogen_receptor + 
                         progesteron_receptor +
                         ns(age, 3) +
                         regional_nodes_positive +
                         regional_nodes_examined +
                         ns(tumor_size,3),
                         data = training_20y)
fitted_values <- fit_weibull$linear.predictors
residuals <- (log(fit_weibull$y[, 1]) - fitted_values) / fit_weibull$scale
#plot(residuals)

# We then compute the Kaplan-Meier estimate of these residuals, and
# we plot it.
residuals_KM <- survfit(Surv(residuals, status) ~ 1, data = training_20y)
plot(residuals_KM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
x_axis <- seq(min(residuals), max(residuals), length.out = 35)
y_axis <- exp(- exp(x_axis))
lines(x_axis, y_axis, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

### Let's try the same, for compariso, for AFT models assuming the log-normal
# (also the one to be used for AFT models assuming the Gaussian distribution)
## After fitting the model of interest
fitted_values_lnorm <- fit_lnorm$linear.predictors
residuals_lnorm <- (log(fit_lnorm$y[, 1]) - fitted_values_lnorm) / fit_lnorm$scale

# We then compute the Kaplan-Meier estimate of these residuals, and
# we plot it.
residuals_KM_lnorm <- survfit(Surv(residuals_lnorm, status) ~ 1, data = training_20y)
plot(residuals_KM_lnorm, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
x_axis_lnorm <- seq(min(residuals_lnorm), max(residuals_lnorm), length.out = 35)
y_axis_lnorm <- pnorm(x_axis_lnorm, lower.tail = FALSE)
lines(x_axis_lnorm, y_axis_lnorm, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

### For AFT models assuming the log-logistic or the logistic distribution
# we use instead:
fitted_values_llogis <- fit_llogis$linear.predictors
residuals_llogis <- (log(fit_llogis$y[, 1]) - fitted_values_llogis) / fit_llogis$scale

# We then compute the Kaplan-Meier estimate of these residuals, and
# we plot it.
residuals_KM_llogis <- survfit(Surv(residuals_llogis, status) ~ 1, data = training_20y)
plot(residuals_KM_llogis, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
x_axis_llogis <- seq(min(residuals_llogis), max(residuals_llogis), length.out = 35)
y_axis_llogis <- plogis(x_axis_llogis, lower.tail = FALSE)
lines(x_axis_llogis, y_axis_llogis, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

##########################################################################
#################### Cox Proportional Hazards Models #####################
##########################################################################
###### Model Fitting ######
### We fit the Cox model to our breast cancer dataset
fit_cox <- coxph(Surv(survival_months, status) ~ 1,
                 data = training_20y)
summary(fit_cox)

###### Effect plots #######
### We want to produce an effect plot for the following Cox model for
# the dataset that contains the surgery of primary_international site, the nonlinear effect of age
# and its interaction with sex, and the nonlinear effect of the surg. of. prim. site:
fit_cox_breast_ef2 <- coxph(Surv(survival_months, status) ~ ethnicity +
                              grade +
                              estrogen_receptor + 
                              progesteron_receptor +
                              ns(age, 3) +
                              regional_nodes_positive +
                              regional_nodes_examined +
                              ns(tumor_size,3),
                              data = training_20y)
fit_cox_breast_ef2

## Result: In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
# Loglik converged before variable  5 ; coefficient may be infinite. 
# Interpretation at link: https://stat.ethz.ch/pipermail/r-help/2008-September/174201.html

# Using function expand.grid() we build the dataset that gives
# the combinations of covariate values based on which the plot will
# be constructed.

min(training_20y$tumor_size)
max(training_20y$tumor_size)

breast_grid2 <- with(training_20y, expand.grid(
  age = seq(18, 80, length.out = 10),
  surgery_of_primary_site_labeled = levels(training_20y$surgery_of_primary_site_labeled)
))
head(breast_grid2)

## We will do the plot, for increasing age from 30 to 101 years old,
## and for surgery of primary_international site
## Next, we obtain the predictions from the model (these are in the
# log-hazard scale), and we calculate the 95% confidence interval. 
breast_predictions2 <- predict(fit_cox_breast_ef2,
                               newdata = breast_grid2,
                               type = "lp",
                               se.fit = TRUE)
breast_grid2$pred <- breast_predictions2[[1]]
breast_grid2$se <- breast_predictions2[[2]]
breast_grid2$lo <- breast_grid2$pred - 1.96 * breast_grid2$se
breast_grid2$up <- breast_grid2$pred + 1.96 * breast_grid2$se
### Invece di ph.karno metti ad esempio i quantili per
# - numero di linfonodi esaminati
# - numero di linfonodi positivi
# - numero di tumori maligni
# - numero di tumori benigni

###### Hypothesis testing ######
### We can use the likelihood ratio test by comparing the models
# under the null and alternative hypothesis. This is done using
# the anova() function. 
# To test for the effect of the surgery of primary_international site, we run
# the test under the null hypothesis of the surgery.
# The alternative hypothesis contains instead the surgery.

fit_null_cox <- coxph(Surv(survival_months, status) ~ age, data = training_20y)
summary(fit_cox_breast_ef2)

anova(fit_null_cox, fit_cox_breast_ef2)

## There is a confirmation that age is always significan

############## Con l'età (variabile continua) -> metti delle funzioni delle età (età più età al quadrato)
# splines!!!!
# termplot -> funzione per visualizzare le spline
# COME SPLINE METTI AGE ,3 O ,2

## We will continue our analysis with a small sample of the dataset
#breast4 <- sample_frac(breast3, size = 0.1, replace = FALSE)
#View(training_20y)

# Grade
fit_cox_grade <- coxph(Surv(survival_months, status) ~ grade,
                       data=training_20y)
summary(fit_cox_grade)

# Laterality
fit_cox_lat <- coxph(Surv(survival_months, status) ~ laterality,
                     data=training_20y)
summary(fit_cox_lat)

# ICD-O-3 Hist/behav, malignant
fit_cox_hist <- coxph(Surv(survival_months, status) ~hist_recode,
                      data=training_20y)
summary(fit_cox_hist)

# Regional nodes examined
fit_cox_reg_examined <- coxph(Surv(survival_months, status) ~ regional_nodes_examined,
                              data=training_20y)
summary(fit_cox_reg_examined)

# Regional nodes positive
fit_cox_reg_pos <- coxph(Surv(survival_months, status) ~ regional_nodes_positive,
                         data=training_20y)
summary(fit_cox_reg_pos)

# Estrogen receptor
fit_cox_es <- coxph(Surv(survival_months, status) ~ estrogen_receptor,
                    data=training_20y)
summary(fit_cox_es)

# Progesteron receptor
fit_cox_pr <- coxph(Surv(survival_months, status) ~ progesteron_receptor,
                    data=training_20y)
summary(fit_cox_pr)

# Breast AJCC classification T
fit_cox_T <- coxph(Surv(survival_months, status) ~ breast_ajcc_T,
                   data=training_20y)
summary(fit_cox_T)

# Breast AJCC classification N
fit_cox_N <- coxph(Surv(survival_months, status) ~ breast_ajcc_N,
                   data=training_20y)
summary(fit_cox_N)

# Breast AJCC classification M
fit_cox_M <- coxph(Surv(survival_months, status) ~ breast_ajcc_M,
                   data=training_20y)
summary(fit_cox_M)

# age+lymph_nodes+laterality+seer_summary_stage+hist_recode+breast_ajcc_N

### ATTENZIONE: 1) quanti eventi abbiamo = non possiamo avere più variabili che eventi| 2) 
# - The column coef in the output denotes the log hazard ratios
# - exp(coef) denotes the the corresponding hazard ratios.
# fai un elenco delle variabili e vedi le correlazioni!!!

anova(fit_cox, fit_cox_mult)
# link: https://stats.stackexchange.com/questions/66591/coxph-ran-out-of-iterations-and-did-not-converge

#####################################################################################################
###### Hypothesis testing ######
### When a single categorical covariate is included in the Cox model,
# the summary() method returns:
# - the Wald,
# - the score (in this case equivalent to the log-rank test),
# - the likelihood ratio tests.

fit_cox_surg <- coxph(Surv(survival_days, status) ~ surgery_of_primary_site_labeled, data = training_20y)
summary(fit_cox_surg)

anova(fit_cox, fit_cox_mult, fit_cox_surg)

fit_alt_spline <- coxph(Surv(years, status2) ~ drug + ns(age, 3), data = pbc2.id)

###################################################################################
### Drop the following covariates
#training_20y <- training_20y[ , -which(names(training_20y) %in% c('patient_id','primary_site_labeled','type_follow_up_expected','sequence_number','primary_international','breast_ajcc','breast_ajcc_T','breast_ajcc_N',
#                                                                  'breast_ajcc_M','seer_hist_stage','cod','survival_days','state_county','laterality','marital_status','hist_recode','year_birth','year_diagnosis','month_diagnosis','cause_of_death'))]

#training_20y <- training_20y[ , -which(names(training_20y) %in% c('fage','fsize'))]

##################################################################################################
###### Proportional Hazards Assumption ######
### Kaplan-Meier estimate in log log scale
km_fit_copy <- survfit(Surv(survival_months, status) ~ 1,
                       data = training_20y)
km_fit_copy
plot(km_fit_copy, xlab = "Months", ylab = "Survival")
#legend("topright", levels(training_20y$surgery_of_primary_site_labeled), lty = 1, col = 1:2, bty = "n")

plot(km_fit_copy, fun = function (s) -log(-log(s)), xlab = "Month", 
     ylab = "-log(- log(Survival))", col = 1:2)

### Schoenfeld residuals test
fit_breast_schoen <- coxph(Surv(survival_months, status) ~ ethnicity +
                             grade +
                             estrogen_receptor + 
                             progesteron_receptor +
                             ns(age, 3) +
                             regional_nodes_positive +
                             regional_nodes_examined +
                             ns(tumor_size,3),
                             data = training_20y)

# The function that calculates the Schoenfeld residuals is cox.zph().
# The two primary_international arguments of this function are the fitted Cox
# model and the transformation of time to be used.
# The code below does the calculations for the KM scale.

check_PH <- cox.zph(fit_breast_schoen, transform = "km")
check_PH

# It is more useful to inspect the PH assumption by plotting the
# Schoenfeld residuals graphically.
par(mfrow = c(1, 1))
plot(check_PH, var = 1)
abline(h = coef(fit_breast_schoen)[1], col = "blue", lwd = 2)

#plot(check_PH, var = 2)
#abline(h = coef(fit_breast_schoen)[2], col = "red", lwd = 2)

#######################################################################
##################### Extension of the Cox Model ######################
#######################################################################
###### Survival probabilities from Cox Models #######
### We would like to derive survival probabilities from the
# following Cox model from the following variables
## Let us first define the variables
survival_months <- as.numeric(training_20y$survival_months)
ethnicity <- factor(training_20y[,"ethnicity"]) 
grade <- factor(training_20y[,'grade'])
regional_nodes_examined <- as.numeric(training_20y[,'regional_nodes_examined'])
regional_nodes_positive <- as.numeric(training_20y[,'regional_nodes_positive'])
estrogen_receptor <- factor(training_20y[,'estrogen_receptor'])
progesteron_receptor <- factor(training_20y[,'progesteron_receptor'])
status <- factor(training_20y[,"status"])
age <- as.numeric(training_20y[,"age"])
tumor_size <- as.numeric(training_20y[,"tumor_size"])
class(regional_nodes_examined)

## Fitting the model
fit_breast_cox_lin_age <- coxph(Surv(survival_months, status) ~ ethnicity +
                                  grade +
                                  estrogen_receptor + 
                                  progesteron_receptor +
                                  age +
                                  regional_nodes_positive +
                                  regional_nodes_examined +
                                  pspline(tumor_size,3),
                                  data = training_20y) # ok

fit_breast_cox1 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     ns(age, 3) +
                                     ns(tumor_size,3),
                                     data = training_20y) # ok

fit_breast_cox4 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     ns(age, 3) +
                                     #breast_ajcc +
                                     regional_nodes_examined +
                                     regional_nodes_positive +
                                     ns(tumor_size,3),
                                     data = training_20y,
                                     x= TRUE) # ok

summary(fit_breast_cox4)

training_sample <- sample_frac(training_20y,0.20)

fit_breast_cox5 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     pspline(age,3) +
                                     regional_nodes_positive +
                                     regional_nodes_examined +
                                     pspline(tumor_size,3),
                                     data = training_20y,
                                     x = TRUE) # ok

summary(fit_breast_cox5)

fit_breast_cox6 <- coxph(Surv(survival_months, status) ~ ethnicity +
                           grade +
                           estrogen_receptor + 
                           progesteron_receptor +
                           ns(age, 3) +
                           #breast_ajcc +
                           regional_nodes_examined +
                           regional_nodes_positive +
                           ns(tumor_size,3),
                           data = training_20y,
                           x = TRUE) # ok

summary(fit_breast_cox6)

concordance(fit_breast_cox6, newdata = test_20y)

prop_hazards6 <- cox.zph(fit_breast_cox6)
print(prop_hazards6)

anova(fit_breast_cox5, fit_breast_cox_lin_age)
AIC(fit_breast_cox5, fit_breast_cox_lin_age)
BIC(fit_breast_cox5, fit_breast_cox_lin_age)

fit_breast_cox4
fit_breast_cox4$concordance[[6]]
fit_breast_cox4$loglik
perror4=pec(list(Cox=fit_breast_cox4),Hist(survival_months,status)~ethnicity +
              # grade +
              estrogen_receptor + 
              progesteron_receptor +
              ns(age, 3) +
              breast_ajcc +
              regional_nodes_examined +
              regional_nodes_positive +
              ns(tumor_size,3),
              data = training_20y,
              #x= TRUE
            )

brier4 <- pec(object=fit_breast_cox4,
               formula=Surv(survival_months,status) ~ ethnicity+estrogen_receptor+progesteron_receptor+ns(age, 3)+breast_ajcc+regional_nodes_examined+regional_nodes_positive+ns(tumor_size,3),
               data=training_20y,
               exact=TRUE,
               cens.model="marginal",
               splitMethod="none",
               B=0,
               verbose=TRUE)


print(brier4,times=seq(min(training_20y$survival_months),max(training_20y$survival_months),12))
summary(brier4)
brier4
plot(brier4,xlim=c(0,max(training_20y$survival_months)))

fit_breast_cox5
fit_breast_cox5$concordance[[6]]
fit_breast_cox5$loglik
perror5=pec(list(Cox=fit_breast_cox5),
              Hist(survival_months,status) ~ ethnicity +
              grade +
              estrogen_receptor + 
              progesteron_receptor +
              pspline(age,3) +
              regional_nodes_positive +
              regional_nodes_examined +
              pspline(tumor_size,3),
              data = training_20y,
              #x= TRUE
              )

brier5 <- pec(object=fit_breast_cox5,
              formula=Surv(survival_months,status) ~ ethnicity+grade+estrogen_receptor+progesteron_receptor+pspline(age,3)+regional_nodes_positive+regional_nodes_examined+pspline(tumor_size,3),
              data=training_20y,
              exact=TRUE,
              cens.model="marginal",
              splitMethod="none",
              B=0,
              verbose=TRUE)

print(brier5,times=seq(min(training_20y$survival_months),max(training_20y$survival_months),12))
brier5
plot(brier5,xlim=c(0,max(training_20y$survival_months)))

fit_breast_cox6
fit_breast_cox6$concordance[[6]]
fit_breast_cox6$loglik
perror6=pec(list(Cox=fit_breast_cox6),
              Hist(survival_months,status) ~ ethnicity+grade+estrogen_receptor+progesteron_receptor+ns(age, 3)+regional_nodes_examined+regional_nodes_positive+ns(tumor_size,3),
              data = training_20y,
              x = TRUE)

brier6 <- pec(object=fit_breast_cox6,
              formula=Surv(survival_months, status) ~ ethnicity+grade+estrogen_receptor+progesteron_receptor+ns(age, 3)+regional_nodes_examined+regional_nodes_positive+ns(tumor_size,3),
              data=training_20y,
              exact=TRUE,
              cens.model="marginal",
              splitMethod="none",
              B=0,
              verbose=TRUE)

print(brier6,times=seq(min(training_20y$survival_months),max(training_20y$survival_months),12))
brier6
plot(brier6,xlim=c(0,max(training_20y$survival_months)))

## Checking the linearity assumption
plot(predict(fit_breast_cox4), residuals(fit_breast_cox4, type="martingale"),
     xlab="fitted values", ylab="Martingale residuals",
     main="Residual plot", las=1)
# There does not seem to be linearity

## Checking the PH assumption
ph <- cox.zph(fit_breast_cox4)
ph
plot <- ggcoxzph(ph)
plot

ggcoxdiagnostics(fit_breast_cox4,
                 type="schoenfeld",
                 ox.scale="time")

ggforest(fit_breast_cox4, data=training_20y)
ggforest(fit_breast_cox5, data=training_20y)
ggforest(fit_breast_cox6, data=training_20y)

# Plot the baseline survival function
ggsurvplot(survfit(fit_breast_cox4, data = training_20y),
           palette = "#2E9FDF")

# split input and output
head(test_20y)
copy_test <- test_20y
testing_labels <- as.numeric(copy_test$survival_months)
drops <- c("survival_months","survival_days")
copy_test <- copy_test[ , !(names(copy_test) %in% drops)]

## Fitting the model on the test set:
pred_validation = predict(fit_breast_cox5, newdata = copy_test)
#pred_validation2 <- predict (fit_breast_cox_lin_age,newdata = test_20y)

model_evaluation <- cbind(testing_labels, pred_validation)
colnames(model_evaluation) <- c('Actual', 'Predicted')
model_evaluation <- as.data.frame(model_evaluation)

mse <- mean((model_evaluation$Actual - model_evaluation$Predicted)^2)
rmse <- sqrt(mse)

data.frame(
  R2 = R2(pred_validation, testing_labels),
  RMSE = RMSE(pred_validation, testing_labels),
  MAE = MAE(pred_validation, testing_labels)
)

length(pred_validation)
length(testing_labels)
pred_validation
testing_labels

## Plotting smooth terms
fit_breast_cox_spline <- coxph(Surv(survival_months, status) ~ ethnicity +
                           grade +
                           estrogen_receptor + 
                           progesteron_receptor +
                           pspline(age, df=3) +
                           regional_nodes_positive +
                           regional_nodes_examined +
                           ns(tumor_size,3),
                           data = training_20y) # ok
ptemp <- termplot(fit_breast_cox_spline, se=TRUE, plot=FALSE)
attributes(ptemp)

ageterm <- ptemp$age # this will be a data frame
center <- with(ageterm, y[x==50])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96), '*')
matplot(ageterm$x, exp(ytemp - center), log='y',
          type='l', lty=c(1,2,2), col=1,
          xlab="Age at diagnosis", ylab="Relative death rate")

#######################################################################
#######################################################################
#######################################################################
# Sul dataset di training_20y ancora sottocampionamento sempre più grande
# per vedere se la dimensione campionaria impatta sul risultato
train_sample1 <- sample_frac(training_20y, size = size, replace = FALSE)

hist(training_20y$regional_nodes_examined, exclude=NA)
hist(training_20y$regional_nodes_positive, exclude=NA)
table(regional_nodes_positive)

table(grade, seer_summary_stage)

boh <- sample_frac(training_20y, 0.1)
View(boh)
cox_mod_fit

sample = seq(0.2 , 1, by = 0.1)
concordance_vector <- c()
             
### Da fare singolarmente
     
for (n in sample){
  print(n)
  data1 <- sample_frac(training_20y,n)
  cox_mod_fit <- coxph(Surv(survival_months, status) ~ ethnicity +
                         grade +
                         estrogen_receptor + 
                         progesteron_receptor +
                         ns(age, 3) +
                         #breast_ajcc +
                         regional_nodes_examined +
                         regional_nodes_positive +
                         ns(tumor_size,3),
                         data = data1,
                         x = TRUE)
  concordance_vector[n*10-1] <- cox_mod_fit$concordance[[6]]
  print(cox_mod_fit$concordance[[6]])
  brier <- pec(object=cox_mod_fit,
                formula=Surv(survival_months, status) ~ ethnicity+
                                                        grade+
                                                        estrogen_receptor+
                                                        progesteron_receptor+
                                                        ns(age,3)+
                                                        regional_nodes_examined+
                                                        regional_nodes_positive+
                                                        ns(tumor_size,3),
                data=data1,
                exact=TRUE,
                cens.model="marginal",
                splitMethod="none",
                B=0,
                verbose=TRUE)
  print(brier6,times=seq(min(data1$survival_months),max(data1$survival_months),60))
  ipa <- IPA(cox_mod_fit, formula=Surv(survival_months,status!=0)~1,times=60)
  print(ipa)}

print(concordance_vector)
SD(concordance_vector, na.rm = TRUE)

set.seed(271)
# For 0,1 converged
cox_mod_fit <- coxph(Surv(survival_months, status) ~ ethnicity +
                       grade +
                       breast_ajcc +
                       estrogen_receptor + 
                       progesteron_receptor +
                       ns(age, 3) +
                       breast_ajcc +
                       regional_nodes_positive +
                       regional_nodes_examined +
                       ns(tumor_size,3),
                     data = sample_frac(training_20y,1))


#########################################################################
## In order to understand the importance

# Do the same for a randomSurvivalForest model
##' library(randomForestSRC)
##' rsfmodel <- rfsrc(Surv(time,status)~X1+X2,data=d)
##' predictSurvProb(rsfmodel,newdata=ndat,times=ttt)
##' 
##' ## Cox with ridge option
##' f1 <- coxph(Surv(time,status)~X1+X2,data=learndat,x=TRUE,y=TRUE)
##' f2 <- coxph(Surv(time,status)~ridge(X1)+ridge(X2),data=learndat,x=TRUE,y=TRUE)
##' plot(predictSurvProb(f1,newdata=valdat,times=10),
##'      pec:::predictSurvProb.coxph(f2,newdata=valdat,times=10),
##'      xlim=c(0,1),
##'      ylim=c(0,1),
##'      xlab="Unpenalized predicted survival chance at 10",
##'      ylab="Ridge predicted survival chance at 10")
##'      
##'      
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
### Alternative code in case of convergence:
# we change the variables regional_nodes_examined"
# and "regional_nodes_positive" into categorical ones
copy_20y <- data.frame(training_20y)
table(copy_20y$regional_nodes_examined)
copy_20y$reg_nodes_ex <- 'None' # default
copy_20y$reg_nodes_ex <- with(copy_20y, replace(reg_nodes_ex, regional_nodes_examined == 1, "One"))
copy_20y$reg_nodes_ex <- with(copy_20y, replace(reg_nodes_ex, regional_nodes_examined == 2, 'Two'))
copy_20y$reg_nodes_ex <- with(copy_20y, replace(reg_nodes_ex, regional_nodes_examined >= 3, "Equal or More than three"))
table(copy_20y$reg_nodes_ex, exclude=NULL)

# Same thing for the regiona nodes positive
table(copy_20y$regional_nodes_positive)
copy_20y$reg_nodes_pos <- 'None' # default
copy_20y$reg_nodes_pos <- with(copy_20y, replace(reg_nodes_pos, regional_nodes_positive == 1, "One"))
copy_20y$reg_nodes_pos <- with(copy_20y, replace(reg_nodes_pos, regional_nodes_positive == 2, 'Two'))
copy_20y$reg_nodes_pos <- with(copy_20y, replace(reg_nodes_pos, regional_nodes_positive >= 3, "Equal or More than Three"))
table(copy_20y$reg_nodes_pos, exclude=NULL)

copy_20y <- copy_20y[, -which(names(copy_20y) %in% c("regional_nodes_examined","regional_nodes_positive"))]

fit_breast_cox_categ <- coxph(Surv(survival_days, status) ~ ethnicity +
                                grade +
                                breast_ajcc +
                                estrogen_receptor + 
                                progesteron_receptor +
                                ns(age, 3) +
                                reg_nodes_pos +
                                reg_nodes_ex +
                                ns(tumor_size,3),
                              data = copy_20y)

fit_breast_cox_categ
fit_breast_cox_categ$concordance[[6]]
fit_breast_cox_categ$loglik
ph_cox_categ <- cox.zph(fit_breast_cox_categ)
print(ph_cox_categ)
plot(ph_cox_categ)
ggforest(fit_breast_cox_categ, data = copy_20y)

help(ris)

### Codice prof ###
### IPA ###
data(pbc)
pbc = na.omit(pbc)
pbc$status[pbc$status == 2] <- 1
pbctest = (1:NROW(pbc)) %in% sample(1:NROW(pbc), size=0.632*NROW(pbc))
pbclearn = pbc[pbctest,]
pbctest

train_ind <- sample(seq_len(NROW(pbc)), size = 0.632*NROW(pbc))

train_pbc <- pbc[train_ind, ]
test_pbc <- pbc[-train_ind, ]

cox1 = coxph(Surv(time,status!=0) ~ age+
                                    sex+
                                    log(bili)+
                                    log(albumin)+
                                    log(protime),
                                    data=train_pbc,
                                    x=TRUE)

print(brier,times=seq(min(pbclearn$time),max(pbc$time),365))
IPA(cox1,
    formula=Surv(time,status!=0)~1,
    times=1000)

IPA(cox1,
    formula=Surv(time,status!=0)~1,
    times=500, newdata=test_pbc)


IPA(fit_breast_cox6,
    formula=Surv(survival_months,status!=0)~1,
    times=60, newdata=test_20y)

IPA(fit_breast_cox6,
    formula=Surv(survival_months,status!=0)~1,
    times=60)

### Prova con il training_20y

## [18:05] Federico Ambrogi
predictEventProb.nnetcr <- function(object, newdata, times, ...){
  newdata <- setDT(newdata)
  tms <- all.vars(object$call$formula)[-c(1:2)]
  newdata <- newdata[, ..tms]
  newdataNN <-   do.call(rbind, replicate(length(times), newdata, simplify=FALSE))
  newdataNN <- data.frame(h=rep(times, each=nrow(newdata)), newdataNN)
  pred <- matrix(predict(object$fit, newdataNN), nrow(newdata), length(times))
  return(pred)
}

fit <- survfit(fit_breast_cox4, newdata = test_20y)
fit$logse

test_C <- survival::survConcordance(fit_breast_cox4, newdata = test_20y)
test_C

concordance(fit_breast_cox4, newdata = test_20y)

###################################################################
library(rms)
cox_rms <- rms::cph(Surv(survival_months, status) ~ ethnicity +
                      grade +
                      estrogen_receptor + 
                      progesteron_receptor +
                      ns(age, 3) +
                      regional_nodes_positive +
                      regional_nodes_examined +
                      ns(tumor_size,3),
                      data = training_20y,
                      #x = T, 
                      #y = T
                    )

actuals <- Surv(external_validation$survival_months, external_validation$status)
estimates <- survest(cox_rms, newdata = external_validation, times = 1000)$surv
test_C_v1 <- as.numeric(rcorr.cens(x = estimates, S = actuals)[1])
round(test_C_v1, 4) # 0.6977
