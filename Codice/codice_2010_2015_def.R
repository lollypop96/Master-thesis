# Install and load relevant packages
#library(installr)
#updateR()
#remotes::install_github("mlr-org/mlr3extralearners")
#install_pycox(pip = TRUE, install_torch = TRUE)
#install_keras(pip = TRUE, install_tensorflow = TRUE)
library(riskRegression)
library(neuralnet)
library(PAmeasures)
library(sm)
library(randomForestSRC)
library(survivalmodels)
library(xgboost)
library(glmnet)
library(caret)
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
library(modelr)
library(broom)
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
library(splines)
library(lattice)
library(JM)
library(cmprsk)
library(sqldf)
library(sm)
library(condsurv)
library(params)
library(mlr3) # https://mlr3book.mlr-org.com
library(data.table)
library(pec)
library(prodlim)
library(Publish)
library(Tmisc)
library(caret)
library(corrplot)
library(ROCR)
library(mlr3proba)
library(mlr3learners)
library(mlr3pipelines)
library(mlr3viz)
library(MASS)
library(rcompanion)
library(lsr)
library(vcd)
library(DescTools)
library(survminer)
library(survPen)
library(GGally)
library(asaur)
library(rms)
library(pec)

# Load dataset
breast <- read.csv("C:/Users/liry9/Desktop/Tesi/Dataset/def_breast_1995_2015.txt", header=TRUE, sep = "\t", dec = ".")
sum(is.na(breast)) ## 0 nan values

### Missing data
#gg_na(breast)

breast %>%
  summarize(across(everything(), ~sum(is.na(.))/n()))

breast <- breast %>% drop_na()
sum(is.na(breast))

## We keep only the patients diagnosed with breast cancer between 2010 and 2015
breast2 <- data.frame(breast)

#View(breast2)
## We rename the variables in order for them to be more manageable
# Patient id
colnames(breast2)[1] <- "patient_id"

# Ethnicity
colnames(breast2)[2] <- "ethnicity"
breast2$ethnicity[breast2$ethnicity == "Unknown"] <- NA
table(breast2$ethnicity, exclude=NULL)
barplot(table(breast2$ethnicity), xlab = "Ethnicity", ylab = "Count", 
        names.arg = c("Black", "Other", "White"))
breast2$ethnicity <- ifelse(breast2$ethnicity == 'Other (American Indian/AK Native, Asian/Pacific Islander)',
                             'Other', breast2$ethnicity)

# State-County
colnames(breast2)[3] <- "state_county"
table(breast2$state_county, exclude=NULL)
unique(breast2$state_county)

# Primary site labeled
colnames(breast2)[4] <- "primary_site_labeled"
breast2$primary_site_labeled[breast2$primary_site_labeled == "C50.9"] <- NA
table(breast2$primary_site_labeled, exclude=NULL)
barplot(table(breast2$primary_site_labeled), xlab = "Primary site", ylab = "Count", 
        names.arg = c("C50.0", "C50.1", "C50.2", "C50.3", "C50.4",
                      "C50.5", "C50.6", "C50.8", "C50.9"), exclude=NULL)

# Grade
colnames(breast2)[5] <- "grade"
breast2$grade[breast2$grade == "Unknown"] <- NA
barplot(table(breast2$grade, exclude=NA))
barplot(table(breast2$grade), xlab = "Grade", ylab = "Count", 
        names.arg = c("Grade II",
                      "Grade III",
                      "Grade IV ",
                      "Grade I"))

# Laterality
colnames(breast2)[6] <- "laterality"
breast2$laterality[breast2$laterality == "Paired site, but no information concerning laterality"] <- NA
table(breast2$laterality, exclude=NULL)

# Histology Recode
colnames(breast2)[7] <- "hist_recode"
table(breast2$hist_recode, exclude=NULL)

# Surgery of primary_international site
colnames(breast2)[8] <- "surgery_of_primary_site"
breast2$surgery_of_primary_site[breast2$surgery_of_primary_site == 99] <- NA
breast2$surgery_of_primary_site[breast2$surgery_of_primary_site == "Blank(s)"] <- NA
table(breast2$surgery_of_primary_site, exclude=NULL)

# Regional nodes examined (link: http://web2.facs.org/cstage0205/breast/Breast_gpa.html)
colnames(breast2)[9] <- "regional_nodes_examined" 
breast2$regional_nodes_examined[breast2$regional_nodes_examined %in% c(95:126)] <- NA
table(breast2$regional_nodes_examined, exclude=NULL)
breast2 <- transform(breast2, regional_nodes_examined = as.numeric(regional_nodes_examined))

# Positive regional nodes (link: http://web2.facs.org/cstage0205/breast/Breast_fab.html)
colnames(breast2)[10] <- "regional_nodes_positive"
breast2$regional_nodes_positive[breast2$regional_nodes_positive %in% c(95:126)] <- NA
table(breast2$regional_nodes_positive, exclude=NULL)
breast2 <- transform(breast2, regional_nodes_positive = as.numeric(regional_nodes_positive))

# Estrogen receptor
colnames(breast2)[11] <- "estrogen_receptor"
breast2$estrogen_receptor[breast2$estrogen_receptor %in% c("Unknown", "Not 1900+ Breast")] <- NA
table(breast2$estrogen_receptor, exclude=NULL)

# Progesteron receptor
colnames(breast2)[12] <- "progesteron_receptor"
breast2$progesteron_receptor[breast2$progesteron_receptor %in% c("Unknown", "Not 1900+ Breast")] <- NA
table(breast2$progesteron_receptor, exclude=NULL)

# Number of concordant pairs
length(which(breast2$estrogen_receptor == breast2$progesteron_receptor))

# Number of discordant pairs
length(which(breast2$estrogen_receptor != breast2$progesteron_receptor))

prova <- breast2 %>% dplyr::select(estrogen_receptor, progesteron_receptor)
final <- prova[complete.cases(prova), ]
dmy <- dummyVars(" ~ .", data = final)
trsf <- data.frame(predict(dmy, newdata = final))

chisq.test(final$estrogen_receptor, final$progesteron_receptor)
# Since we get a p-value of less than the significance level of
# 0.05, we can reject the null hypothesis and conclude that the
# two variables are, indeed, independent.

# Cause of death
colnames(breast2)[13] <- "cause_of_death"
table(breast2$cause_of_death, exclude=NULL)

# Survival months
colnames(breast2)[14] <- "survival_months"

# Type of follow up expected
colnames(breast2)[15] <- "type_follow_up_expected"
table(breast2$type_follow_up_expected, exclude=NULL)

# Tumor sequence number
colnames(breast2)[16] <- "sequence_number"
table(breast2$sequence_number, exclude=NULL)

# Primary by international rules
colnames(breast2)[17] <- "primary_international"
table(breast2$primary_international, exclude=NULL)

# Number of total in situ malignant tumors
colnames(breast2)[18] <- "number_of_insitu_malignant_tumors"
breast2$number_of_insitu_malignant_tumors[breast2$number_of_insitu_malignant_tumors == 99] <- NA
table(breast2$number_of_insitu_malignant_tumors, exclude=NULL)
breast2 <- transform(breast2, number_of_insitu_malignant_tumors = as.numeric(number_of_insitu_malignant_tumors))

# Number of total borderline benign tumors
colnames(breast2)[19] <- "number_of_benign_tumors"
breast2$number_of_benign_tumors[breast2$number_of_benign_tumors == 99] <- NA
table(breast2$number_of_benign_tumors, exclude=NULL)
breast2 <- transform(breast2, number_of_benign_tumors = as.numeric(number_of_benign_tumors))

# Year of birth
colnames(breast2)[20] <- "year_birth"

# Month of diagnosis
colnames(breast2)[21] <- "month_diagnosis"
breast2$month_diagnosis[breast2$month_diagnosis == "Blank(s)"] <- NA
table(breast2$month_diagnosis)

# Year of diagnosis
colnames(breast2)[22] <- "year_diagnosis"
table(breast2$year_diagnosis)

# Breast AJCC classification
colnames(breast2)[23] <- "breast_ajcc"
table(breast2$breast_ajcc, exclude=NULL)
breast2$breast_ajcc[breast2$breast_ajcc == "UNK Stage"] <- NA
table(breast2$breast_ajcc, exclude=NULL)

# Breast AJCC classification T
colnames(breast2)[24] <- "breast_ajcc_T"
table(breast2$breast_ajcc_T, exclude=NULL)

# Breast AJCC classification N
colnames(breast2)[25] <- "breast_ajcc_N"
table(breast2$breast_ajcc_N, exclude=NULL)

# Breast AJCC classification M
colnames(breast2)[26] <- "breast_ajcc_M"
table(breast2$breast_ajcc_M, exclude=NULL)

# Seer historic stage
colnames(breast2)[27] <- "seer_hist_stage"
table(breast2$seer_hist_stage, exclude=NULL)

# Marital status
colnames(breast2)[28] <- "marital_status"
table(breast2$marital_status, exclude=NULL)
#breast2$marital_status[breast2$marital_status == "Unknown"] <- NA
#table(breast2$marital_status, exclude=NULL)

# Tumor size (2004-2015)
colnames(breast2)[29] <- "tumor_size"
table(breast2$tumor_size, exclude=NULL)
breast2$tumor_size[breast2$tumor_size %in% c(888, 990, 996:999)] <- NA
breast2$tumor_size[breast2$tumor_size == 991] <- 10
breast2$tumor_size[breast2$tumor_size == 992] <- 20
breast2$tumor_size[breast2$tumor_size == 993] <- 30
breast2$tumor_size[breast2$tumor_size == 994] <- 40
breast2$tumor_size[breast2$tumor_size == 995] <- 50

# Tumor size (1988-2003)
colnames(breast2)[30] <- "tumor_size_temp"
table(breast2$tumor_size_temp, exclude=NULL)
breast2$tumor_size_temp[breast2$tumor_size_temp %in% c(001, 996:999)] <- NA

#table(is.na(breast2$tumor_size), is.na(breast2$tumor_size_temp),exclude=NULL)
#breast2$boh <- breast2$tumor_size
#breast2$boh[breast2$boh == 990 & !is.na(breast2$boh)] <- NA
#breast2$boh[breast2$boh %in% c(996:999) & !is.na(breast2$boh)] <- NA
#breast2$boh[breast2$boh == 991 & !is.na(breast2$boh)] <- 10
#breast2$boh[breast2$boh == 992 & !is.na(breast2$boh)] <- 20
#breast2$boh[breast2$boh == 993 & !is.na(breast2$boh)] <- 30
#breast2$boh[breast2$boh == 994 & !is.na(breast2$boh)] <- 40
#breast2$boh[breast2$boh == 995 & !is.na(breast2$boh)] <- 50
#table(breast2$boh, exclude=NULL)

#breast2$boh2 <- breast2$tumor_size_temp
#table(breast2$boh2, exclude=NULL)
#breast2$boh2[breast2$boh2 == 990 & !is.na(breast2$boh2)] <- NA
#breast2$boh2[breast2$boh2 %in% c(996:999) & !is.na(breast2$boh2)] <- NA
#table(breast2$boh2, exclude=NULL)

#breast2$boh3 <- ifelse(is.na(breast2$boh), breast2$boh2, breast2$boh)
#table(breast2$boh3, exclude = NULL)

## We replace the Blank(s) value of the "tumor size" variable with
# the values of the "tumor size temp" variable:
breast2$tumor_size <- with(breast2, ifelse(tumor_size == "Blank(s)", tumor_size_temp, tumor_size))
table(breast2$tumor_size, exclude=NULL)

## We can now drop the"tumor size temp" variable:
breast2 <- breast2[ , -which(names(breast2) %in% c("tumor_size_temp"))]

breast2$tumor_size <- as.numeric(breast2$tumor_size)

## We create the variable age as the difference between the year of diagnosis and that of birth
breast2$age <- with(breast2, (year_diagnosis-year_birth))

max(breast2$age, na.rm = TRUE)
75 %in% breast2$age

summary(breast2$age) # no NAs
hist(breast2$age)

### We modify the following variables: 
## "surgery of primary_international site" variable
# breast2$surgery_of_primary_site2 <- with(breast2, ifelse(surgery_primary_site== "0", "None", ifelse(surgery_primary_site>="10"|surgery_primary_site<="19"),"external","internal"))
breast2$surgery_of_primary_site_labeled <- 'None' # default
breast2$surgery_of_primary_site_labeled <- with(breast2, replace(surgery_of_primary_site_labeled, surgery_of_primary_site >= 10 &  surgery_of_primary_site <= 19, 'Tumor destruction'))
breast2$surgery_of_primary_site_labeled <- with(breast2, replace(surgery_of_primary_site_labeled, surgery_of_primary_site >= 20 & surgery_of_primary_site <= 80, 'Resection'))
breast2$surgery_of_primary_site_labeled <- with(breast2, replace(surgery_of_primary_site_labeled, surgery_of_primary_site == 90, 'Surgery but no info'))
breast2$surgery_of_primary_site_labeled <- with(breast2, replace(surgery_of_primary_site_labeled, surgery_of_primary_site == 98, 'Resection'))
breast2$surgery_of_primary_site_labeled <- with(breast2, replace(surgery_of_primary_site_labeled, surgery_of_primary_site >= 99, 'Ill-defined'))
table(breast2$surgery_of_primary_site_labeled, exclude=NULL)

# Breast AJCC classification
table(breast2$breast_ajcc, exclude=NULL)
breast2$breast_ajcc[breast2$breast_ajcc %in% c("IIA", "IIB")] <- "II"
breast2$breast_ajcc[breast2$breast_ajcc %in% c("IIIA", "IIIB", "IIIC", "IIINOS")] <- "III"

# ER
table(breast2$estrogen_receptor, exclude=NULL)
breast2$estrogen_receptor[breast2$estrogen_receptor %in% c("Positive")] <- "ER_Positive"
breast2$estrogen_receptor[breast2$estrogen_receptor %in% c("Negative")] <- "ER_Negative"
breast2$estrogen_receptor[breast2$estrogen_receptor %in% c("Borderline")] <- "ER_Borderline"

# PR
table(breast2$progesteron_receptor, exclude=NULL)
breast2$progesteron_receptor[breast2$progesteron_receptor %in% c("Positive")] <- "PR_Positive"
breast2$progesteron_receptor[breast2$progesteron_receptor %in% c("Negative")] <- "PR_Negative"
breast2$progesteron_receptor[breast2$progesteron_receptor %in% c("Borderline")] <- "PR_Borderline"

# Creation of a bivariate COD (cause of status variable)
breast2$cod <- 'Other' # default
breast2$cod <- with(breast2, replace(cod, cause_of_death == "Breast", 'Breast'))
table(breast2$cod, exclude=NULL)

### Create a new column for the Cause of death variable where:
# 0 = alive
# 1 = death
breast2$status <- ifelse(breast2$cause_of_death=="Alive", 0, 1)
barplot(table(breast2$status))

# Remove surgery_of_primary_site column
breast2 <- breast2[ , -which(names(breast2) %in% c("surgery_of_primary_site"))]

## Survival months are non numerical in the table so we need to tell
# R that these data are acutally numbers. Furthermore
# we remove those rows having survival_months = 0.
breast2$survival_months <- as.numeric(breast2$survival_months)
na_values <- breast2[rowSums(is.na(breast2)) > 0,]
breast2 <- breast2[complete.cases(breast2), ]

### Survival times equal to 0
0 %in% breast2$survival_months # -> will return TRUE

## We replace the 00 values with a value equal to 0,01
breast2$survival_months <- with(breast2, replace(survival_months, survival_months == 00, 0.01))

# We run again this line of code: 
0 %in% breast2$survival_months # -> will return FALSE

## Summary of the age variable
summary(breast2$age) # no NAs
hist(breast2$age)

as.numeric(breast2$age)
breast2<-breast2[!(breast2$age<=18),]
breast2<-breast2[!(breast2$age>=85),]

summary(breast2$age)
hist(breast2$age, main = "Age distribution", xlab = "Age", ylab = "Count") 

na_values <- breast2[rowSums(is.na(breast2)) > 0,]
breast2 <- breast2[complete.cases(breast2), ]
sum(is.na(breast2))

### If instead we wanted to remove altoghether the survival months equal to zero
# we could run the following code:
# sqldf("select COUNT(*)
#        from breast2
#        where survival_months = 0")

## Let's remove the rows having 0 as survival time
# breast2 <- sqldf("SELECT *
#                  FROM breast2
#                  WHERE survival_months <> 0")

# sqldf("SELECT COUNT(*)
#      FROM breast2
#      WHERE survival_months = 0")

#0 %in% breast2$survival_months # -> will now return FALSE

### We have the survival time expressed in moths. Let's assume
# that each month has 30 days (in absence of any other information).
# We will therefore add a column called "survival_days" which is the
# result of the product of the values in the "survival_month" column
# times 30.
breast2$survival_days <- as.numeric(breast2$survival_months*30.436875) ### termine di conversione 30.625 per i bisestili 

## Create column for state variable
breast2$state <- substr(breast2$state_county,1,2)
table(breast2$state)
barplot(table(breast2$state, exclude=NA), xlab = "State", ylab = "Count")

######################################################################
### Exploratory Data Analysis ###

## Plot grade per ethnicity ##
qplot(grade, data = breast2, facets = ethnicity ~ .) +
  coord_flip()

######################################################################################
######################################################################################
## Define the variables
#survival_days <- as.numeric(breast2[,"survival_days"])
#ethnicity <- factor(breast2[,"ethnicity"]) 
#state_county <- factor(breast2[,"state_county"]) 
#grade <- factor(breast2[,'grade']) 
#laterality <- factor(breast2[,'laterality'])
#hist_recode <- factor(breast2[,'hist_recode'])
#seer_summary_stage <- factor(breast2[,'seer_combined_summary_stage'])
#surgery_of_primary_site_labeled <- factor(breast2[,'surgery_of_primary_site_labeled'])
#regional_nodes_examined <- factor(breast2[,'regional_nodes_examined'])
#regional_nodes_positive <- factor(breast2[,'regional_nodes_positive'])
#estrogen_receptor <- factor(breast2[,'estrogen_receptor'])
#progesteron_receptor <- factor(breast2[,'progesteron_receptor'])
#status <- factor(breast2[,"status"])
#survival_months <- as.numeric(breast2[,"survival_months"])
#type_follow_up_expected <- factor(breast2[,'type_follow_up_expected'])
#sequence_number <- factor(breast2[,'sequence_number'])
#primary_international <- factor(breast2[,'primary_international'])
#number_of_insitu_malignant_tumors <- as.numeric(breast2[,'number_of_insitu_malignant_tumors'])
#number_of_benign_tumors <- as.numeric(breast2[,'number_of_benign_tumors'])
#marital_status <- factor(breast2[,'marital_status'])
#seer_hist_stage <- factor(breast2[,'seer_hist_stage'])
#breast_ajcc_T <- factor(breast2[,'breast_ajcc_T'])
#breast_ajcc_N <- factor(breast2[,'breast_ajcc_N'])
#breast_ajcc_M <- factor(breast2[,'breast_ajcc_M'])
#cod <- factor(breast2[,'cod'])
#age <- as.numeric(breast2[,"age"])

breast_5y <- breast2[breast2$year_diagnosis %in% c("2010","2011","2012","2013","2014","2015"),]
breast_10y <- breast2[breast2$year_diagnosis %in% c("2005","2006","2007","2008","2009","2010",
                                                    "2011","2012","2013","2014","2015"),]
breast_15y <- breast2[breast2$year_diagnosis %in% c("2000","2001","2002","2003","2004","2005",
                                                    "2006","2007","2008","2009","2010","2011",
                                                    "2012","2013","2014","2015"),]


##################################################################################
########################### 2010-2015 Analysis ###################################
## We split our dataset into three subsample: train, validation and
# test sample
set.seed(127)
spec = c(train = .6, testing = .2, validate = .2)

g = sample(cut(
  seq(nrow(breast_5y)), 
  nrow(breast_5y)*cumsum(c(0,spec)),
  labels = names(spec)
))

res = split(breast_5y, g)

training <- as.data.frame(res$train)
test <- as.data.frame(res$testing)
valid <- as.data.frame(res$validate)

#######################################################################
##################### Plotting a Kaplan-Meier curve ################### 
#######################################################################
### 1.a Generate the survival curve, for the survival_months
#km_fit_months <- survfit(Surv(survival_months, status) ~ 1, data=training)
km_fit_months <- survfit(Surv(survival_months, status) ~ 1, data=training)
plot(km_fit_months, ylim=c(0,1))
plot(km_fit_months, ylim=c(0.8,1),
     xlab = "Days", 
     ylab = "Overall survival probability")
summary(km_fit_months)

km_fit_months
## Mean survival time
print(km_fit_months, print.rmean=TRUE)

## Median survival time
training %>% 
  filter(status == 1) %>% 
  summarize(median_surv = median(survival_months))

# https://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Advanced-statistics/SurvivalAnalysis.html
km_fit_prodlim <- prodlim(Hist(survival_days,status)~1, data=training)
plot(km_fit_prodlim, ylim=c(0.7,1),
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
km_fit_days <- survfit(Surv(survival_days, status) ~ 1, data=training)
plot(km_fit_days,
     ylim=c(0.7,1),
     xlab = "Days",
     ylab = "Overall survival probability",
     main="K-M Breast cancer (2010-2015)")
summary(km_fit_days)

## Mean survival time
print(km_fit_days, print.rmean=TRUE)

## Median survival time
training %>% 
  filter(status == 1) %>% 
  summarize(median_surv = median(survival_days))

### 1.c: Fit Kaplan Meyer for tumor grade
km_fit_grade <- prodlim(Hist(survival_days, status)~grade, data=training)
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

### 1.d: Fit the Kaplan-Meyer model for the estrogen receptor
km_fit_estrogen <- prodlim(Hist(survival_days, status)~estrogen_receptor, data=training)
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

### 1.e Fit km for the follwing variables:
# Tumor grade
fit_km_grade <- survfit(Surv(survival_days, status) ~ grade, data=training)
ggsurvplot(fit_km_grade, ylim=c(0.6,1))

ggsurvplot(fit_km_grade, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Grade I","Grade II","Grade III","Grade IV"), legend.title="Grade",  
           palette=c("dodgerblue2", "orchid2", "deeppink", "lightblue4"), 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.4,
           ylim=c(0.45,1))

# Ethnicity
fit_km_etn <- survfit(Surv(survival_days, status) ~ ethnicity, data=training)
ggsurvplot(fit_km_etn, ylim=c(0.6,1))

ggsurvplot(fit_km_etn, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("White","Other","Black"), legend.title="Ethnicity",  
           palette=c("dodgerblue2", "orchid2", "deeppink"), 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.4,
           ylim=c(0.65,1), 
           ggtheme = theme_bw())

# Combined summary stage
fit_km_sum_stage <- survfit(Surv(survival_days, status) ~ seer_summary_stage, data=training)
ggsurvplot(fit_km_sum_stage,
           legend = "right", 
           ggtheme = theme_bw())

# Surgery of primary_international site
fit_km_surg <- survfit(Surv(survival_days, status) ~ surgery_of_primary_site_labeled, data=training)
ggsurvplot(fit_km_surg,
           legend = "right", 
           ggtheme = theme_bw())

# Estrogen receptor
fit_km_er <- survfit(Surv(survival_days, status) ~ estrogen_receptor, data=training)
ggsurvplot(fit_km_er,
           legend = "right", 
           ggtheme = theme_bw())

# Progesteron receptor
fit_km_pr <- survfit(Surv(survival_days, status) ~ progesteron_receptor, data=training)
ggsurvplot(fit_km_pr,
           legend = "right", 
           ggtheme = theme_bw())

## TNM classification
# T
fit_km_T <- survfit(Surv(survival_days, status) ~ breast_ajcc_T, data=training)
ggsurvplot(fit_km_T,
           legend="left",
           ylim=c(0.2,1),
           ggtheme = theme_bw())

# N
fit_km_N <- survfit(Surv(survival_days, status) ~ breast_ajcc_N, data=training)
ggsurvplot(fit_km_N,
           legend="bottom",
           ylim=c(0.4,1),
           ggtheme = theme_bw())

# M
fit_km_M <- survfit(Surv(survival_days, status) ~ breast_ajcc_M, data=training)
ggsurvplot(fit_km_M,
           legend="bottom",
           ylim=c(0.15,1), 
           ggtheme = theme_bw())

### 1.f Since the plot is not exactly what we would expect we try with
# the following: we select a small sample of our dataset, in particular
# 10 rows and repeat the exact procedure.
### Mini Example ###
popol_example <- training[sample(nrow(training),10),]

## KM curve
fit_example <- survfit(Surv(survival_months, status) ~ 1, data=popol_example)
plot(fit_example, ylim=c(0.7,1))
#text(fit_example$time, 0, format(fit_example$n.risk), cex = 0.7 )
summary(fit_example)

# If we do not want the confidence interval we can write:
plot(fit_example, conf.int=F, ylim=c(0.72,1))
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
     main = "Kaplan-Meier Estimate of S(t) for the training cancer dataset", ylim=c(0.7,1))

### 4.a) Quantiles for the age variable
training$fage <- with(training, cut(age, quantile(age), include = TRUE))

## don't use `attach()`; use the `data` argument of model fitting routine
fit_quantile <- survfit(Surv(survival_days, status) ~ training$fage, data = training, type="kaplan-meier")
summary(fit_quantile)
fit_quantile

fit <- survfit(Surv(survival_days, status) ~ 1, data=training)
quantile(fit)

cfit <- coxph(Surv(time, status) ~ age + strata(ph.ecog), data=lung)
csurv<- survfit(cfit, newdata=data.frame(age=c(40, 60, 80)),
                conf.type ="none")
temp <- quantile(csurv, 1:5/10)
temp[2,3,]  # quantiles for second level of ph.ecog, age=80
quantile(csurv[2,3], 1:5/10)  # quantiles of a single curve, same result

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
training$fsize <- with(training, cut(tumor_size, quantile(tumor_size), include = TRUE))

## don't use `attach()`; use the `data` argument of model fitting routine
fit_quantile_size <- survfit(Surv(survival_days, status) ~ training$fsize, data = training, type="kaplan-meier")
summary(fit_quantile_size)
fit_quantile_size

fit_quantile2 <- survfit(Surv(survival_months, status) ~ findInterval(training$age, quantile(training$age)[-5]), 
                         data=training, type = "kaplan-meier")
fit_quantile2

#######################################################################################################
####### Breslow Estimator ########
br_fit <- survfit(Surv(survival_days, status) ~ 1, data = training, 
                  type = "fleming-harrington")
br_fit

plot(br_fit, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     ylim=c(0.7,1),
     main = "Breslow Estimate of S(t) for the SEER training Cancer Data")

## We try the same procedure for our small subset ##
br_fit_example <- survfit(Surv(survival_months, status) ~ 1, data = popol_example, 
                          type = "fleming-harrington")
br_fit_example

plot(br_fit_example, mark.time = FALSE, xlab = "Time to Death (days)", 
     ylab = "Survival Probability",
     main = "Breslow Estimate of S(t) for the Example dataset",
     ylim=c(0.7,1))

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

logrank_breast_ethnicity <- survdiff(Surv(survival_days, status == 1) ~ ethnicity, data = training)
logrank_breast_ethnicity

# There is an important difference between the caucasian ethnicity and
# the black and the other ethnicities.

#### Peto & Peto Gehan-Wilcoxon test ####
peto_peto_breast <- survdiff(Surv(survival_days, status == 1) ~ ethnicity, data = training, rho = 1)
peto_peto_breast

# The conclusion remains the same with also the same p-value.

#############################################################################
############# Checking for Proportional-Hazard assumption ###################
km_fit_breast_pha <- survfit(Surv(survival_days, status == 1) ~ ethnicity, data = training)
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
fit_weibull <- survreg(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training)
summary(fit_weibull)
coef(fit_weibull)
vcov(fit_weibull)
fit_weibull$scale

### Fit the same model but with the exponential distribution,
# the code is:
fit_exp <- survreg(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training, dist = 'exponential')
summary(fit_exp)
coef(fit_exp)
vcov(fit_exp)
fit_exp$scale

### To fit the same model but with the log-normal distribution,
# the code is:
fit_lnorm <- survreg(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training, dist = "lognormal")
summary(fit_lnorm)
coef(fit_lnorm)
vcov(fit_lnorm)
fit_lnorm$scale

### To fit the same model but with the log-logistic distribution,
# the code is:
fit_llogis <- survreg(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training, dist = "loglogistic")
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
# the surgery of primary_international site variable, the linear effect of age,
# the quadratic effect of age. We assume that the survival times follow
# the Weibull distribution.
# The code to fit the model is:

fit_ef1 <- survreg(Surv(survival_months, status) ~ (surgery_of_primary_site_labeled)*(age + I(age^2)), 
                   data = training)

summary(fit_ef1)

#### expand grid (surgery type con le età) -> dovrei trovare una prediction
# a tutti i tempi per quei surgery type e per quelle età

### Create the dataset with expand.grid 
breast_grid <- with(training, expand.grid(
  age = seq(30, 101, length.out = 10),
  surgery_of_primary_site_labeled = levels(surgery_of_primary_site_labeled)
))

breast_grid

#require(rms)
#cfit <- cph(Surv(survival_months, status) ~ surgery_of_primary_site_labeled, data = training, surv=TRUE)
#survest(cfit, newdata=expand.grid(x=levels(training$surgery_of_primary_site_labeled)) , 
#        times=seq(30,75,by=10)
#)$surv

## - Function seq() creates a sequence of numbers, in our example from 30
# to 101 years old, with the length of the sequence equal to 10.
## - Function levels() extracts from the training database the possible
# values of the surgery of primary_international site variables.

head(breast_grid, 10) # to print the first ten rows of this dataset.

## Note: In the definition of the data frame breast_grid, all variables that
# appear in the model need to have their values specified. Also, you
# need to respect the class of these variables, e.g., in the example
# above and in the original data frame pbc2.id, variable surgery of primary_international sity is a
# factor. This is why we define its levels in ND.

### Now the predictions (in the log survival times scale) and the corresponding
# standard errors.

breast_prediction <- predict(fit_ef1, breast_grid, se.fit = TRUE, type = "lp")
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
xyplot(pred + lo + up ~ age | (surgery_of_primary_site_labeled), data = breast_grid, type = "l", 
       col = "black", lwd = 2, xlab = "Age", 
       ylab = "Log Survival Time")

anova(fit_ef1, fit_weibull, fit_exp, fit_lnorm, fit_llogis)

### Residuals ###
### To confirm that our Weibull model is correct we could also use Residuals
# Next, we construct the residuals based on their definition

## First we fit the model of interest, which is the fit_weibull one
fit_weibull <- survreg(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training)
fitted_values <- fit_weibull$linear.predictors
residuals <- (log(fit_weibull$y[, 1]) - fitted_values) / fit_weibull$scale

# We then compute the Kaplan-Meier estimate of these residuals, and
# we plot it.

residuals_KM <- survfit(Surv(residuals, status) ~ 1, data = training)
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

residuals_KM_lnorm <- survfit(Surv(residuals_lnorm, status) ~ 1, data = training)
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

residuals_KM_llogis <- survfit(Surv(residuals_llogis, status) ~ 1, data = training)
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
fit_cox <- coxph(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training)
summary(fit_cox)

###### Effect plots #######
### We want to produce an effect plot for the following Cox model for
# the dataset that contains the surgery of primary_international site, the nonlinear effect of age
# and its interaction with sex, and the nonlinear effect of the surg. of. prim. site:
fit_cox_breast_ef2 <- coxph(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + ns(age, 3), data = training)
fit_cox_breast_ef2

## Result: In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
# Loglik converged before variable  5 ; coefficient may be infinite. 
# Interpretation at link: https://stat.ethz.ch/pipermail/r-help/2008-September/174201.html

# Using function expand.grid() we build the dataset that gives
# the combinations of covariate values based on which the plot will
# be constructed.

breast_grid2 <- with(training, expand.grid(
  age = seq(30, 101, length.out = 10),
  surgery_of_primary_site_labeled= levels(surgery_of_primary_site_labeled)
))
head(breast_grid2)

## We will do the plot, for increasing age from 30 to 101 years old,
## and for surgery of primary_international site
## Next, we obtain the predictions from the model (these are in the
# log-hazard scale), and we calculate the 95% confidence interval. 
breast_predictions2 <- predict(fit_cox_breast_ef2, newdata = breast_grid2, type = "lp", se.fit = TRUE)
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

fit_null_cox <- coxph(Surv(survival_months, status) ~ age, data = training)
fit_alt_cox <- coxph(Surv(survival_months, status) ~ surgery_of_primary_site_labeled + age, data = training)
summary(fit_alt_cox)

anova(fit_null_cox, fit_alt_cox)

## There is a confirmation that age is always significan

############## Con l'età (variabile continua) -> metti delle funzioni delle età (età più età al quadrato)
# splines!!!!
# termplot -> funzione per visualizzare le spline
# COME SPLINE METTI AGE ,3 O ,2

# Remove patiend_id, cod, primary_interntional, surgery of primary_international site
# type of follow up expected, state_county, survival_days, seer_combined_summary_stage,
# cs_tumor size, cs_extension, cs_lymph nodes, sequence_number
drops <- c("patient_id","cod", "primary_international", "surgery_of_primary_site_labeled","survival_days","type_follow_up_expected", "state_county",
           "seer_combined_summary_stage", "cs_tumor_size", "cs_tumor_extension", "cs_lymph_nodes", "sequence_number")
breast3 <- training[ , !(names(training) %in% drops)]
## We will continue our analysis with a small sample of the dataset
#breast4 <- sample_frac(breast3, size = 0.1, replace = FALSE)
#View(training)

# Grade
fit_cox_grade <- coxph(Surv(survival_months, status) ~ grade, data=breast3)
summary(fit_cox_grade)

# Laterality
fit_cox_lat <- coxph(Surv(survival_months, status) ~ laterality, data=breast3)
summary(fit_cox_lat)

# ICD-O-3 Hist/behav, malignant
fit_cox_hist <- coxph(Surv(survival_months, status) ~ hist_recode, data=breast3)
summary(fit_cox_hist)

# Regional nodes examined
fit_cox_reg_examined <- coxph(Surv(survival_months, status) ~ regional_nodes_examined, data=breast3)
summary(fit_cox_reg_examined)

# Regional nodes positive
fit_cox_reg_pos <- coxph(Surv(survival_months, status) ~ regional_nodes_positive, data=breast3)
summary(fit_cox_reg_pos)

# Estrogen receptor
fit_cox_es <- coxph(Surv(survival_months, status) ~ estrogen_receptor, data=breast3)
summary(fit_cox_es)

# Progesteron receptor
fit_cox_pr <- coxph(Surv(survival_months, status) ~ progesteron_receptor, data=breast3)
summary(fit_cox_pr)

# Breast AJCC classification T
fit_cox_T <- coxph(Surv(survival_months, status) ~ breast_ajcc_T, data=breast3)
summary(fit_cox_T)

# Breast AJCC classification N
fit_cox_N <- coxph(Surv(survival_months, status) ~ breast_ajcc_N, data=breast3)
summary(fit_cox_N)

# Breast AJCC classification M
fit_cox_M <- coxph(Surv(survival_months, status) ~ breast_ajcc_M, data=breast3)
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

fit_cox_surg <- coxph(Surv(survival_days, status) ~ surgery_of_primary_site_labeled, data = training)
summary(fit_cox_surg)

anova(fit_cox, fit_cox_mult, fit_cox_surg)

fit_alt_spline <- coxph(Surv(years, status2) ~ drug + ns(age, 3), data = pbc2.id)
##################################################################################################
###### Proportional Hazards Assumption ######
### Kaplan-Meier estimate in log log scale
km_fit_copy <- survfit(Surv(survival_months, status) ~ surgery_of_primary_site_labeled, data = training)
km_fit_copy
plot(km_fit_copy, xlab = "Months", ylab = "Survival")
legend("topright", levels(training$surgery_of_primary_site_labeled), lty = 1, col = 1:2, bty = "n")

plot(km_fit_copy, fun = function (s) -log(-log(s)), xlab = "Month", 
     ylab = "-log(- log(Survival))", col = 1:2)

### Schoenfeld residuals test
fit_breast_schoen <- coxph(Surv(survival_months, status) ~ ethnicity+estrogen_receptor+progesteron_receptor+breast_ajcc_T, data = training)

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
fit_breast_cox_extension1 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     ns(age, 3) +
                                     ns(tumor_size,3), data = training) # converged

fit_breast_cox_extension2 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     seer_summary_stage +
                                     ns(age, 3) +
                                     ns(tumor_size,3), data = training) # converged

fit_breast_cox_extension3 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     seer_summary_stage +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     ns(age, 3) +
                                     breast_ajcc_T +
                                     breast_ajcc_N +
                                     ns(tumor_size,3), data = training) # converged

fit_breast_cox_extension4 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     seer_summary_stage +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     ns(age, 3) +
                                     breast_ajcc +
                                     regional_nodes_examined +
                                     ns(tumor_size,3), data = training) # not converged

fit_breast_cox_extension5 <- coxph(Surv(survival_months, status) ~ ethnicity +
                                     grade +
                                     seer_summary_stage +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     ns(age, 3) +
                                     breast_ajcc_T +
                                     breast_ajcc_N +
                                     regional_nodes_positive +
                                     ns(tumor_size,3), data = training) # not converged

fit_breast_cox_extension <- coxph(Surv(survival_months, status) ~ ethnicity +
                                    grade +
                                    seer_summary_stage +
                                    regional_nodes_examined +
                                    regional_nodes_positive +
                                    estrogen_receptor + 
                                    progesteron_receptor +
                                    ns(age, 3) +
                                    ns(tumor_size,3), data = training)
## Since they still converge we change the variables regional_nodes_examined"
# and "regional_nodes_positive" into categorical ones
table(training$regional_nodes_examined)
training$reg_nodes_ex <- 'None' # default
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined == 1, "One"))
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined == 2, 'Two'))
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined == 3, "Three"))
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined == 4, 'Four'))
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined == 5, "Five"))
training$reg_nodes_ex <- with(training, replace(reg_nodes_ex, regional_nodes_examined > 5, 'More than five'))
table(training$reg_nodes_ex, exclude=NULL)

# Same thing for the regiona nodes positive
table(training$regional_nodes_positive)
training$reg_nodes_pos <- 'None' # default
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive == 1, "One"))
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive == 2, 'Two'))
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive == 3, "Three"))
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive == 4, 'Four'))
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive == 5, "Five"))
training$reg_nodes_pos <- with(training, replace(reg_nodes_pos, regional_nodes_positive > 5, 'More than five'))
table(training$reg_nodes_pos, exclude=NULL)

rm("regional_nodes_examined","regional_nodes_positive")

regional_nodes_examined <- factor(training[,'reg_nodes_ex'])
regional_nodes_positive <- factor(training[,'reg_nodes_pos'])

drops2 <- c("patient_id","cod", "primary_international", "surgery_of_primary_site_labeled","survival_days","type_follow_up_expected", "state_county",
            "seer_combined_summary_stage", "sequence_number")
breast4 <- training[ , !(names(training) %in% drops)]

fit_breast_cox_extension2 <- coxph(Surv(survival_days, status) ~ ethnicity +
                                     grade +
                                     seer_summary_stage +
                                     regional_nodes_examined +
                                     regional_nodes_positive +
                                     estrogen_receptor + 
                                     progesteron_receptor +
                                     ns(age, 3) +
                                     ns(tumor_size,3), data = training)

######################################################################
# Sul dataset di training ancora sottocampionamento sempre più grande
# per vedere se la dimensione campionaria impatta sul risultato
training <- sample_frac(breast, size = size, replace = FALSE)

hist(training$regional_nodes_examined, exclude=NA)
hist(training$regional_nodes_positive, exclude=NA)
table(regional_nodes_positive)

table(grade, seer_summary_stage)