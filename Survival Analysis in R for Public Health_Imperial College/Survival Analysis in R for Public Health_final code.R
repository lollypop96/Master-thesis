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
library(ggfortify) 

# Load dataset
g <- read.csv("C:/Users/liry9/Desktop/Tesi/Survival Analysis in R for Public Health_Imperial College/Week 1/simulated_HF.csv", header=TRUE, sep=',')
#g <- na.omit(g)
dim(g)
head(g)

# Define variables 
gender <- factor(g[,"gender"]) 
fu_time <- g[,"fu_time"]  
death <-  g[,"death"] 
age <- g[,"age"] 
copd <- factor(g[,"copd"]) 
ethnicgroup <- factor(g[,"ethnicgroup"]) 
quintile <- factor(g[,"quintile"]) 
ihd <- factor(g[,'ihd']) 
valvular <- factor(g[,'valvular_disease']) 
pvd <- factor(g[,'pvd']) 
stroke <- factor(g[,'stroke']) 
pneumonia <- factor(g[,'pneumonia']) 
renal <- factor(g[,'renal_disease']) 
ca <- factor(g[,'cancer']) 
mets <- factor(g[,'metastatic_cancer']) 
mental_health <- factor(g[,'mental_health']) 
ht <- factor(g[,"hypertension"]) 
cog_imp <- factor(g[,"senile"]) 
prior_dnas <- g[,"prior_dnas"]

# Plotting a Kaplan-Meier curve 
###################### 

# 1. Generate the survival curve 
km_fit <- survfit(Surv(fu_time, death) ~ 1)

# Output the probability of survival at certain times after hospital admission 
summary(km_fit, times = c(1:7,30,60,90*(1:10)))

# 1b. Alternative plot with ggplot2 
autoplot(km_fit) + theme_bw() # theme_bw() is a predesigned "theme" which makes the plot prettier 

# Plotting a Kaplan-Meier curve by gender 
###################### 

# 2. Generate the survival curve 
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

# 2. Plot the curve 
plot(km_gender_fit)

# 2b. Alternative plot with ggplot2 
autoplot(km_gender_fit) + theme_bw()


######################################################################
# Perform log rank test to see whether survival varies by gender 
survdiff(Surv(fu_time, death) ~ gender, rho = 0)

# Testing whether those over the age of 65 have different survival to those under it 
###################### 

# 1. Dichotomise age into categorical (binary in this case) variable 
age_65plus <- ifelse(g[,'age']>=65, 1, 0) 

# 2. Perform log rank test 
survdiff(Surv(fu_time, death) ~ age_65plus, rho = 0)

# Plot survival curve by age above or below 65 
###################### 

# 1. Generate survival curve 
km_old_fit <- survfit(Surv(fu_time, death) ~ age_65plus) 

# 2. Plot 
plot(km_old_fit)

# 2b. Alternative plot in ggplot2 
autoplot(km_old_fit) + theme_bw()

######################################################################
# Run Cox regression model with age as predictor (continuous variable) 
###################### 

# 1. Generate model 
cox <- coxph(Surv(fu_time, death) ~ age, data = g) 

# 2. Summarise model 
summary(cox)

######################
# Run Cox regression model with quintile as predictor (categorical variable) 
# Changing the reference group to first quintile 
# Removing the zero quintile altogether

# 1. Summarise the variable 
table(quintile, exclude = NULL)

# 2. Check levels 
levels(quintile)

# 3. Generate model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning

# 4. Summarise model
summary(cox)

# 5. Make the first quintile the reference group 
quintile <- relevel(quintile, ref = "1") 

# 6. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning
summary(cox) # still an issue where quintile = 0

# 7. Inspecting quintile variable 
table(quintile, g$death) # Only 4 entries for quintile = 0 and 100% didn't die 

# 8. Removing quintile = 0 entries as there are only 4 of them 
quintile_5groups <- quintile 
quintile_5groups[quintile_5groups == 0] <- NA # set the zeroes to missing 
quintile_5groups <- factor(quintile_5groups) # this removes 0 as a level as it is an empty category 

# 9. Regenerating the model and summarising 
cox <- coxph(Surv(fu_time, death) ~ quintile_5groups) 
summary(cox) # still an issue where quintile = 0

######################################################################
# Run Cox regression model with ethnic group as predictor (categorical variable) 
# Including missing values as another category 
###################### 

# 1. Summarise variable 
table(ethnicgroup, exclude = NULL)

# 2. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox)

# 3. Add another category (8) 
levels(ethnicgroup) <- c(levels(ethnicgroup),"8")  

# 4. Redefine NA as another group, 8 
ethnicgroup[is.na(ethnicgroup)] <- "8" 

# 5. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox)

######################################################################
###################### 



# Investigating our variables in order to best perform a Cox model with multiple predictors 
# Checking for missing values 
# Running a multiple Cox regression 
###################### 

# 1. Summarising age 
summary(g$age) # no NAs
hist(g$age)

# 2. Gender 
gender_table <- table(gender, exclude = NULL) 
addmargins(gender_table) # no NAs

round(100 * prop.table(gender_table), digits = 1) # Percentages rounded to 1 decimal place

# 3. Chronic Obstructive Pulmonary Disease (COPD) 
copd_table <- table(copd, exclude = NULL)  
addmargins(copd_table) # no NAs

round(100 * prop.table(copd_table), digits = 1) # Percentages rounded to 1 decimal place 

# 4. Prior OPD appointments missed  
prior_dnas_table <- table(prior_dnas, exclude = NULL)  
addmargins(prior_dnas_table) # no NAs

round(100 * prop.table(prior_dnas_table), digits = 1) # Percentages rounded to 1 decimal place 
# 5. Ethnic group 
ethnicgroup_table <- table(ethnicgroup, exclude = NULL)  
addmargins(ethnicgroup_table) # 4.3% NA

round(100 * prop.table(ethnicgroup_table), digits = 1) # Percentages rounded to 1 decimal place 

# 6. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + ethnicgroup) 
summary(cox)

######################################################################
# Investigating whether the assumptions of the Cox model are being broken 
# Testing for proportional hazards assumption (with gender as predictor variable) 
###################### 

# 1. Generate model fit 
fit <- coxph(Surv(fu_time, death) ~ gender) 

# 2. Apply the test to the model 
temp <- cox.zph(fit)     

# 3. Display results 
print(temp) 

# 4. Plot the curves 
plot(temp)

# 4b. Alternative plot in ggplot 
ggcoxzph(temp)

######################################################################
# Generating other diagnostic plots for Cox Proportional Hazards model 
###################### 

# 1. Define model 
res.cox <- coxph(Surv(fu_time, death) ~ age) 

# Generate diagnostic plots 

# 2. Plotting the estimated changes in the regression coefficients on deleting each patient 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

# 3. Plotting deviance residuals 
ggcoxdiagnostics(res.cox, type = "deviance", 
                 linear.predictions = FALSE, ggtheme = theme_bw())

# 4. Plotting Martingale residuals 
fit <- coxph(Surv(fu_time, death) ~ age + log(age) + sqrt(age)) 
ggcoxfunctional(fit, data = g) # note we must specify original dataframe

######################################################################
# Testing proportionality assumption 
# Testing for a statistical relationship between gender and time 
###################### 

# 1. Generate model with time-transform function (tt) 
fit <- coxph(Surv(fu_time, death) ~ gender + tt(gender))  

# 2. Summarise 
summary(fit)

######################################################################
# Backwards elimination to choose predictors for Cox regression 
###################### 

# 1. Run the full model with all of your predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + 
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               ca + mets + mental_health + cog_imp + los + prior_dna) 
summary(cox)

# 2. Run the model with only significant predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + mets + cog_imp) 
summary(cox)

# 3. Test proportionality assumption on these predictors 
cox.zph(cox) 

######################################################################
# 1. Run the full model with all of your predictors
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd +                 valvular + pvd + stroke + copd + pneumonia + ht + renal +                 ca + mets + mental_health + cog_imp + los + prior_dna)
summary(cox)
