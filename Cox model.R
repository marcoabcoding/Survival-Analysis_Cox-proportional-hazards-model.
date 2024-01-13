# Load the survival package
require("survival")

# Define the survival times, censoring status, and group membership
tempos=c(28, 89, 175, 195, 309, 377, 393, 421, 447, 462, 709, 744, 770, 1106, 1206, # tumor grande
         34, 88, 137, 199, 280, 291, 299, 300, 309, 351, 358, 369, 369, 370, 375, 382, 392, 429, 451, 1119)
cens=c(1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,# tumor grande
       1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0)
grupo=c(rep(1,15),rep(0,20))

# Fit a Cox proportional hazards model
ajuste<-coxph(Surv(tempos,cens)~factor(grupo), method="breslow")
summary(ajuste)$coef

# Calculate the confidence interval for the model coefficients
beta=ajuste$coefficients
var=ajuste$var
alfa=0.05
zalf=qnorm(1-alfa/2)
ic= beta + c(-1,1) * zalf * sqrt(var)
ic

# Display various tests from the Cox model summary
summary(ajuste)$waldtest
summary(ajuste)$logtest
summary(ajuste)$sctest
exp(ic)

# Perform and display the logrank test using Cox model output
summary(ajuste)$sctest # Logrank test using coxph output
teste=survdiff(Surv(tempos,cens)~factor(grupo),rho=0)
teste # Logrank test

# Section 2: Analysis of post-surgical treatments for ovarian cancer
require("survival")
require("survminer")
require("ggplot2")
require("ggfortify")
require("gridExtra")

# Load the ovarian cancer dataset
ovario = read.csv2('C:\\Users\\mkm\\Downloads\\cancer de ovario.csv')
head(ovario)

# Fit a Cox model with multiple covariates
fit<-coxph(Surv(tempo,falha)~factor(tratamento)+idade+factor(residuo)+factor(status), data=ovario,x = T, method="breslow")
summary(fit)$coef

# Analyze Schoenfeld residuals for the model with interaction
# Global test and by factor combination of Pearson correlation for Schoenfeld standardized residuals
ftest<-cox.zph(fit,transform="identity",terms=FALSE) # g(t)=t
ftest
ggcoxzph(ftest)

# Global test of Pearson correlation for Schoenfeld standardized residuals
ftest<-cox.zph(fit,transform="identity",terms=TRUE) # g(t)=t
ftest
ggcoxzph(ftest)

# Analyze Schoenfeld residuals for the model without interaction
# Global test and by factors of Pearson correlation for Schoenfeld standardized residuals
ftestse<-cox.zph(fit3,transform="identity",terms=FALSE) # g(t)=t
ftestse
ggcoxzph(ftestse)

# Global test of Pearson correlation for Schoenfeld standardized residuals
ftestse<-cox.zph(fit3,transform="identity",terms=TRUE) # g(t)=t
ftestse
ggcoxzph(ftestse)

# Display summary of the Cox model fit
summary(fit4)

# Section 3: Data Insertion and Model Fitting
# Load necessary packages
require("survival")
require("survminer")

# Insert data for survival times, censoring status, and group membership
tempos = c( 1, 4, 5, 6, 7, 7, 3, 5, 5, 5, 6, 1, 3, 4, 7, 7, 7, 3, 5, 7, 7, 7, 3, 5, 5, 7, 7, 7)
# Define the censoring status and group membership for the study
censuras = c( rep(1,4), 0, 0, 1, 1, 0, 0, rep(1,3), 0, 1, 1, 0, 1, 0, 1, 0, 0, rep(1,3), rep(0,3) )
grupos = c(rep("A",6), rep("B",5), rep("C",6), rep("D",5), rep("Controle",6))
grupos=as.factor(grupos)
levels(grupos)=c("Controle","A","B","C","D") # Ordering the factor levels

# Fit a Cox proportional hazards model for the groups
ajuste<-coxph(Surv(tempos,censuras)~grupos, method="breslow") # Adjusting the model
summary(ajuste) # Summary of model fitting
summary(ajuste)$coef

# Analyze Schoenfeld residuals without interaction
# Global test of Pearson correlation for Schoenfeld standardized residuals
ftest<-cox.zph(ajuste,transform="identity",terms=FALSE) # g(t)=t
ftest
ggcoxzph(ftest)

# Chi-square statistics for various factors in the model
# Format: factor name, chi-square statistic, degrees of freedom, p-value
# factor(tratamento) 0.323  1 0.570
# idade              0.415  1 0.520
# factor(residuo)    1.490  1 0.222
# factor(status)     2.817  1 0.093
# GLOBAL             4.505  4 0.342

#Based on the chi-square statistics, none of the factors - treatment, age, residual tumor, or patient status - 
#show a statistically significant impact on survival times. The global test also supports this finding, 
#indicating no significant collective effect of these factors.
