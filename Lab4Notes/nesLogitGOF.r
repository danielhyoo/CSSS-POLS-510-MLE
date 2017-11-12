
## Example code for binary logit MLE:  Model Fitting
## Voting behavior example
##
## Christopher Adolph    faculty.washington.edu/cadolph
## 23 October 2016
##
##  Estimation by ML using optim() or by glm() on reduced dataset:
##     Model 1:  Age, Age^2, HS, College
##     Model 2:  Age, Age^2, HS, College, Married
##
##  Goodness of fit tests shown here:
##
##  Likelihood ratio test
##  Akaike Information Criterion
##  Bayesian Information Criterion
##  Deviance
##  Percent Correctly Predicted
##  Separation plots
##  Actual vs Predicted Plots
##  Error vs Predicted Plots
##  ROC plots
##  Residual vs Leverage Plots
##  Cross-validation

# Clear memory
rm(list=ls())

# Load libraries
library(simcf)
library(MASS)
library(nlme)
library(boot)            # For cv.glm()
library(separationplot)  # For separation plot
library(pscl)            # Alternative PCP code
library(verification)    # For ROC area
library(tile)            # For some graphics; used by plot.binPredict()
library(RColorBrewer)    # For nice colors
source("binaryGOF.R")    # Percent correctly predicted and concordance indexes
source("binPredict.R")   # Code for making predicted vs actual plots

# Get nice colors
col <- brewer.pal(5, "Set1")
blue <- col[2]
orange <- col[5]

# Models in R formula format
m1 <- vote00 ~ age + I(age^2) + hsdeg + coldeg
m2 <- vote00 ~ age + I(age^2) + hsdeg + coldeg + marriedo

# Note:  the variable marriedo is current marrieds,
#        the variable married is ever-marrieds

# Load data
setwd("~/desktop")
file <- "nes00a.csv"    
fulldata <- read.csv(file,header=TRUE)

# Keep only cases observed for all models
data <- extractdata(m2, fulldata, na.rm = TRUE)
attach(data)

# Construct variables and model objects
y <- vote00
x1 <- cbind(age,age^2,hsdeg,coldeg)
x2 <- cbind(age,age^2,hsdeg,coldeg,marriedo)

# Likelihood function for logit
llk.logit <- function(param,y,x) {
  os <- rep(1,length(x[,1]))
  x <- cbind(os,x)
  b <- param[ 1 : ncol(x) ]
  xb <- x%*%b
  sum( y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb)))
               # optim is a minimizer, so min -ln L(param|y)
}

# Fit logit model using optim
ls.result <- lm(y~x1)  # use ls estimates as starting values
stval <- ls.result$coefficients  # initial guesses
logit.m1 <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x1)
                   # call minimizer procedure
pe.m1 <- logit.m1$par   # point estimates
vc.m1 <- solve(logit.m1$hessian)  # var-cov matrix
se.m1 <- sqrt(diag(vc.m1))    # standard errors
ll.m1 <- -logit.m1$value  # likelihood at maximum

# Alternative estimation technique:  GLM
glm.m1 <- glm(m1, data=data, family="binomial")
print(summary.glm(glm.m1))


# Fit logit model with added covariate:  married
ls.result <- lm(y~x2)  # use ls estimates as starting values
stval <- ls.result$coefficients  # initial guesses
logit.m2 <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x2)
                   # call minimizer procedure
pe.m2 <- logit.m2$par   # point estimates
vc.m2 <- solve(logit.m2$hessian)  # var-cov matrix
se.m2 <- sqrt(diag(vc.m2))    # standard errors
ll.m2 <- -logit.m2$value  # likelihood at maximum

# GLM estimation of model with married
glm.m2 <- glm(m2, data=data, family="binomial")
print(summary.glm(glm.m2))


## Goodness of fit of model 1 and model 2

# Check number of parameters in each model
k.m1 <- length(pe.m1)
k.m2 <- length(pe.m2)

k.m1
k.m2

# Likelihood ratio (LR) test
lr.test <- 2*(ll.m2 - ll.m1)
lr.test.p <- pchisq(lr.test,df=(k.m2 - k.m1),lower.tail=FALSE)

lr.test.p

# Bayesian Information Criterion (BIC)
bic.m1 <- log(nrow(x1))*k.m1 - 2*ll.m1
bic.m2 <- log(nrow(x2))*k.m2 - 2*ll.m2
bic.test <- bic.m2 - bic.m1
bic.test

# Akaike Information Criterion (AIC)
aic.m1 <- 2*k.m1 - 2*ll.m1
aic.m2 <- 2*k.m2 - 2*ll.m2
aic.test <- aic.m2 - aic.m1
aic.test

# Deviance (the "-0" terms refer to the log-likelihood of the saturated model,
# which is zero for categorical outcomes)
deviance.m1 <- -2*(ll.m1 - 0)
deviance.m2 <- -2*(ll.m2 - 0)

# Percent correctly predicted (using glm result and my source code)

pcp.glm

pcp.null <- pcp.glm(glm.m1, vote00, type="null")
pcp.m1 <- pcp.glm(glm.m1, vote00, type="model")
pcp.m2 <- pcp.glm(glm.m2, vote00, type="model")
pcpi.m1 <- pcp.glm(glm.m1, vote00, type="improve")
pcpi.m2 <- pcp.glm(glm.m2, vote00, type="improve")

pcp.null
pcp.m1
pcp.m2
pcpi.m1
pcpi.m2

## Another way to cumpute PCP with the pscl package
#library(pscl)
#hitmiss(glm.m1)
#hitmiss(glm.m1, k=.3) #change the threshold

## Still another way with the DAMisc package 
#pre(glm.m1)

# Separation plots
separationplot(pred=glm.m1$fitted.values, actual=glm.m1$y)

separationplot(pred=glm.m2$fitted.values, actual=glm.m2$y)

# binPredict for Actual vs Predicted plots, Error vs Predicted plots, and ROC plots
# From binPredict.R source code

# We use a helper function binPredict() to compute bins and ROC curves for us.  
# The we can plot one or more models using the plot function

# Other options for binPredict():
#   bins = scalar, number of bins (default is 20)
#   quantiles = logical, force bins to same # of observations (default is FALSE)
#   sims = scalar, if sim=0 use point estimates to compute predictions;
#                  if sims>0 use (this many) simulations from predictive distribution
#                            to compute predictions (accounts for model uncertainty)
#          default is 100 simulations; note: ROC curves always use point estimates only

binnedM1 <- binPredict(glm.m1, col=blue, label="M1: Age, Edu", quantiles=TRUE)

binnedM2 <- binPredict(glm.m2, col=orange, label="M2: Age, Edu, Married", quantiles=TRUE)

## To make bins of equal probability width instead of equal # obs: 
#binnedM1b <- binPredict(glm.m1, col=blue, label="M1: Age, Edu", quantiles=FALSE)
#binnedM2b <- binPredict(glm.m2, col=orange, label="M2: Age, Edu, Married", quantiles=FALSE)


## Some options for plot.binPredict (more in source code)
##   together = logical, plot models overlapping on same plot (default is TRUE)
##   display = character, avp:  plot predicted actual vs predicted probs
##                        evr:  plot actual/predicted vs predicted probs
##                        roc:  plot receiver operator characteristic curves
##             default is c("avp", "evp", "roc") for all three
##   thresholds = numeric, show these thresholds on ROC plot (default is NULL)
##   hide = logical, do not show number of observations in each bin (default is TRUE)
##   ignore = scalar, do not show bins with fewer observations than this (default = 5)
##   totalarea = scalar, total area of all circles for a model relative to plot (default=0.1)
##   cex = scalar, size of numeric labels
##   showbins = logical, show bin boundaries
##   file = character, save result to a pdf with this file name

# Show actual vs predicted of M1 on screen
plot(binnedM1, display="avp", hide=TRUE, labx=0.35)

# Show actual vs predicted of M1 and M2 to file
plot(binnedM1, binnedM2, display="avp", hide=TRUE, labx=0.35)

# Send error vs predicted of M1 and M2 to file
plot(binnedM1, binnedM2, display="evp", hide=TRUE, labx=0.35)

# Send ROC plots for M1 and M2 to file
plot(binnedM1, binnedM2, display="roc", thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.35)

# Send actual vs predicted, error rate vs predicted, and ROC to file
plot(binnedM1, binnedM2, thresholds=c(0.9, 0.8, 0.7, 0.6),
     hide=TRUE, labx=0.35)


# Also see ROCR package for ROC curves and many other prediction metrics
# and the verification package for a rudimentary roc plot function roc.plot()

# Concordance Indexes / AUC (using glm result and my source code)
concord.null <- concord.glm(glm.m1, vote00, type="null")
concord.m1 <- concord.glm(glm.m1, vote00, type="model")
concord.m2 <- concord.glm(glm.m2, vote00, type="model")
concordi.m1 <- concord.glm(glm.m1, vote00, type="improve")
concordi.m2 <- concord.glm(glm.m2, vote00, type="improve")

concord.null
concord.m1
concord.m2
concordi.m1
concordi.m2

###  Residuals using glm version
hatscore.m1 <- hatvalues(glm.m1)/mean(hatvalues(glm.m1))
rstu.m1 <- rstudent(glm.m1)

hatscore.m2 <- hatvalues(glm.m2)/mean(hatvalues(glm.m2))
rstu.m2 <- rstudent(glm.m2)

usr <- c(0,10,-3,3)

plot.new()
par(usr=usr, tcl=-0.1, mgp=c(2,0.35,0))
axis(2, las=1)
par(usr=usr, tcl=-0.1, mgp=c(2,0.15,0))
axis(1, at=c(0,1,2,3,4,5,6,7,8,9,10))

title(xlab="Standardized hat-values",
      ylab="Studentized residuals")
points(hatscore.m1, rstu.m1,col = blue)
lines(c(usr[1], usr[2]), c(-2,-2), lty="dashed")
lines(c(usr[1], usr[2]), c(2,2), lty="dashed")
lines(c(3,3), c(usr[3], usr[4]), lty="dashed")

plot.new()
par(usr=usr,tcl=-0.1,mgp=c(2,0.35,0))
axis(2,las=1)
par(usr=usr,tcl=-0.1,mgp=c(2,0.15,0))
axis(1,at=c(0,1,2,3,4,5,6,7,8,9,10))

title(xlab="Standardized hat-values",
      ylab="Studentized residuals")
points(hatscore.m2, rstu.m2, col=orange)
lines(c(usr[1], usr[2]), c(-2,-2), lty="dashed")
lines(c(usr[1], usr[2]), c(2,2), lty="dashed")
lines(c(3,3), c(usr[3], usr[4]), lty="dashed")

### Cross-validation (takes a few minutes to run)

## A precent-correctly-predicted-style cost function
## r is actual y, pi is expected y
## Rate of inaccuracy: mean(vote00!=round(yp))
costpcp <- function(r, pi=0) mean(r!=round(pi))

## an alternative cost function for binary data
## cost <- function(r, pi=0) mean(abs(r-pi)>0.5)

cv.m1 <- cv.glm(data, glm.m1, costpcp)
cvPCP.m1 <- 1 - cv.m1$delta[2]

cv.m2 <- cv.glm(data, glm.m2, costpcp)
cvPCP.m2 <- 1 - cv.m2$delta[2]

cvPCP.m1
cvPCP.m2

#### More cross-validation

## A simple leave-one-out cross-validation function for logit glm; returns predicted probs
loocv <- function (obj) {
    data <- obj$data
    m <- dim(data)[1]
    form <- formula(obj)
    fam <- obj$family$family
    loo <- rep(NA, m)
    for (i in 1:m) {
        i.glm <- glm(form, data = data[-i, ], family = fam)
        loo[i] <- predict(i.glm, newdata = data[i,], family = fam, type = "response")
    }
    loo
}

# LOOCV for models 1 and 2
predCVm1 <- loocv(glm.m1)
predCVm2 <- loocv(glm.m2)

predCVm1
predCVm2

# Make cross-validated AVP and ROC plots; note use of newpred input in binPredict
binnedM1cv <- binPredict(glm.m1, newpred=predCVm1, col=blue, label="M1: LOO-CV", quantiles=TRUE)

plot(binnedM1cv, display=c("avp","roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25)

binnedM2cv <- binPredict(glm.m2, newpred=predCVm2, col=orange, label="M2: LOO-CV", quantiles=TRUE)

plot(binnedM2cv, display=c("avp","roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25)

plot(binnedM1cv, binnedM2cv, display=c("avp","roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25)






