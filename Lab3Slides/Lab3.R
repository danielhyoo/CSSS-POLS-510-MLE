

#################################################################
# Lab Session 3: Maximum Likelihood Estimation in R
# Using simcf and tile to explore an estimated logistic regression
# Voting example using 2000 NES data after King, Tomz, and Wittenberg
# Chris Adolph  www.chrisadolph.com
#################################################################

#rm(list = ls()) # clear up the memory
# Penalize Scientific Notation if you want
options(scipen=20)
#setwd("~/Desktop") 

#install and load the packages needed

installed.packages()[,1] # see what packages are currently installed 

#from CRAN: install.packages("MASS", dependencies = TRUE)
library(MASS)
library(RColorBrewer)

#download simcf and tile packages from- http://faculty.washington.edu/cadolph/software
# don't unzip the archive (tar) file 
# from your local source: install.packages("~/YourDrive/YourFolder/tile_0.4.7_R_x86_64-apple-darwin9.8.0.tar", repos = NULL)
library(simcf)
library(tile)


# Load data
file <- "nes00a.csv"
data <- read.csv(file, header=TRUE)

# attach(data)  if you "attach" data, you don't have to specify data name every time you call a variable from the data (ex. data$vote00). However, if you are working on several different datasets at the same time, attaching data might cause confusion. You can use detach(data) when you are done with the dataset

# Estimate logit model using optim()
# Construct variables and model objects
y <- data$vote00
x <- cbind(data$age,data$hsdeg,data$coldeg)

# Likelihood function for logit
llk.logit <- function(param,y,x) {
  os <- rep(1,length(x[,1])) # constant 
  x <- cbind(os,x) # constant+covariates
  b <- param[ 1 : ncol(x) ] # number of parameters to be estimated equals number of columns in x (i.e, one for constant and one for each covariates : total 4)
  xb <- x%*%b  
  sum( y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb))) # log-likelihood function for logit model (based on our choice of standard logistic cdf as the systematic component)
               # optim is a minimizer, so use -lnL here
}

# Fit logit model using optim
ls.result <- lm(y~x)  # use ls estimates as starting values (for convenience)
stval <- ls.result$coefficients  # initial guesses
logit.result.opt <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x)
                   # call minimizer procedure or max by adding control=list(fnscale=-1)
pe.opt <- logit.result.opt$par   # point estimates
vc.opt <- solve(logit.result.opt$hessian)  # var-cov matrix
se.opt <- sqrt(diag(vc.opt))    # standard errors
ll.opt <- -logit.result.opt$value  # likelihood at maximum

logit.optim<-data.frame(cbind(round(pe.opt,3), round(se.opt,3)))
rownames(logit.optim)<-c("intercept", "age", "highschool" , "college")
colnames(logit.optim)<-c("pe", "std.err")
logit.optim

#p-value based on t-statistics
2*pt(abs(logit.optim$pe/logit.optim$std.err), df=length(y)-length(pe.opt) , lower.tail = FALSE)

# Estimate logit model using glm()

# Run logit & extract results using glm. GLM solves the likelihood equations with a common numeric algorithm called iteratively re-weighted least squares (IRWLS).

logit.result <- glm(vote00~age +hsdeg+coldeg, family=binomial, data=data) # family "binomial" calls logit transformation (log(pi/1-pi)) as a "link function" corresponding to logistic distribution. Link function transforms pi so that it follows a linear model. (So although pi itself is dependent on covariates in a non-linear way, logit transformed pi is dependent on covariates in a linear way.) 

summary(logit.result)

# now a new model adding age^2
model <- vote00 ~ age + I(age^2) + hsdeg + coldeg
mdata <- extractdata(model, data, na.rm=TRUE) # needs library(simcf)

logit.result <- glm(model, family=binomial, data=mdata) 

pe <- logit.result$coefficients  # point estimates
vc <- vcov(logit.result)         # var-cov matrix

# Simulate parameter distributions
sims <- 10000
simbetas <- mvrnorm(sims, pe, vc) #needs library(MASS) # draw 10000 sets of simulated parameter (beta) estimates from a multivariate normal distribution with mean pe and variance-covariance vc

# Now let's plan counterfactuals: We will have three sets of counterfactuals based on education lavel (less than hs edu, hs edu, college or higher edu), and for each set we will make age varies between 18 years old and 97 years old. 

# Set up counterfactuals:  all ages, each of three educations
xhyp <- seq(18,97,1)  # create age vector 
nscen <- length(xhyp)  # we will have total 80 different age scenarios for each education level

nohsScen <- hsScen <- collScen <- cfMake(model, mdata, nscen)  #this is just to initialize 80 scenarios for each education level. As default, all covariate values are set at the mean.

# Create three sets of education counterfactuals


for (i in 1:nscen) {
  # No High school scenarios (loop over each age, total 80 scenarios)
  nohsScen <- cfChange(nohsScen, "age", x = xhyp[i], scen = i)
  nohsScen <- cfChange(nohsScen, "hsdeg", x = 0, scen = i) # no hs degree
  nohsScen <- cfChange(nohsScen, "coldeg", x = 0, scen = i) # no college degree

  # HS grad scenarios (loop over each age, total 80 scenarios)
  hsScen <- cfChange(hsScen, "age", x = xhyp[i], scen = i)
  hsScen <- cfChange(hsScen, "hsdeg", x = 1, scen = i) # has hs degree
  hsScen <- cfChange(hsScen, "coldeg", x = 0, scen = i)  # no college degree

  # College grad scenarios (loop over each age, total 80 scenarios)
  collScen <- cfChange(collScen, "age", x = xhyp[i], scen = i)
  collScen <- cfChange(collScen, "hsdeg", x = 1, scen = i) # has hs degree
  collScen <- cfChange(collScen, "coldeg", x = 1, scen = i) # has college degree
}

# # Now given the counterfactual covariates (nohsScen/hsScen/collScen) and simulated parameters (simbetas), we can calculate expected value of the response. In this case, expected probability of voting!

head(nohsScen$x) #we will fit the counterfactual data
nohsScen$model # in the model specification
nohsSims <- logitsimev(nohsScen, simbetas, ci=0.95) # using simulated betas to get expected values. Built-in function "logitsimev" calculates the expected value for every individual scenario you created.


# don't run 
for (i in 1:nscen) {
  
  simmu <- simbetas %*% scenario[i, ] # matrix of simmulated beta coefficients (10000*5) x each set of counterfactual covariates including 1 for constant (5*1), so we get 10000 simmulated version of xb per scenario  
 
  simy <- 1/(1 + exp(-simmu)) # standard logistic distribution as the functional form of pi, we get 10000 simulated pi per scenario (which is also the expected value in case of bernoulli)
  
  
# For each set of simy (total 80 sets), logitsimev output provides its mean as the "pe", 2.5% and 97.5% quantiles as "lower" and "upper" confidence intervals at 0.95.    
  
  }

nohsSims #reports lower and upper confidence intervals as well as expected probabilities


# same thing for two other sets of sceneraios
hsSims <- logitsimev(hsScen, simbetas, ci=0.95)
collSims <- logitsimev(collScen, simbetas, ci=0.95)

  

# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")

# Set up lineplot traces of expected probabilities

#Traces are elements of tile : lines, labels, legends...

# no hs 
nohsTrace <- lineplot(x=xhyp, # age on x-axis
                      y=nohsSims$pe, #expected probability on y-axis
                      lower=nohsSims$lower, # lower confidence interval
                      upper=nohsSims$upper, #upper confidence interval
                      col=col[1], #color choice 
                      extrapolate=list(data=mdata[,2:ncol(mdata)], #actual covariates (i.e., values in your data)
                               cfact=nohsScen$x[,2:ncol(hsScen$x)],#counterfactual covariates
                        omit.extrapolated=FALSE), #don't show extrapolated values 
                      plot=1)

# hs but no college
hsTrace <- lineplot(x=xhyp,
                    y=hsSims$pe,
                    lower=hsSims$lower,
                    upper=hsSims$upper,
                    col=col[2],
                    extrapolate=list(data=mdata[,2:ncol(mdata)],
                                       cfact=hsScen$x[,2:ncol(hsScen$x)],
                                       omit.extrapolated=FALSE),
                    plot=1)

#college 
collTrace <- lineplot(x=xhyp,
                      y=collSims$pe,
                      lower=collSims$lower,
                      upper=collSims$upper,
                      col=col[3],
                      extrapolate=list(data=mdata[,2:ncol(mdata)],
                                       cfact=collScen$x[,2:ncol(hsScen$x)],
                                       omit.extrapolated=FALSE),
                      plot=1)

# Set up traces with labels

labelTrace <- textTile(labels=c("Less than HS", "High School", "College"),
                       x=c( 55,    49,     30),
                       y=c( 0.26,  0.56,   0.87),
                       col=col,
                       plot=1)

# For legend

legendTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                        x=c(82, 82, 82),
                        y=c(0.2, 0.16, 0.12),
                        plot=1)


#options(device="quartz")

# Plot traces using tile
voting<-tile(nohsTrace,
     hsTrace,
     collTrace,
     labelTrace,
     legendTrace,
     limits=c(18,94,0,1),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5),
     output=list(file="educationEV",width=5.5)
)



################################################################
#
# Now consider a new specification adding the variable
# "ever married", or marriedo
#
# We will estimate this new model with glm(), then
# simulate new scenarios for marrieds and non-marrieds


# Estimate logit model using glm()


# Set up a new model formula and model specific data frame

model2 <- vote00 ~ age + I(age^2) + hsdeg + coldeg + marriedo
mdata2 <- extractdata(model2, data, na.rm=TRUE)

# Run logit & extract results
logit.m2 <- glm(model2, family=binomial, data=mdata2)
pe.m2 <- logit.m2$coefficients  # point estimates
vc.m2 <- vcov(logit.m2)         # var-cov matrix


# Simulate parameter distributions
sims <- 10000
simbetas.m2 <- mvrnorm(sims, pe.m2, vc.m2)


# Set up counterfactuals:  all ages

xhyp <- seq(18,97,1)
nscen <- length(xhyp)
marriedScen <- notmarrScen <- cfMake(model2, mdata2, nscen)
for (i in 1:nscen) {
  

  #  - we will use the marriedScen counterfactuals in FDs and RRs as well as EVs
  # Note below the careful use of before scenarios (xpre) and after scenarios (x) :i.e., use of the same age range (18-97) for both x and xpre, only marriedo values differ.
  
  # Married (loop over each age)
  marriedScen <- cfChange(marriedScen, "age", x = xhyp[i], xpre= xhyp[i], scen = i)
  marriedScen <- cfChange(marriedScen, "marriedo", x = 1, xpre= 0, scen = i)
  
  # Not Married (loop over each age)
  notmarrScen <- cfChange(notmarrScen, "age", x = xhyp[i], scen = i)
  notmarrScen <- cfChange(notmarrScen, "marriedo", x = 0, scen = i)
}

# Simulate expected probabilities for all age scenarios for married and not married respectively
marriedSims <- logitsimev(marriedScen, simbetas.m2, ci=0.95) 
notmarrSims <- logitsimev(notmarrScen, simbetas.m2, ci=0.95) 

# Simulate first difference of voting wrt marriage: E(y|married)-E(y|notmarried)
marriedFD <- logitsimfd(marriedScen, simbetas.m2, ci=0.95)

# Simulate relative risk of voting wrt marriage: E(y|married)/E(y|notmarried)
marriedRR <- logitsimrr(marriedScen, simbetas.m2, ci=0.95) 


## Make plots using tile

# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")

# Set up lineplot traces of expected probabilities
marriedTrace <- lineplot(x=xhyp,
                         y=marriedSims$pe,
                         lower=marriedSims$lower,
                         upper=marriedSims$upper,
                         col=col[1],
                         extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                                          cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                                          omit.extrapolated=TRUE),
                         plot=1)

notmarrTrace <- lineplot(x=xhyp,
                         y=notmarrSims$pe,
                         lower=notmarrSims$lower,
                         upper=notmarrSims$upper,
                         col=col[2],
                         ci = list(mark="dashed"),
                         extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                                          cfact=notmarrScen$x[,2:ncol(notmarrScen$x)],
                                          omit.extrapolated=TRUE),
                         plot=1)


# Set up traces with labels and legend
labelTrace <- textTile(labels=c("Currently Married", "Not Married"),
                       x=c( 35,    53),
                       y=c( 0.8,  0.56),
                       col=col,
                       plot=1)

legendTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                        x=c(80, 80, 80),
                        y=c(0.2, 0.15, 0.10),
                        cex=0.9,
                        plot=1)

# Plot traces using tile
tile(marriedTrace,
     notmarrTrace,
     labelTrace,
     legendTrace,
     limits=c(18,94,0,1),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5),
     output=list(file="marriedEV",width=5.5)
)



# Plot First Difference

# Set up lineplot trace of first difference

marriedFDTrace <- lineplot(x=xhyp,
                           y=marriedFD$pe,
                           lower=marriedFD$lower,
                           upper=marriedFD$upper,
                           col=col[1],
                           extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                                            cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                                            omit.extrapolated=TRUE),
                           plot=1)


# Set up baseline: for first difference, this is 0
baseline <- linesTile(x=c(18,94),
                      y=c(0,0),
                      plot=1)

# Set up traces with labels and legend
labelFDTrace <- textTile(labels=c("Married compared \n to Not Married"),
                         x=c( 40),
                         y=c( 0.20),
                         col=col[1],
                         plot=1)

legendFDTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                          x=c(80, 80, 80),
                          y=c(-0.02, -0.05, -0.08),
                          cex=0.9,
                          plot=1)

# Plot traces using tile
tile(marriedFDTrace,
     labelFDTrace,
     legendFDTrace,
     baseline,
     limits=c(18,94,-0.1,0.5),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Difference in Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5),
     output=list(file="marriedFD",width=5.5)
)


# Plot Relative Risk

# Set up lineplot trace of relative risk
marriedRRTrace <- lineplot(x=xhyp,
                           y=marriedRR$pe,
                           lower=marriedRR$lower,
                           upper=marriedRR$upper,
                           col=col[1],
                           extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                                            cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                                            omit.extrapolated=TRUE),
                           plot=1)


# Set up baseline: for relative risk, this is 1

baseline <- linesTile(x=c(18,94),
                      y=c(1,1),
                      plot=1)

# Set up traces with labels and legend
labelRRTrace <- textTile(labels=c("Married compared \n to Not Married"),
                         x=c( 55),
                         y=c( 1.25),
                         col=col[1],
                         plot=1)

legendRRTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                          x=c(80, 80, 80),
                          y=c(0.98, 0.95, 0.92),
                          cex=0.9,
                          plot=1)

# Plot traces using tile
tile(marriedRRTrace,
     labelRRTrace,
     legendRRTrace,
     baseline,
     limits=c(18,94,0.9,1.5),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Relative Risk of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5),
     output=list(file="marriedRR",width=5.5)
)





















# Example code for binary logit MLE:  Model Fitting
# Voting behavior example

#
#  Estimation by ML using optim() on reduced dataset
#     Model 1:  Age, HS, College
#     Model 2:  Age, HS, College, Married
#     Model 3:  Age, Age^2, HS, College, Married
#
#  Likelihood ratio test
#  Akaike Information Criterion
#  Bayesian Information Criterion
#  Average vs Predicted Plots
#  ROC plots
#  Residual vs Leverage Plots


# Load libraries
library(MASS)
library(nlme)
library(verification)
source("avp.r")          # Average vs Predicted plotting code

# Load data
file <- "nes00a.csv";    # Reduced to only shared obs
data <- read.csv(file,header=TRUE);
attach(data)

# Construct variables and model objects
y <- vote00
x <- cbind(age,hsdeg,coldeg)

# Likelihood function for logit
llk.logit <- function(param,y,x) {
  os <- rep(1,length(x[,1]))
  x <- cbind(os,x)
  b <- param[ 1 : ncol(x) ]
  xb <- x%*%b
  sum( y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb)));
               # optim is a minimizer, so min -ln L(param|y)
}

# Fit logit model using optim
ls.result <- lm(y~x)  # use ls estimates as starting values
stval <- ls.result$coefficients  # initial guesses
logit.result <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x)
                   # call minimizer procedure
pe.1 <- logit.result$par   # point estimates
vc.1 <- solve(logit.result$hessian)  # var-cov matrix
se.1 <- sqrt(diag(vc.1))    # standard errors
ll.1 <- -logit.result$value  # likelihood at maximum


# Fit logit model with added covariate:  married
x <- cbind(age,hsdeg,coldeg,marriedo)
ls.result <- lm(y~x)  # use ls estimates as starting values
stval <- ls.result$coefficients  # initial guesses
logit.result <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x)
                   # call minimizer procedure
pe.2 <- logit.result$par   # point estimates
vc.2 <- solve(logit.result$hessian)  # var-cov matrix
se.2 <- sqrt(diag(vc.2))    # standard errors
ll.2 <- -logit.result$value  # likelihood at maximum

# GOF of added variable

# LR test
lr.test <- 2*(ll.2 - ll.1)
lr.test.p <- pchisq(lr.test,df=1,lower.tail=FALSE)

# BIC
bic.test <- - 2*(ll.2 - ll.1) + 1*log(nrow(x))

# AIC
aic.test <- 2*(ll.2 - ll.1) - 1*2

# Deviance
deviance <- -2*(ll.2)

# Act v Pred plot
avp(y,
    x=cbind(rep(1,nrow(x)),x),
    beta=pe.2,
    fnform="logit",
    cutpoints=seq(0,1,0.1),
    usr=c(0,1,0,1),
    sizefactor=.25,
    color = "blue",
    output = list(outfile="nes_logit_AvP_0.eps",high=6,wide=5.5,epsi=FALSE),
    lab = list(x = .7, y=0.175, str="Model with married", col="blue", cex=1),
    ylab = "Actual Voting Rate, by bins of predicted",
    xlab = "Predicted Voting Rate, binned",
    closeplot=T)

# Act v Pred plot
avp(y,
    x=cbind(rep(1,nrow(x)),x),
    beta=pe.2,
    fnform="logit",
    cutpoints=seq(0,1,0.1),
    usr=c(0,1,0,1),
    sizefactor=.25,
    color = "blue",
    output = list(outfile="nes_logit_AvP.eps",high=6,wide=5.5,epsi=FALSE),
    lab = list(x = .7, y=0.175, str="Model with married", col="blue", cex=1),
    ylab = "Actual Voting Rate, by bins of predicted",
    xlab = "Predicted Voting Rate, binned",
    closeplot=F)

x <- cbind(age,hsdeg,coldeg)
avp(y,
    x=cbind(rep(1,nrow(x)),x),
    beta=pe.1,
    fnform="logit",
    cutpoints=seq(0,1,0.1),
    usr=c(0,1,0,1),
    sizefactor=.25,
    color = "red",
    output = list(outfile="nes_logit_AvP.eps",high=6,wide=5.5,epsi=FALSE),
    lab = list(x = .7, y=0.1, str="Model w/o married", col="red", cex=1),
    addtoplot=T,
    closeplot=T)

# ROC
yhat.1 <- 1/(1+exp(-cbind(rep(1,nrow(x)),x)%*%pe.1))
postscript("roc_m1.eps",paper="letter",pointsize = 14,width=6,height=5.5,horizontal = FALSE, onefile = TRUE)
roc.plot(y,cbind(yhat.1))
dev.off()

yhat.2 <- 1/(1+exp(-cbind(rep(1,nrow(x)),x,marriedo)%*%pe.2))
postscript("roc_m1m2.eps",paper="letter",pointsize = 14,width=6,height=5.5,horizontal = FALSE, onefile = TRUE)
roc.plot(y,cbind(yhat.1,yhat.2))
dev.off()


###  Residuals
glm.result <- glm(y~age + hsdeg + coldeg + marriedo,family="binomial")
summary.glm(glm.result)
hatscore <- hatvalues(glm.result)/mean(hatvalues(glm.result))
rstu <- rstudent(glm.result)

#x11()
#plot(hatscore,rstu)

usr <- c(0,10,-3,3)

postscript("nes_resid1.eps",paper="letter",pointsize = 14,width=5.5,height=5,horizontal = FALSE, onefile = TRUE);
plot.new()
par(usr=usr,tcl=-0.1,mgp=c(2,0.35,0))
axis(2,labels = T, tick = T, line = 0, outer = F, font = NA, fg = "black",las=1);
par(usr=usr,tcl=-0.1,mgp=c(2,0.15,0))
axis(1,at=c(0,1,2,3,4,5,6,7,8,9,10),labels = T, tick = T, line = 0, outer = F, font = NA, fg = "black");

title(xlab="Standardized hat-values",ylab="Studentized residuals")
points(hatscore,rstu,col = "red")
lines(c(usr[1],usr[2]),c(-2,-2),lty="dashed")
lines(c(usr[1],usr[2]),c(2,2),lty="dashed")
lines(c(3,3),c(usr[3],usr[4]),lty="dashed")

dev.off()


# Fit logit model with added covariate:  married & age^2
x <- cbind(age,hsdeg,coldeg,marriedo,age^2)
ls.result <- lm(y~x)  # use ls estimates as starting values
stval <- ls.result$coefficients  # initial guesses
logit.result <- optim(stval,llk.logit,method="BFGS",hessian=T,y=y,x=x)
                   # call minimizer procedure
pe.3 <- logit.result$par   # point estimates
vc.3 <- solve(logit.result$hessian)  # var-cov matrix
se.3 <- sqrt(diag(vc.2))    # standard errors
ll.3 <- -logit.result$value  # likelihood at maximum

###  Residuals
age2 <- age*age
glm.result <- glm(y~age + hsdeg + coldeg + marriedo + age2,family="binomial")
summary.glm(glm.result)
hatscore2 <- hatvalues(glm.result)/mean(hatvalues(glm.result))
rstu2 <- rstudent(glm.result)
#x11()
#plot(hatscore2,rstu2)

usr <- c(0,10,-3,3)

postscript("nes_resid2.eps",paper="letter",pointsize = 14,width=5.5,height=5,horizontal = FALSE, onefile = TRUE);
plot.new()
par(usr=usr,tcl=-0.1,mgp=c(2,0.35,0))
axis(2,labels = T, tick = T, line = 0, outer = F, font = NA, fg = "black",las=1);
par(usr=usr,tcl=-0.1,mgp=c(2,0.15,0))
axis(1,at=c(0,1,2,3,4,5,6,7,8,9,10),labels = T, tick = T, line = 0, outer = F, font = NA, fg = "black");

title(xlab="Standardized hat-values",ylab="Studentized residuals")
points(hatscore,rstu,col = "red")
points(hatscore2,rstu2, col = "blue")
lines(c(usr[1],usr[2]),c(-2,-2),lty="dashed")
lines(c(usr[1],usr[2]),c(2,2),lty="dashed")
lines(c(3,3),c(usr[3],usr[4]),lty="dashed")

dev.off()
