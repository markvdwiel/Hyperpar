#######################################################################################

# Title: Estimation of variance components, heritability and the ridge penalty in high-dimensional generalized linear models
# Authors: Jurre Veerman (linear model, plotting), Gwenael Leday (Bayes), Mark van de Wiel (linear model, GLM; mark.vdwiel@vumc.nl)
# Code for 'fridge' kindly provided by Kristoffer Hellton
# Date: 1/2/2019

# README: This code reproduces the simulations and plots of the corresponding manuscript 
# Content
# 1: Linear model simulations
#    - Random effects (ridge regression)
#    - Non-Gaussian priors + errors
#    - Data-based covariates
#    - Mixed effects
# 2: Data analysis example: weight gain data
# 3: Bayesian estimators: re-analysis of simulation by Moran, Rockova, George (2018)
# 4: GLM: Poisson and Binomial regression
# 5: Instructions for batch and parallel computations


#######################################################################################




#######################################################################################
################################ LINEAR MODEL #########################################
#######################################################################################

# Empty work space
rm(list=ls())

# Libraries required for estimation
library(mvtnorm)
library(glmnet)
library(HiLMM)
library(corpcor)

# Libraries required for plotting
library(ggplot2)
library(gridExtra)
library(reshape2)

# Set working directory 
wd <- "C:/Synchr/Jurre/Hyperpar/Github/SourceAndDataFiles"
#wd <- "/home/gwenael/Documents/project_MvdW/Jurre/2019_02_04_paper/HyperparNoKIRC"
setwd(wd)

# Functions needed
source("HyperparEst_linearSource.R")
source("plots_HyperparEst.R")

# Number of simulated data sets
nr.iter = 100

############
# Simulation 1: independent X, n=100, p=1000, gaussian prior beta
pmt <- proc.time(); set.seed(6326) #for reproducibility
simulate1 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i)); 
  wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0 = 0, cor = 0, data = F)})
proc.time()-pmt
save(simulate1, file="simulate1.Rdata")

# Plot Fig 1 (a)
#load("simulate1.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate1) 

############
# Simulation 2: independent X, large setting, n=1000, p=15000, gaussian prior beta
#FOR SPEEDING UP COMPUTATIONS: CONSIDER CODE FOR PARALLEL COMPUTING BELOW
pmt <- proc.time(); set.seed(6326) 
simulate2 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i)); 
  wrapper_basic(tau2_true = 0.01, sigma2_true = 150, n = 1000, p = 15000, p0 = 0, cor = 0, data = F)})
proc.time()-pmt
save(simulate2, file="simulate2.Rdata")

# Plot Fig 1 (b)
#load("simulate2.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate2) 

############
# Simulation 3: multicollinear X (blocks of size 10, rho=0.5), gaussian prior beta
pmt <- proc.time(); set.seed(6326) 
simulate3 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0 = 0, cor = 0.5, data = F, pblock=10)})
proc.time()-pmt
save(simulate3, file="simulate3.Rdata")

# Plot Fig 2 (a)
#load("simulate3.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate3) 


############
# Simulation 4: Spike-and-slab prior for beta with variance as for gaussian
pmt <- proc.time(); set.seed(6326) 
simulate4 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0 = 0.9, cor = 0, data = F)})
proc.time()-pmt
save(simulate4, file="simulate4.Rdata")

# Plot Supp Fig 3 (a)
#load("simulate4.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate4) 

############
# Simulation 5: Laplace prior for beta with variance as for gaussian
pmt <- proc.time(); set.seed(6326) 
simulate5 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0=0, cor = 0, data = F, simmodel= "laplace")})
proc.time()-pmt
save(simulate5, file="simulate5.Rdata")

# Plot Fig 1 (c)
#load("simulate5.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate5) 

############
# Simulation 6: Uniform prior for beta with variance as for gaussian
pmt <- proc.time(); set.seed(6326)
simulate6 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0=0, cor = 0, data = F, simmodel= "uniform")})
proc.time()-pmt
save(simulate6, file="simulate6.Rdata")

# Plot Supp Fig 3 (b)
#load("simulate6.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate6) 


############
# Simulation 7: Gaussian for beta, but t4 errors, scaled to have variance equal to sigma2_true
pmt <- proc.time(); set.seed(6326)
simulate7 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
 wrapper_basic(tau2_true = 0.01, sigma2_true = 10, n = 100, p = 1000, p0=0, cor = 0, data = F, simmodel= "t4error")})
proc.time()-pmt
save(simulate7, file="simulate7.Rdata")

# Plot Supp Fig 3 (c)
#load("simulate7.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate7) 


############
# Simulation 8: Data covariates large, response simulated. KIRC data (n=71, p=18391)
#FOR SPEEDING UP COMPUTATIONS: CONSIDER CODE FOR PARALLEL COMPUTING BELOW
load("KIRCdat.Rdata")
X <- t(KIRCdat)
tau2_true = 0.01; p<- ncol(X)
sigma2_true = p*tau2_true  #such that heritability =0.5
sigma2_true
pmt <- proc.time(); set.seed(6326) 
simulate8 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = sigma2_true, data = X)})
proc.time()-pmt
save(simulate8, file="simulate8.Rdata")

# Plot Fig 2(b)
#load("simulate8.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate8) 


############
# Simulation 9: Data covariates small, response simulated. TCGA Protein data (n=408, p=224)
tcpaOV <- read.table(file="TCGA-OV-L3-S35.csv", sep=",", header=TRUE) ## Read data
datamatOV <- data.matrix(tcpaOV[,-(1:3)]) # Remove annotations
rownames(datamatOV) <- tcpaOV$Sample_ID #set row names
X <- datamatOV
tau2_true = 0.01; p<- ncol(X)
sigma2_true = p*tau2_true  #such that heritability =0.5
sigma2_true
pmt <- proc.time(); set.seed(6326) 
simulate9 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = sigma2_true, data = X)})
proc.time()-pmt
save(simulate9, file="simulate9.Rdata")

# Plot Fig 2 (c)
load("simulate9.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate9) 


############
# Simulation 10: Basic mixed model, p=1000,n=100,m=10. Betas sampled from spike and slab: p0=0.9; Fixed effects alpha also 
#from spike-and-slab with p0=0.5 and larger variance
pmt <- proc.time(); set.seed(6326) 
simulate10 = lapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_mixed(tau2_true = 0.01, sigma2_true = 10, tau2_truefixed=0.10, p0=0.9, nr.fixed=10)})
proc.time()-pmt
save(simulate10, file="simulate10.Rdata")

# Plot Fig 4
#load("simulate10.Rdata")
set.seed(6326)
pl <- plotHyperpar(simulate10, mixed=TRUE)


#######################################################################################
###################### DATA ANALYSIS: WEIGHT GAIN DATA ################################
#######################################################################################

# Empty work space
rm(list=ls())

# Load data
load('E-GEOD-33070_gene.RData')
load('E-GEOD-33070_clinical.RData')

# Functions needed (Code for fridge by Kristoffer Hellton)
source("HyperparEst_linearSource.R")

# Include gene selection in outer LOOCV loop? ngene: number of genes to select a priori (1000)
# CVgenesel <- FALSE reproduces results of main document
# CVgenesel <- TRUE  reproduces results in Supplementary Material
CVgenesel <- TRUE; ngene <- 1000

# Script "weight_gain.R" performs prior gene selection and runs ridge + lambdaCV, ridge + lambdaMML, fridge + lambdaCV, fridge + lambdaMML
# Performs loocv to assess predictive performance
# Returns predictions and quadratic errors
source("weight_gain.R")

# All predictions, and true weight gain value
cbind(true=y0.vec, y0cv, y0.fridge, y0mml, y0mml.fridge)

# All quadratic errors, and mean errors
errors <- cbind(errcv, errcv.fridge, errmml, errmml.fridge)
colMeans(errors)

# Relative differences in error to ridge+lambdagcv: fridge, ridge + lambdamml, fridge + lambdamml
(mean(errcv)-mean(errcv.fridge))/mean(errcv)
(mean(errcv)-mean(errmml))/mean(errcv)
(mean(errcv)-mean(errmml.fridge))/mean(errcv)

# Absolute errors: ridge + lambdaGCV, fridge, ridge + lambdamml
errcv2 <- abs(y0cv-y.raw); errf2 <- abs(y0.fridge-y.raw); errm2 <- abs(y0mml-y.raw)

# Sort according to absolute errors ridge lambdaGCV
ordcv <- order(errcv2)

# Plot absolute prediction errors (Fig 3 and Supp Fig 2) 
n <- length(y.raw)
par(mar=c(2.5,2.5,1,1)) #b,l,t,r
plot(errcv2[ordcv],pch=15, xlab="")
points(1:n+0.05,errf2[ordcv]+0.05,col="red",pch=16)
points(1:n+2*0.05,errm2[ordcv],col="blue",pch=17)
legend(x=1,y=10,legend=c("GCV","fridge", "MML"),pch=c(15,16,17),col=c("black","red","blue")) 


#######################################################################################
## BAYESIAN ESTIMATORS, Re-analysis of simulation by Moran, Rockova, George (2018) ####
#######################################################################################

# Sparse setting from Moran, Rockova, George (2018)
# Compute estimates of sigma2, using conjugate Bayes model with co-estimation of tau2 [by MML], using maximum likelihood, and using a fixed, but misspecified tau2 = 100 (as in Moran et al)

# Empty work space
rm(list=ls())

# Simulation settings as in the Moran et al paper
tausqfix <- 100;CorX=0;n<-100;vareps<-3; sigmatrue<-sqrt(vareps); p0 <- 84;p1 <- 6; p=p0+p1; betanonzero <- c(-2.5,-2,-1.5,1.5,2,2.5)
n.iter <- 100

# Libraries
library(corpcor)
library(ggplot2)
library(reshape2)

# Function needed
source("HyperparEst_linearSource.R")  

# Simulate covariates and response, n.iter times
sigmasB <- sigmasML <-  sigmastau2fix <- c()
set.seed(6326)
for(k in 1:n.iter){
  # k<-1
  X <- matrix(rnorm(n*p),n,p)
  X <- t((t(X) - apply(t(X),1,mean))/apply(t(X),1,sd))
  Beta <- c(rep(0,p0), betanonzero)
  Y <- as.numeric(X%*%Beta ) + rnorm(n,0,sd=sqrt(vareps))
  Y <- matrix(Y,nrow=1)
  Yn <- as.numeric(Y)
  BR <- BayesRidge(Yn,X,plotML=FALSE)
  sigmaB <- sqrt(BR$postsig[2]/(BR$postsig[1]-1))
  sigmasB <- c(sigmasB,sigmaB)
  
  # ML estimator sigma^2
  n <- nrow(X)
  p <- ncol(X)
  ls <- lm(Yn ~ 0 + X)
  sigmaML <- sqrt(sum(ls$residuals^2)/(n-p-1))
  sigmasML <- c(sigmasML,sigmaML)
  
  # sigma^2 estimator for fixed prior variance of beta  
  inv <- solve(t(X) %*% X + tausqfix^(-1)*diag(p))
  H <- X %*% inv %*% t(X)
  Yr <- matrix(Y,nrow=n)
  newsigma <- sqrt((t(Yr) %*% (diag(n)-H) %*% Yr)/(n-1))
  sigmastau2fix <- c(sigmastau2fix, newsigma)
}
allres <- data.frame(ML=sigmasML,BayesFix=sigmastau2fix,BayesEB=sigmasB)
save(allres, file="resMoran_sim.Rdata")

# Plot Supp Fig 1 
allres2 <- suppressMessages(melt(allres))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
color.index_sigma = ifelse(is.wholenumber(as.numeric(allres2$variable)/2), color <- "red", color<-"blue")

set.seed(6326)
ggplot(data = allres2, aes(x = variable,  y = value)) +
  geom_boxplot(aes(fill = color.index_sigma), outlier.shape=NA) +
  geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
  geom_hline(yintercept = sigmatrue, color = "red", size = 1, linetype = "dashed") + 
  labs(y = expression(hat(sigma)), x = " ") + theme_classic() + 
  scale_fill_brewer(palette = "BuGn")+
  theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1,size=10),
        plot.margin=unit(c(10,2,-3,2),"mm"), #t,r,b,l
        panel.background = element_rect(colour = "black", size=0.5))


#######################################################################################
####################### GLM: POISSON and BINOMIAL REGRESSION ##########################
#######################################################################################

# Empty work space
rm(list=ls())

# Libraries
library(mgcv)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(glmnet)
library(penalized)

# Function needed
source("HyperparEst_linearSource.R") #source for simulating linear predictor
source("MMLGLMridge_poisson.R") #source for Poisson ridge
source('MMLGLMridge_binomial.R') #source for Binomial ridge
source("plots_HyperparEst.R") #source for gg-plots

############
# Poisson Simulation 1: independent X, n=100, p=1000, gaussian prior beta
nr.iter <- 100
set.seed(6326) 
pmt <- proc.time()
simulatep1 = lapply(1:nr.iter, function(i){wrapper_poisson(lambda_true = 100, n = 100, p = 1000, p0 = 0, cor = 0)})
proc.time()-pmt
save(simulatep1, file="simulatep1.Rdata")

# Plot Fig 6 (a)
#load("simulatep1.Rdata")
set.seed(6326)
pl <- plotHyperparlam(simulatep1, plotlog=TRUE)

############
# Poisson Simulation 2: multicollinear X (blocks of size 10, rho=0.5), gaussian prior beta
nr.iter <- 100
set.seed(6326) 
pmt <- proc.time()
simulatep2 = lapply(1:nr.iter, function(i){wrapper_poisson(lambda_true = 100, n = 100, p = 1000, p0 = 0, cor = 0.5)})
proc.time()-pmt
save(simulatep2, file="simulatep2.Rdata")

# Plot Fig 6 (b)
#load("simulatep2.Rdata")
set.seed(6326)
pl <- plotHyperparlam(simulatep2, plotlog=TRUE)

############
# Binomial(N,p) Simulation 1: independent X, n=100, p=1000, gaussian prior beta. N=5.
nr.iter <- 100
set.seed(6326) 
pmt <- proc.time()
simulateb1N5 = lapply(1:nr.iter, function(i){wrapper_binomial(lambda_true = 100, n = 100, p = 1000, p0 = 0, cor = 0, N=5)})
proc.time()-pmt
save(simulateb1N5, file="simulateb1N5.Rdata")

# Plot Supp Fig 4 (a)     # COULD REPRODUCE
#load("simulateb1N5.Rdata")
set.seed(6326)
pl <- plotHyperparlam(simulateb1N5, plotlog=TRUE, methinclude = c("MML_mgcv", "GCV_mgcv", "CV_glmnet"))

############
# Binomial(N,p) Simulation 2: multicollinear X (blocks of size 10, rho=0.5), gaussian prior beta. N=5.
nr.iter <- 100
set.seed(6326) 
pmt <- proc.time()
simulateb2N5 = lapply(1:nr.iter, function(i){wrapper_binomial(lambda_true = 100, n = 100, p = 1000, p0 = 0, cor = 0.5, N=5)})
proc.time()-pmt
save(simulateb2N5, file="simulateb2N5.Rdata")

# Plot Supp Fig 4 (b)    
set.seed(6326)
pl <- plotHyperparlam(simulateb2N5, plotlog=TRUE, methinclude = c("MML_mgcv", "GCV_mgcv", "CV_glmnet"))


#######################################################################################
########################## BATCH AND PARALLEL COMPUTING ###############################
#######################################################################################

# R CMD BATCH HyperparEst_Simulationmaster.R output.txt &

# Empty work space
rm(list=ls())

# Initializing snowfall:
library(snowfall)
ncpus = 8
sfInit(parallel = TRUE, cpus = ncpus)

sfLibrary(mvtnorm)
sfLibrary(glmnet)
sfLibrary(HiLMM)
sfLibrary(corpcor)

# Function needed
source("HyperparEst_linearSource.R") #source for simulating linear predictor

sfExportAll()

# number of simulated data sets
nr.iter = 100

############
# Simulation 2: independent X, large setting, n=1000, p=15000, gaussian prior beta
pmt <- proc.time(); set.seed(6326) 
simulate2 = sfLapply(1:nr.iter, function(i){
  wrapper_basic(tau2_true = 0.01, sigma2_true = 150, n = 1000, p = 15000, p0 = 0, cor = 0, data = F)})
proc.time()-pmt
save(simulate2, file="simulate2.Rdata")

############
# Simulation 8: Data covariates large, response simulated. KIRC data (n=71, p=18391)
load("KIRCdat.Rdata")
X <- t(KIRCdat)
tau2_true = 0.01; p<- ncol(X)
sigma2_true = p*tau2_true  #such that heritability =0.5
sigma2_true
sfExport("sigma2_true", "X")
pmt <- proc.time(); set.seed(6326) 
simulate8 = sfLapply(1:nr.iter, function(i){print(paste("ITERATION",i));
  wrapper_basic(tau2_true = 0.01, sigma2_true = sigma2_true, data = X)})
proc.time()-pmt
save(simulate8, file="simulate8.Rdata")
load("simulate8.Rdata")
pl <- plotHyperpar(simulate8)

