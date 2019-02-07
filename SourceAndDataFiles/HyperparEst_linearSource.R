########################################################
####Estimation of error and prior variance functions####
########################################################

########################################################
######################packages##########################
########################################################
# library(nlme)
# library(mvtnorm)
# library(HiLMM)
# library(glmnet)
# library(reshape2)
# library(ggplot2)
# library(gridExtra)

########################################################
##################Functions#############################
########################################################

#source("C:/Users/jurre/Documents/StatisticalScience/Year2/MasterThesis/Simulations/linRegMM.R")


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

SimpleMoM = function(Y, X){
  Y <- matrix(Y,nrow=1)
  out <- t(Y) %*% Y
  
  Sigma <- X %*% t(X)
  
  #sum(Sigma)-sum(diag(Sigma))
  
  
  tau2 <- (sum(out)-sum(diag(out)))/(sum(Sigma)-sum(diag(Sigma)))
  
  
  if(tau2<0) tau2 <- 0
  
  sigmaest <- max(0,mean(((as.numeric(Y)-mean(Y))^2 - diag(Sigma)*tau2)))
  
  return(c(sigmaest, tau2))
}

tauest <- function(X, Y, optl=1,innfold=10,compareCV=TRUE,trues2=NULL,SVD=NULL){
  #X=X;Y=Y;optl=1;innfold=10;compareCV=FALSE
  nsam <- nrow(X)
  if(compareCV){
    sd_y <- sqrt(var(Y)*(nsam-1)/nsam)
    # glmnlam1 <- cv.glmnet(X,Y,alpha=0,nfolds=innfolds[1],standardize=FALSE,intercept=FALSE,lambda.min.ratio=0.0001)
    # glmntaus1 <- (glmnlam1$lambda.min*(nsam/sd_y))^(-1)
    pmt <- proc.time()
    glmnlam2 <- cv.glmnet(X,Y,alpha=0,nfolds=innfold,standardize=FALSE,intercept=FALSE,lambda.min.ratio=0.0001)
    pmt1 <- proc.time()-pmt
    pmt1
    # penlam2 <- optL2(Y,X,unpenalized= ~0,fold=innfold,standardize=FALSE)
    glmntaus2 <- (glmnlam2$lambda.min*(nsam/sd_y))^(-1)
  }
  
  
  pmt <- proc.time()
  nfeat <- ncol(X)
  if(is.null(SVD)) SVD <- svd(X)
  constlam <- 1
  leftmat <- SVD$v %*% diag(1/((SVD$d)^2+constlam*optl)) %*% diag(SVD$d) %*% t(SVD$u)
  coeffsvd <- leftmat %*% Y
  coeff2svd <- (coeffsvd)^2
  #  coeff2 <- (solve(t(X) %*% X + diag(rep(optl,nfeat))) %*% t(X) %*% Y)^2
  if(is.null(trues2)){
    Hatm <- X %*% leftmat 
    df <- nsam - sum(diag(2*Hatm - Hatm %*% t(Hatm)))
    preds <- X %*% coeffsvd
    VarRes <- sum((Y - preds)^2)/df
    print(paste("Sigma^2 estimate:",VarRes))
    vars3 <- VarRes * rowSums(leftmat^2)
  } else  {vars3 <- trues2* rowSums(leftmat^2)}
  cii2 <- (rowSums(leftmat * t(X)))^2
  leftmatsq <- leftmat/sqrt(vars3)
  lefts2 <-  sum(coeff2svd/vars3)-nfeat
  #rightmat =  X
  cfm <-  sum((t(leftmatsq) %*% leftmatsq) * (X %*% t(X)))
  cfm2 <-  sum((t(leftmat) %*% leftmat) * (X %*% t(X)))
  
  est1 <- lefts2/cfm
  pmt2 <- proc.time()-pmt
  pmt2
  #est2 <- (sum(coeff2svd)-sum(vars3))/nfeat
  est2 <- (sum(coeff2svd)-sum(vars3))/cfm2
  if(compareCV) allest <- c(cv10fold=glmntaus2,estvdW=est1,estsimple=est2) else {
    allest <- c(estvdW=est1,estsimple=est2)
  }
  return(allest)
}

mylinearRidge = function (Y, X,  lambda = "automatic", nPCs = NULL, SVD=NULL, scaling ="corrForm") 
{
 # Y <- Y ; X = X; SVD = svd(X); lambda = "automatic";nPCs=NULL;scaling="corrForm"
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(SVD)) Xs  <- svd(X) else Xs <- SVD
  Q <- Xs$v
  Z <- X %*% Q
  Lambda <- Xs$d^2
  if (!is.null(lambda) && lambda == "automatic") {
    automatic <- TRUE
    if (is.null(nPCs)) {
      propVar <- cumsum(Lambda)/sum(Lambda) * 100
      nPCs <- which.max(propVar[propVar < 90]) + 1
    }
  }
  ahat <- diag(1/Lambda) %*% t(Z) %*% Y
  if (!is.null(lambda) && lambda == "automatic" && !is.null(nPCs)) {
    ks.vector <- sig2hat.vector <- vec.df <- numeric(nPCs)
    flag <- TRUE
    P <- 0
    while ((P < nPCs) && flag) {
      P <- P + 1
      sig2hat <- ifelse(P == 1, as.numeric(crossprod(Y - 
                                                       (Z[, 1]) * ahat[1])/(n - 1)), as.numeric(crossprod(Y - 
                                                                                                            (Z[, 1:P]) %*% ahat[1:P])/(n - P)))
      ahatsum <- ifelse(P == 1, ahat[1]^2, sum(ahat[1:P]^2))
      ks.vector[P] <- P * sig2hat/ahatsum
      if (is.finite(ks.vector[P])) {
        vec.df[P] <- sum(Lambda^2/(Lambda + ks.vector[P])^2)
      }
      if (!is.finite(ks.vector[P])) {
        flag <- FALSE
        ks.vector <- ks.vector[1:(P - 1)]
      }
    }
    nPCs.dof <- which.min(abs(vec.df - seq(nPCs)))
    lambda <- ks.vector
    chosen.nPCs <- nPCs.dof
    max.nPCs <- nPCs
  }
  else if (!is.null(nPCs)) {
    P <- nPCs
    sig2hat <- ifelse(P == 1, as.numeric(crossprod(Y - (Z[, 
                                                          1]) * ahat[1])/(n - 1)), as.numeric(crossprod(Y - 
                                                                                                          (Z[, 1:P]) %*% ahat[1:P])/(n - P)))
    ahatsum <- ifelse(P == 1, ahat[1]^2, sum(ahat[1:P]^2))
    lambda <- P * sig2hat/ahatsum
    chosen.nPCs <- nPCs
  }

  res <- list(lambda=lambda, sig2hat = sig2hat)
  if (!is.null(nPCs)) {
    if (automatic) {
      res$max.nPCs <- max.nPCs
    }
    res$chosen.nPCs <- chosen.nPCs
  }
  res
}



MargLikelihood = function(ts, XXT, Y){
  logtausq<-ts[1] ; logsigmasq<-ts[2]; n = nrow(XXT)
  varY = XXT * exp(logtausq) + diag(rep(1,n))*exp(logsigmasq)
  mlk <- -dmvnorm(Y,mean=rep(0,n),sigma=varY,log=TRUE)
  return(mlk)
}

# INPUT:
# y : n-dimensional response vector
# X : n by p design matrix
# plotML : logical (to plot log ML as a function of hyperparameter)
# a : Gamma shape prior parameter for error noise
# b : Gamma rate prior parameter for error noise
BayesRidge <- function(y, X, plotML=FALSE, a=0.001, b=0.001, SVD=NULL){
  #y <- Y; X<-X; a=0.001; b=0.001
  ## SVD decomposition
  if(is.null(SVD)) mysvd <- fast.svd(X) else mysvd <- SVD
  eigval <- mysvd$d^2
  
  ## Needed to calculate log-ML
  alphahat <- diag(1/mysvd$d)%*%t(mysvd$u)%*%y
  q <- length(eigval)
  yTy <- as.numeric(crossprod(y, y))
  
  ## log-ML function
  logML <- function(lambda, q, di, a, b, n, ahat, yy){
    part1 <- (q/2)*log(lambda)-0.5*sum(log(lambda+di))-(a+n/2)*log(b+0.5*yy-0.5*sum(((ahat^2)*(di^2))/(di+lambda)))
    part2 <- a*log(b) + lgamma(a+n/2) - lgamma(a) - (n/2)*log(pi)
    part1 + part2
  }
  myF <- function(lambda, q, di, a, b, n, ahat, yy){
    part1 <- 0.5*sum(log(lambda))-0.5*sum(log(lambda+di))-(a+n/2)*log(b+0.5*yy-0.5*sum(((ahat^2)*(di^2))/(di+lambda)))
    part2 <- a*log(b) + lgamma(a+n/2) - lgamma(a) - (n/2)*log(pi)
    -(part1 + part2)
  }
  
  ## Determine the log-ML optimal lambda
  # Step 1: find interval that contain global maximum
  lambdas <- exp(seq(from=-5, to=20, length=200))
  mygrid <- matrix(NA, length(lambdas), 2)
  colnames(mygrid) <- c("lambda", "logML")
  mygrid[,1] <- lambdas
  vals <- sapply(lambdas, logML, q=q, di=eigval, a=a, b=b, n=length(y), ahat=alphahat, yy=yTy)
  mygrid[,2] <- vals
  idx <- which.max(vals)
  
  if(idx==1){
    lb <- 1e-6
  }else{
    lb <- lambdas[idx-1]
  }
  if(idx==length(lambdas)){
    ub <- lambdas[idx]*2
  }else{
    ub <- lambdas[idx+1]
  }
  # Step 2: Optimization routine to find the exact global maximizer
  opt <- optimize(f = logML, interval = c(lb,ub), maximum = TRUE, q=q, di=eigval, a=a, b=b, n=length(y), ahat=alphahat, yy=yTy)
  
  # Get Optimum
  optLambda <- opt$maximum
  
  # Plot logML
  if(plotML){
    plot(log(lambdas), mygrid[,2], type="l", ylab="log marginal likelihood", xlab="log lambda")
    abline(v=log(optLambda), col="red")
  }
  
  # Calculate posterior means and variances
  #postMean <- as.vector(mysvd$v %*% solve(diag(q)+optLambda*diag(1/mysvd$d^2)) %*% alphahat)
  postsiga <-  (a+length(y)/2)
  postsigb <- 0.5*(yTy-sum(((alphahat^2)*eigval^2)/(eigval+optLambda)))
  # sigmaminus2hat <- postsiga/postsigb
  # varalpha <- sigmaminus2hat*solve(diag(eigval)+optLambda*diag(q))
  # myf <- function(ll, mm1, mm2) mm1[ll,]%*%mm2[,ll]
  # mymm1 <- mysvd$v %*% varalpha
  # mymm2 <- t(mysvd$v)
  # postSd <- sqrt(sapply(1:nrow(mysvd$v), myf, mm1=mymm1, mm2=mymm2))
  
  # Output
  out <- list("postsig"= c(postsiga,postsigb),"optLambda"=optLambda)
  return(out)
}



MomentsMethod = function(Y, X, lambda, SVD=NULL){
  if(is.null(SVD)) SVD.x = svd(X) else SVD.x <- SVD
  R = SVD.x$u %*% diag(SVD.x$d)
  solve. = solve(t(R) %*% R + lambda * diag(min(ncol(X), nrow(X))))
  
  H = R %*% solve. %*% t(R)
  
  C = SVD.x$v %*% solve(t(R) %*% R + lambda * diag(min(ncol(X), nrow(X)))) %*% t(R) 
  sigma_tildejj = rowSums(C*C)
  # left.sigma = SVD.x$v %*% solve(t(R)%*%R + lambda * diag(min(ncol(X), nrow(X)))) %*% t(R) %*% 
  #   R %*% solve(t(R)%*%R + diag(min(ncol(X), nrow(X))))
  # right.sigma = t(SVD.x$v)
  #sigma_tilde =  rowSums(left.sigma * t(right.sigma)) #diagonal of \tilde{\Sigma}
  
  L.star = SVD.x$v %*% solve. %*% t(R)
  L = (1/sqrt(sigma_tildejj)) * L.star
  sum.C2 = sum((t(L) %*% L) * (X %*% t(X)))
  
  #Start using Y here:
  beta = C %*% Y
  Xbeta = H %*% Y
  p = ncol(X)
  
  a1 = t(Y - Xbeta) %*% (Y - Xbeta)
  a2 = sum(beta^2/sigma_tildejj)/sum.C2
  c1 = length(Y) + sum(H * t(H)) - 2*sum(diag(H))
  c2 = sum((t(X) * t(X))) + sum((t(X) %*% H %*% H) * t(X)) - 2 * sum((t(X) %*% H) * t(X))
  d1 = p/sum.C2
  d2 = 1
  
  A = matrix(c(c1, c2, d1, d2), nrow = 2, byrow = T)
  b = matrix(c(a1, a2), nrow = 2)
  Estimates = solve(A, b)
  
  if(Estimates[1] < 0){  return(c(0, a2))}
  else(if(Estimates[2] <0) {return(c(a1/c1 , 0))}
       else(return(Estimates)))
}

SimpleEstimate = function(response=data$y, X = as.matrix(data[,-1]), lambda = 1, SVD=NULL){
  
  if(is.null(SVD)) SVD <- svd(X)
  leftmat <- SVD$v %*% diag(1/((SVD$d)^2+1*lambda)) %*% diag(SVD$d) %*% t(SVD$u)
  n = length(response)
  Hatm <- X%*% leftmat 
  preds = Hatm %*% response
  df <- n - sum(diag(2*Hatm - Hatm %*% t(Hatm)))
  VarRes <- sum((response - preds)^2)/df
  
  return(VarRes)
}

SimpleEstimate2 = function(response=data$y, X = as.matrix(data[,-1]), lambda = 1, SVD=NULL){
  if(is.null(SVD)) SVD <- svd(X)
  leftmat <- SVD$v %*% diag(1/((SVD$d)^2+1*lambda)) %*% diag(SVD$d) %*% t(SVD$u)
  n = length(response)
  Hatm <- X%*% leftmat 
  preds = Hatm %*% response
  df <- n - sum(diag(Hatm))
  VarRes <- sum((response - preds)^2)/df
  
  return(VarRes)
}

Golub_V = function(lambda, Y, X, SVD=NULL){
  #Equation 1.4 from Golub, Heath & Wahba 1979
  Y = as.matrix(Y)
  n = nrow(X); p = ncol(X)
  
  if(is.null(SVD)) SVD = svd(X)
  R = SVD$u %*% diag(SVD$d); V = SVD$v
  A = R %*% solve(t(R) %*% R + n * lambda * diag(min(n,p))) %*% t(R)
  
  V = 1/n * t((diag(n) - A)%*%Y) %*% ((diag(n) - A)%*%Y) /
    (1/n * sum(diag((diag(n) - A))))^2
  return(V)
}

Simulatedata = function(tau2_true, sigma2_true, n, p, p0, cor, pblock=10, data, simmodel = c("normal","laplace", "uniform", "t4error")){
  #data = F if no data is used otherwise add a matrix with p columns and n rows to be used as X matrix.
  #tau2_true = true value of beta prior variance
  #sigma2_true = true value of error variance
  #n = number of samples
  #p = number of variables
  #cor = correlation coefficient. cor = 0 means "no" correlations in X matrix.
  #p0 = Proportion of null coefficients. if p0 > 0, a spike and slab prior is used
  #simmodel: simulation model
  if(is.vector(simmodel)) simmodel <- simmodel[1]
  
  if(length(dim(data)) != 2){ #Check for data = F
    
    # if(cor == 0){
    #   x = matrix(rnorm(n*p, 0, 1), ncol = p)
    #   x <- t((t(x) - apply(t(x),1,mean))/apply(t(x),1,sd))
    #             } #if no correlation..
    # else{
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    Nblock <- p/pblock
    if(!is.wholenumber(Nblock)){
      print("ERROR: p should be a multiple of pblock")
      return(NULL) 
    }
    
    if(cor==0) x <- matrix(rnorm(n*p),n,p) else x <- Reduce(cbind,lapply(1:Nblock, function(z) matrix(rep(rnorm(n,sd=sqrt(cor/(1-cor))),times=pblock),n,pblock))) + matrix(rnorm(n*p),n,p)
    x <- t((t(x) - apply(t(x),1,mean))/apply(t(x),1,sd))
    if(cor != 0) print(paste("Simulating covariate data with cor = ",cor," and block size = ", pblock, sep = "")) else
      print(paste("Simulating independent covariate data"))
    
  }
  else{x = scale(data, center = TRUE, scale = TRUE) #If there is data supplemented..
  n = nrow(data)
  p = ncol(data)
  print("Simulating using known covariates")
  }
  if((simmodel=="normal" && p0==0) || simmodel=="t4error"){
    print(paste("Simulating betas from N(0,",tau2_true,")",sep=""))
    beta = matrix(rnorm(p, 0, sqrt(tau2_true)), ncol = 1, nrow = p) 
  }
  if(simmodel=="laplace"){
    #var = 2b^2
    
    lambda1 <- sqrt(tau2_true/2)  
    print(paste("Simulating betas from laplace(0,",signif(lambda1,3),")",sep=""))
    #A Laplace(0,b) variate can also be generated as the difference of two i.i.d. 
    #Exponential(1/b) random variables
    beta <-   matrix(rexp(p, 1/lambda1), ncol = 1, nrow = p)  -  matrix(rexp(p, 1/lambda1), ncol = 1, nrow = p)
  }
  if(simmodel=="uniform"){
    #print("Simulating betas from Uniform")
    #Var = (b-a)^2/12 => b-a = sqrt(12*tau2_true) => b = (b-a)/2
    b= sqrt(12*tau2_true)/2
    a =-b
    print(paste("Simulating betas from Uniform(",signif(a,3),",",signif(b,3),")",sep=""))
    beta <-   matrix(runif(p, min=a,max=b), ncol = 1, nrow = p) 
  }
  
  if(simmodel=="normal" & p0 > 0){ #adding spike and slab prior
    tau0 = tau2_true/(1-p0)
    print(paste("Simulating betas from ",p0,"*Spike0+",(1-p0),"*N(0,",signif(tau0,3),")",sep=""))
    beta = ifelse(runif(p) < p0, 0, rnorm(p, 0, sqrt(tau0)))
  }
  
  if(simmodel !="t4error") {
    print("Simulate Gaussian error")
    error = rnorm(n, 0, sqrt(sigma2_true))
  } else {
    mult <- sqrt(sigma2_true/(4/(4-2)))
    error = mult*rt(n,df=4)
    print(paste("Simulate scaled t4 errors. Errors scaled to have variance =", sigma2_true))
  }
  y = x %*% beta + error
  
  #Adding colnames
  data = data.frame(y, x)
  colnames(data)<- paste("V", 0:p, sep = "")
  colnames(data)[1] <- "y"
  return(list("data" = data, "beta" = beta))
}

wrapper_basic = function(tau2_true, sigma2_true, n=NULL, p=NULL, p0=0, cor=0, pblock=10, data=F, simmodel = c("normal","laplace", "uniform", "t4error")){
  ##High Dimensional Wrapper here p may be larger than n  
  #FALSE if data has to be simulated.
  #tau2_true = true prior variance
  #sigma2_true = true error variance
  #n = sample size
  #p = number of variables
  #p0 is proportion of null effects in variables (for sparse simulation)
  #lambda = initial lambda value for MoM and Simple Estimate
  #data = either a two dimensional matrix as dataset where the first column is the y variable and the other columns are X variables, or..
  #tau2_true = 0.1;sigma2_true = 1; n = 50;p = 200;p0 = 0;cor = 0;data = F
  print("Simulating data...")
  if(is.null(n)) {
   n=nrow(data)
   p=ncol(data)
  }
  MCdata = Simulatedata(tau2_true, sigma2_true, n = n, p = p, p0 = p0, cor, pblock=pblock, data = data, simmodel=simmodel)
  cts <- c()
  data = MCdata$data
  
  Y = data$y
  X = as.matrix(data[,-1])
  
  print("SVD...")
  pmt<- proc.time()
  SVDx <- svd(X)
  ctsvd <- (proc.time()-pmt)[3]
  
  #Empirical bayes marginal likelihood maximization
  print("Finding Maximum Marginal Likelihood...")
  X = as.matrix(data[,-1])
  pmt<- proc.time()
  XXT = X %*% t(X)
  optimal <- optim(par = c(log(0.01),log(10)), fn = MargLikelihood, XXT = XXT, Y = data$y)
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  
  print("Conjugate Bayes Estimator...")
  pmt<- proc.time()
  BR <- BayesRidge(Y,X,a=0.001,b=0.001,SVD=SVDx)
  ct <- (proc.time()-pmt)[3] + ctsvd
  cts <- c(cts,ct)
  lamB <- BR$optLambda
  sigmaB <- BR$postsig[2]/(BR$postsig[1]+1)
  tauB <- sigmaB/lamB
  
  #MoM Estimates
  print("Estimating Moments Estimate...")
  pmt<- proc.time()
  #MoM.old = MomentsMethod(Y = Y, X = X,lambda = lambda)
  MoM = SimpleMoM(X = X, Y = Y)
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  
  
  #MoM Estimates 2
  # print("Estimating Moments Estimate 2...")
  # pmt<- proc.time()
  # MoM.old = MomentsMethod(Y = Y, X = X,lambda = lambda, SVD=SVDx)
  # #MoM = SimpleMoM(X = X, Y = Y)
  # ct <- (proc.time()-pmt)[3] + ctsvd
  # cts <- c(cts,ct)
  
  #Cules sigma2 estimate:
  print("Calculating PCR estimate...")
  pmt<- proc.time()
  CULEvar = mylinearRidge(Y = Y, X = X, lambda = "automatic", SVD = SVDx)
  ct <- (proc.time()-pmt)[3] + ctsvd
  cts <- c(cts,ct)
  
  #Hasties sigma2 estimate:
  print("Calculating Basic estimate...")
  pmt<- proc.time()
  VarRes = SimpleEstimate(response = Y, X = X, lambda = 1, SVD=SVDx)
  ct <- (proc.time()-pmt)[3] + ctsvd
  cts <- c(cts,ct)
  
  #GrRidge tauest
  #print("Calculating GRridge tau estimate...")
  #grSim = tauest(X = X, Y = Y, compareCV = F, SVD=NULL)
  
  #Jiang estimate
  #print("calculating REML estimate")
  #REMLsigma = tryCatch(linRegMM(X = X, y = Y), error=function(e){return(data.frame(sb2 = NA, sy2 = NA))})
  # Error in if (abs(lb0[iter] - lb0[iter - 1]) < tol) { : 
  #     missing value where TRUE/FALSE needed --> sy2 goes to infinity and sb2 NaN
  
  #GCV
  print("Estimating lambda with Golub")
  pmt<- proc.time()
  Golub = n*optim(1, Golub_V,  method = "Brent", lower = 0, upper = ((sigma2_true/tau2_true*20)/n), X = X, Y = Y, SVD=SVDx)$par
  ct <- (proc.time()-pmt)[3] + ctsvd
  cts <- c(cts,ct)
  
  
  
  #GLMnet
  print("Estimating lambda with glmnet")
  pmt<- proc.time()
  n = length(Y)
  sy = sqrt(var(Y)*(n-1)/n)
  #if(laplace) alphaglm <- 1 else alphaglm = 0
  alphaglm <- 0
  CV = cv.glmnet(x = X, y = Y, alpha = alphaglm, nfolds = 10, grouped = FALSE, standardize = F)$lambda.min * n / sy
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  
  #EstHer
  print("Estimating Heritability")
  pmt<- proc.time()
  EH = estim_herit(Y,X)$heritability
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  
  #Making output data.frame
  output = data.frame("tau2" = c(tau2_true,
    exp(optimal$par[1]), tauB, MoM[2], 
    #MoM.old[2],
    #REMLsigma$sb2, 
    NA, NA, 
    NA, 
    NA, 
    NA), 
    "sigma2" = c(sigma2_true, 
      exp(optimal$par[2]), sigmaB, MoM[1],
      #MoM.old[1],
      #REMLsigma$sy2, 
      CULEvar$sig2hat, VarRes,  
      NA, 
      NA, 
      NA), 
    "lambda" = c(sigma2_true/tau2_true,exp(optimal$par[2])/exp(optimal$par[1]),
    sigmaB/tauB,MoM[1]/MoM[2],
    #MoM.old[1]/MoM.old[2],
    NA,NA,Golub,CV,NA),
   
     "heritability"= c(p*tau2_true/(sigma2_true+p*tau2_true),
                       p*exp(optimal$par[1])/(p*exp(optimal$par[1]) + exp(optimal$par[2])),
                                              p*tauB/(sigmaB+p*tauB),
                                              p*MoM[2]/(p*MoM[2]+ MoM[1]),
                       #p*MoM.old[2]/(p*MoM.old[2]+ MoM.old[1]),
                                              NA,NA,NA,NA,EH),
    "computing times" = c(0,cts)
  )
   rownames(output) = c("true", 
                        "MML", "Bayes", "MoM", #Joint estimation
                        #"MoMbeta",  
                        #"REML", 
                        "PCR", "Basic",        #Sigma Estimation
                        "GCV",       #Lambda estimation
                        "CV",        
                        "HiLMM")               #Heritability estimation
  return(output)
}

#functions for mixed effects case


Coef_Optim = function(ts, X, Y,nr.fixed=5){
  #Also estimating coefficients
  logtausq<-ts[1] ; logsigmasq<-ts[2];
  beta = ts[3:(nr.fixed+3)]
  n = nrow(X)
  X_fixed = cbind(1,X[,1:nr.fixed]); X_random = X[,-(1:nr.fixed)]
  mean_Y = as.numeric(X_fixed %*% beta)
  varY = X_random %*% t(X_random) * exp(logtausq) + diag(rep(1,n))*exp(logsigmasq)
  mlk <- -dmvnorm(Y,mean= mean_Y,sigma=varY,log=TRUE)
  return(mlk)
}

MML_step = function(ts, X, Y, nr.fixed, coeff){
  if(names(coeff)[1] == "(Intercept)"){check1 = (length(coeff)-1 < nr.fixed)}
  else(check1 = length(coeff) < nr.fixed)
  if(check1){
    if(names(coeff)[1] == "(Intercept)"){X_fixed = as.matrix(cbind(1, as.data.frame(X)[,names(coeff)[-1]]))}
    else{X_fixed = as.data.frame(X)[,names(coeff)]}
  }
  else{X_fixed = cbind(1,X[,1:nr.fixed])}
  #Putting all estimated fixed coefficients into the model
  logtausq<-ts[1] ; logsigmasq<-ts[2]
  n = nrow(X)
  X_random = X[,-(1:nr.fixed)]
  mean_Y = as.numeric(X_fixed %*% coeff)
  varY = X_random %*% t(X_random) * exp(logtausq) + diag(rep(1,n))*exp(logsigmasq)
  mlk <- -dmvnorm(Y,mean= mean_Y,sigma=varY,log=TRUE)
  return(mlk)
}

MargLikelihood_REML= function(ts, X, Y,q=5, coeff){
  #Putting all estimated fixed coefficients into the model and using REML
  logtausq<-ts[1] ; logsigmasq<-ts[2]
  n = nrow(X)
  X_fixed = as.matrix(cbind(1,X[,1:q]));X_random = X[,-(1:q)]
  residual = Y - X_fixed %*% solve(t(X_fixed) %*% X_fixed) %*% t(X_fixed) %*% Y
  #mean_Y = as.numeric(X_fixed %*% coeff)
  varY = X_random %*% t(X_random) * exp(logtausq) + diag(rep(1,n))*exp(logsigmasq)
  mlk <- -dmvnorm(as.numeric(residual),mean= rep(0 ,n),sigma=varY,log=TRUE)
  return(mlk)
}

REMML = function(parameter, Y, X, fix.index){
  logsigma2 = parameter[2]; logtau2 = parameter[1]
  
  k = length(fix.index); n = length(Y)
  X.f = X[,fix.index]
  X.r = X[,-(fix.index)]
  
  
  COV = X.r %*% t(X.r) * exp(logtau2) + exp(logsigma2) * diag(nrow(X.r))
  COV.s = solve(COV)
  
  beta = solve(t(X.f) %*% COV.s %*% X.f) %*% t(X.f) %*% COV.s %*% Y
  
  loglikelihood = #-0.5 * (n - k) * log(2 * pi) + 0.5 * log(det(t(X) %*% X)) - #This line is all constants.
    -log(det(COV)) - log(det(t(X.f) %*% COV.s %*% X.f)) -
    t(Y - X.f %*% beta) %*% COV.s %*% (Y - X.f %*% beta)
  return(-loglikelihood)}

Simulatedata_fixedeffects = function(tau2_true, sigma2_true, tau2_truefixed, n = 1000, p = 100, p0=0.9, p0fixed=0.5, m = nr.fixed){
  x = matrix(rnorm(n*(p+m), 0, 1), ncol = (p+m))
  x = t((t(x) - apply(t(x),1,mean))/apply(t(x),1,sd))
  error = rnorm(n, 0, sqrt(sigma2_true))
  
  
  tau0 = tau2_true/(1-p0)
  print(paste("Simulating betas from ",p0,"*Spike0+",(1-p0),"*N(0,",signif(tau0,3),")",sep=""))
  beta = ifelse(runif(p) < p0, 0, rnorm(p, 0, sqrt(tau0)))
  tau0_2 = tau2_truefixed/(1-p0fixed)
  print(paste("Simulating fixed effects from ",p0fixed,"*Spike0+",p0fixed,"*N(0,",signif(tau0_2,3),")",sep=""))
  alpha = ifelse(runif(m) < p0fixed, 0, rnorm(m, sqrt(tau0_2)))
  coef = c(alpha, beta)
  y = x %*% coef + error
  
  data = data.frame(y, x)
  data$y = as.numeric(y)
  colnames(data)<- paste("V", 0:(p+m), sep = "")
  colnames(data)[1] <- "y"
  return(list("data" = data, beta = beta, alpha = alpha))
}

wrapper_mixed = function(tau2_true = 0.01, sigma2_true = 10, tau2_truefixed= 5*tau2_true, n = 100, p = 1000, p0=0.9, p0fixed=0.5,nr.fixed=10){
  # set.seed(19)
  # print(paste0("seed", i))
  dataset = Simulatedata_fixedeffects(tau2_true = tau2_true, sigma2_true = sigma2_true, tau2_truefixed=tau2_truefixed, n = n, p = p, p0=0.9, 
                                      p0fixed=0.5, m= nr.fixed)
  Data = as.matrix(dataset$data)
  cts <- c()
  
  
  print("OptimCoef")
  #fit fixed effect model first to obtain initial values for optimization
  form_fixed = "Data[,1] ~ V1"
  for(i in 2:nr.fixed){
    form_fixed = paste(form_fixed, " + V", i, sep = "")
  }
  form_fixed = formula(form_fixed)
  lmo = lm(form_fixed, data = as.data.frame(Data))

  coeff = step(lmo, verbose = F,trace=0)$coefficients
  pmt <- proc.time()
  estimates = optim(par = c(log(0.01),log(10), lmo$coef), fn = Coef_Optim, X = Data[,-1], Y = Data[,1], nr.fixed = nr.fixed)
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  # print("MML stepwise")
  # estimates_2 = optim(par = c(log(0.01),log(10)), fn = MML_step, X = Data[,-1], Y = Data[,1] ,coeff = coeff, nr.fixed = nr.fixed)
  # print("Regular MML")
  # pmt <- proc.time()
  # estimates_3 = optim(par = c(log(0.01),log(10)), fn = MargLikelihood, X = Data[,-1], Y = Data[,1])
  # ct <- (proc.time()-pmt)[3]
  # cts <- c(cts,ct)
  
  print("REML")
  pmt <- proc.time()
  REMML = optim(par = c(log(0.01),log(10)), REMML, Y = Data[,1], X = Data[,-1], fix.index = 1:10)
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  # print("Jiang")
  # Jiang = tryCatch(linRegMM(X = Data[,-(1:(nr.fixed+1))], Z = Data[, 2:(nr.fixed+1)],y = Data[,1]), error=function(e){return(data.frame(sb2 = NA, sy2 = NA))})
  # 
  #CV = cv.glmnet(x = Data[,-(1:(nr.fixed+1))], y = Data[,1], alpha = 0, nfolds = nrow(Data), grouped = FALSE)
   
  # pmt <- proc.time()
  # print("GCV")
  # Golub = optim(1, Golub_V,  method = "Brent", lower = 0, upper = ((sigma2_true/tau2_true*20)/nrow(Data)), X = Data[,-(1:(nr.fixed+1))], Y = Data[,1])$par
  # ct <- (proc.time()-pmt)[3]
  # cts <- c(cts,ct)
  
  print("CV")
  Y = Data[,1]; n = length(Y); pplusm <- ncol(Data)-1 #first col = response
  pf <- c(rep(0,nr.fixed),rep(1,pplusm-nr.fixed)) #penalty factor = 0 for fixed effects
  sy = sqrt(var(Y)*(n-1)/n)
  pmt <- proc.time()
  CV = cv.glmnet(x = Data[,-1], y = Data[,1], alpha = 0 ,nfolds = 10, penalty.factor=pf, standardize = F, grouped = FALSE)$lambda.min * length(Y) / sy
  ct <- (proc.time()-pmt)[3]
  cts <- c(cts,ct)
  
  #EH = Selvar(Data[,1], Z = Data[,-(1:(nr.fixed+1))], X = Data[, 2:(nr.fixed+1)], thresh_vect=c(0.7))$heritability
  # pmt <- proc.time()
  # print("HiLMM")
  # EH = estim_herit(Y,Data[,-1])$heritability
  # ct <- (proc.time()-pmt)[3]
  # cts <- c(cts,ct)
  # print("done")
  
  
  output <- data.frame("tau2" = c(tau2_true, 
                                  exp(estimates$par[1]), #MML
                                  #exp(estimates_3$par[1]), 
                        exp(REMML$par[1]), NA), 
             "sigma2" = c(sigma2_true, 
                         exp(estimates$par[2]),
                         # exp(estimates_3$par[2]),
                        exp(REMML$par[2]),NA),
             "lambda"= c(sigma2_true/tau2_true, 
                         exp(estimates$par[2])/exp(estimates$par[1]),
                        # exp(estimates_3$par[2])/exp(estimates_3$par[1]),
                         exp(REMML$par[2])/exp(REMML$par[1]),CV),
             "heritability"= c(p*tau2_true/(sigma2_true + p*tau2_true),
                             p*exp(estimates$par[1])/(exp(estimates$par[2])+ p*exp(estimates$par[1])),
                               #p*exp(estimates_3$par[1])/(exp(estimates_3$par[2])+ p*exp(estimates_3$par[1])),
                               p*exp(REMML$par[1])/(exp(REMML$par[2])+p*exp(REMML$par[1])),
                               NA), 
             "computing times" = c(0,cts)
    )
  rownames(output) = c("true", 
                       "MML", "REML", #Joint estimation
                       #"MoMbeta",  
                       #"REML", 
                       "CV"       #Lambda estimation
                       )               #Heritability estimation
  
  
  return(output)
  # return(data.frame("Optim" = c(estimates$value, exp(estimates$par[1:2])), 
  #              #     "STEP" = c(estimates_2$value, exp(estimates_2$par)), 
  #                   "MML" = c(estimates_3$value, exp(estimates_3$par)),
  #                   "REMML" = c(REMML$value, exp(REMML$par)),
  #             #      "Jiang" = c(NA, Jiang$sb2, Jiang$sy2),
  #                   "CV" = CV,
  #                   "GCV" = Golub,
  #                   "Herit" = EH
  # ))
}

wrapper_basicCV = function(tau2_true, sigma2_true, n, p, p0, cor, lambda = 1, data = F, folds=c(10,n)){  #function to compare several CV results
  ##High Dimensional Wrapper here p may be larger than n  
  #FALSE if data has to be simulated.
  #tau2_true = true prior variance
  #sigma2_true = true error variance
  #n = sample size
  #p = number of variables
  #p0 is proportion of null effects in variables (for sparse simulation)
  #lambda = initial lambda value for MoM and Simple Estimate
  #data = either a two dimensional matrix as dataset where the first column is the y variable and the other columns are X variables, or..
  #tau2_true = 0.1;sigma2_true = 1; n = 50;p = 200;p0 = 0;cor = 0;data = F
  print("Simulating data...")
  MCdata = Simulatedata(tau2_true, sigma2_true, n = n, p = p, p0 = p0, cor, data = data)
  data = MCdata$data
  
  Y = data$y
  X = as.matrix(data[,-1])
  
  
  
  #GLMnet
  cts <- CVs <- c()
  for(k in folds){
  print(paste("Estimating lambda with glmnet, fold=", k))
  pmt<- proc.time()
  n = length(Y)
  sy = sqrt(var(Y)*(n-1)/n)
  #if(laplace) alphaglm <- 1 else alphaglm = 0
  alphaglm <- 0
  CV = cv.glmnet(x = X, y = Y, alpha = alphaglm, nfolds = k, grouped = FALSE, standardize = F)$lambda.min * n / sy
  ct <- (proc.time()-pmt)[3]
  CVs <- c(CVs,CV)
  cts <- c(cts,ct)
  }
  lamtims <- cbind(CVs,cts)
  rownames(lamtims) <- paste(as.character(folds),"-fold",sep="")
  return(lamtims)
}


