##################### source function for binomial ridge regression ############
wrapper_binomial <- function(lambda_true,n,p,p0,cor,data=F, N=5, penalized=FALSE){
  #lambda_true=200;n=200;p=2000; p0=0;cor=0; data=F; N=5;penalized=FALSE
  simdata <- Simulatedata(tau2_true = 1/lambda_true, sigma2_true = 0, n = n, p = p, p0 = p0, cor = cor, data = F)
  data <- simdata$data
  
  #variables ('genes')
  X = as.matrix(data[,-1])
  
  #response, linear predictor
  Y = data$y
  
  #logistic
  #generate TRUE response
  Yl <- rbinom(n,size=N, exp(Y)/(1+exp(Y)))
  Yl <- cbind(Yl,N-Yl)
  if(N==1) {Yl <- 2*Yl;
  #w <- rep(1/sqrt(2),n)
  w <- rep(1/2,n)
  }
  else w <- rep(1,n)
  pmt <- proc.time()
  XXT = X %*% t(X)
  
  #### MOORE-PENROSE INVERSE FROM SVD
  pmt <- proc.time()
  svdXXT <- svd(XXT)  #rank of matrix is n-1 due to scaling...
  svdd <- svdXXT$d
  #reci <- 1/svdd[1:n]
  reci <- c(1/svdd[1:n-1],0)
  XXTi <- svdXXT$v %*% diag(reci) %*% t(svdXXT$u)
  ctsvd <- (proc.time()-pmt)[3]
  
  Xgam <- diag(n)
  PP = list(Xgam=list(XXTi),sp=-1)

  #MML
  pmt <- proc.time()
  pred <- try(gam(Yl ~ 0 + Xgam,family="binomial", weights=w, paraPen=PP,method="ML",scale=1,
                  optimizer=c("outer", "optim")))
  #pred <- try(gam(Yl ~ 0 + Xgam,family="binomial", paraPen=PP,method="ML",scale=1))
  if(class(pred) =="try-error")  lambdaest <- NA else lambdaest <- pred$sp
  print(lambdaest)
  ctmml <- (proc.time()-pmt)[3] + ctsvd
  
  
  #GCV
  pmt <- proc.time()
  pred2 <- gam(Yl ~ 0 + Xgam,family="binomial",paraPen=PP, weigths = w, method="GCV.Cp",scale=1)
  lambdaest2 <- pred2$sp
  ctgcv <- (proc.time()-pmt)[3] + ctsvd
  #lambdaest2
  
  #glmnet CV
  pmt <- proc.time()
  cvg <- cv.glmnet(X,Yl,standardize=FALSE, family="binomial", nfolds=10,alpha=0,intercept=FALSE)
  lambdaglmnet <- cvg$lambda.min * n
  ctglm <-(proc.time()-pmt)[3]
  lambdas <- c(lambda_true,lambdaest,lambdaest2,lambdaglmnet)
  cts <- c(ctmml,ctgcv,ctglm)
  
  #penalized
  if(penalized){
  pmt <- proc.time()
  pen <- try(optL2(Yl ~ 0 + X, model="logistic",fold=10,maxlambda2 = 10^10))
  if(class(pen)=="try-error") lambdapen <- Inf else lambdapen <- pen$lambda
  ctpen <-(proc.time()-pmt)[3]
  lambdas <- c(lambdas,lambdapen)
  cts <- c(cts,ctpen)
  }
  #penmod <- penalized(Yl ~ 0 + X, model="poisson",lambda2=lambdapen)
  #coef <- penmod@penalized
  #penmodMML <- penalized(Yl ~ 0 + X, model="poisson",lambda2=lambdaest)
  # coefMML <- penmodMML@penalized
  print(lambdas)
  if(penalized) names(lambdas)<- c("true","MML_mgcv","GCV_mgcv", "CV_glmnet", "CV_pen") else names(lambdas)<- c("true","MML_mgcv","GCV_mgcv", "CV_glmnet")
  df <- data.frame("lambda"=lambdas, "computing times" = c(NA,cts))
  return(df)
}





