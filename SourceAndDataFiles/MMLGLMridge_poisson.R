##################### source function for poisson ridge regression ############
wrapper_poisson <- function(lambda_true,n,p,p0,cor,data=F, penalized=TRUE){
  simdata <- Simulatedata(tau2_true = 1/lambda_true, sigma2_true = 0, n = n, p = 1000, p0 = 0, cor = 0, data = F)
  data <- simdata$data
  
  #variables ('genes')
  X = as.matrix(data[,-1])
  
  #response, linear predictor
  Y = data$y
  
  #Poisson
  #generate TRUE response
  Yp <- rpois(n, exp(Y))
  Xgam <- diag(n)
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
  pred <- gam(Yp ~ 0 + Xgam,family="poisson",paraPen=PP,method="ML")
  ctmml <- (proc.time()-pmt)[3] + ctsvd
  lambdaest <- pred$sp
  #lambdaest
  
  #GCV
  pmt <- proc.time()
  pred2 <- gam(Yp ~ 0 + Xgam,family="poisson",paraPen=PP,method="GCV.Cp")
  lambdaest2 <- pred2$sp
  ctgcv <- (proc.time()-pmt)[3] + ctsvd
  #lambdaest2
  
  #glmnet CV
  pmt <- proc.time()
  cvg <- cv.glmnet(X,Yp,standardize=FALSE, family="poisson", nfolds=10,alpha=0,intercept=FALSE)
  lambdaglmnet <- cvg$lambda.min * n
  ctglm <-(proc.time()-pmt)[3]
  
  #penalized
  if(penalized){
  pmt <- proc.time()
  pen <- try(optL2(Yp ~ 0 + X, model="poisson",fold=10,maxlambda2 = 10^10))
  if(class(pen)=="try-error") lambdapen <- Inf else lambdapen <- pen$lambda
  ctpen <-(proc.time()-pmt)[3]
  }
  lambdas <- c(lambda_true,lambdaest,lambdaest2,lambdaglmnet)
  cts <- c(ctmml,ctgcv,ctglm)
  if(penalized) {lambdas <- c(lambdas,lambdapen); cts <- c(cts,ctpen)}
  if(penalized) names(lambdas)<- c("true","MML_mgcv","GCV_mgcv", "CV_glmnet", "CV_pen") else names(lambdas)<- c("true","MML_mgcv","GCV_mgcv", "CV_glmnet")
  df <- data.frame("lambda"=lambdas, "computing times" = c(NA,cts))
  return(df)
}





