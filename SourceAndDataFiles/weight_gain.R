### COMMENTS BY KRISTOFFER HELLTON : 
# Analysis of gene expression in Adipose tissue
# Outline
# 1. Select top 1000 genes correlated with outcome for analysis
# 2. Find ridge estimate by LOOCV, inspect the value 
# 3. Estimate the (in-sample) sigma^2 by ridge tuned by LOOCV
# 4. Estimate out-of-sample the LOOCV estimatse and fridge tuning parameters (with the ridge in-sample LOOCV plug-in) 
# 5. Visually inspect each MSE curve to check that the global optimum has been found
# 6. Plot result and calculate MSE 

y.raw <- weightData[colnames(exp),'WeightChangePct']*100

n <- length(y.raw)

#prior gene selection, based on marginal correlation
index <- sort(abs(apply(exp,1,function(x)cor(x,y.raw))),index.return = TRUE,decreasing = TRUE)$ix 

#Add intercept as first variable
dat <- t(exp[index[1:1000],])
X.raw0 <- cbind(rep(1,n),dat)

tuning.cv <- function(lambda,k=n-1){
  H <- svd$u[,1:k]  %*%  diag(svd$d[1:k]^2/(svd$d[1:k]^2 +lambda)) %*% t(svd$u[,1:k])
  e <- (diag(k+1) - H) %*% y
  mean((e/(1-diag(H)))^2)
}

tuning.fridge.ridge <- function(lambda,x0,cv.tuning, sigma2_hat){
  return( max(c((lambda*t(x0)%*%svd$v %*% diag((svd$d)/(svd$d^2 + lambda)/(svd$d^2+cv.tuning)) %*% t(svd$u) %*% y)^2-
                  sigma2_hat*lambda^2*t(x0)%*%svd$v  %*% diag(svd$d^2/(svd$d^2 + cv.tuning)^2/(svd$d^2 + lambda)^2)%*% t(svd$v)%*%x0,0))
          +  sigma2_hat*t(x0)%*%svd$v%*% diag(svd$d^2/(svd$d^2 + lambda)^2) %*% t(svd$v)%*%x0)
}

tuning.fridge.ridge.sim <- function(lambda,x0,cv.tuning, sigma2_hat){
  return((lambda*t(x0)%*%svd$v %*% diag(1/(svd$d^2 + lambda)*(svd$d)/(svd$d^2+cv.tuning)) %*% t(svd$u) %*% y)^2
         +  sigma2_hat*t(x0)%*%svd$v%*% diag(svd$d^2/(svd$d^2 + lambda)^2) %*% t(svd$v)%*%x0)
}

## Estimating the variance
svd <- svd(X.raw0)

p <- dim(X.raw0)[2]
n <- dim(X.raw0)[1]
#mean(y.raw)

#Find the optimal tuning value by LOO cross-validation
y <- y.raw
#curve(sapply(x,tuning.cv,k=n-1),0,10^2,10000)

opt.lambda <- optim(par = 10,tuning.cv,k=n-1,lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e3))$par
# opt.lambda

sigma2_hat <- sum((y - svd$u%*%diag(svd$d^2/(svd$d^2 + opt.lambda))%*%t(svd$u)%*%y)^2)/(n-sum(diag(svd$u%*%diag(svd$d^2/(svd$d^2 + opt.lambda))%*%t(svd$u))))
# sigma2_hat


#Find the optimal tuning value by MML estimation  
XXT = X.raw0 %*% t(X.raw0)
optimal <- optim(par = c(log(0.01),log(10)), fn = MargLikelihood, XXT = XXT, Y = y.raw)
sigma2mml <- exp(optimal$par[2])
lambdamml <- exp(optimal$par[2])/exp(optimal$par[1])
# sigma2mml
# lambdamml

#Comparing predictive perfromance of ridge+lambdaCV; ridge + lambdaMML; fridge + lambdaCV; fridge + lambdaMML
y0.fridge <- c(); y0mml.fridge <- c();y0cv <- c();  y0mml <- c();
errcv.fridge <- c(); errmml.fridge <- c(); errcv <- c(); errmml <- c();
opt.cv <- c(); opt.fridge <- c(); optm.fridge <- c(); optm.cv <- c();
y0.vec <- c()

print("Started LOOCV for predictive performance evaluation")
if(CVgenesel) print("Prior gene selection is included in outer-LOOCV loop")
for(i in 1:n){
  #i<-1
  print(paste("Left out sample:",i))
 
  if(CVgenesel){
  index <- sort(abs(apply(exp[,-i],1,function(x)cor(x,y.raw[-i]))),index.return = TRUE,decreasing = TRUE)$ix 
  #Add intercept as first variable
  dat <- t(exp[index[1:ngene],])
  X.raw <- cbind(rep(1,n),dat)
  } else X.raw <- X.raw0
  
  X <- X.raw[-i,] 
  y <- y.raw[-i]
  svd <- svd(X)
  
  x0 <- X.raw[i,]
  y0 <- y.raw[i]
  ntrain <- n-1
  
  opt.cv[i] <- optim(par = 1,tuning.cv,k=ntrain-1, lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e4))$par
  
  XXT = X %*% t(X)
  optimal <- optim(par = c(log(0.01),log(10)), fn = MargLikelihood, XXT = XXT, Y = y)
  optm.cv[i] <- exp(optimal$par[2])/exp(optimal$par[1])
  
  #Start the optimization at three different values to explore the whole range
  start <- 10^(0:10)
  
  value <- c(); param <- c();
    for(l in 1:length(start)){
    optimfridge <- optim(par = start[l],tuning.fridge.ridge.sim,x0=x0,
                         cv.tuning=opt.cv[i],lower=0,upper = Inf, sigma2_hat=sigma2_hat,
                         method = "L-BFGS-B",control = list(factr=1e-4))
    value[l] <- optimfridge$value   
    param[l] <- optimfridge$par  
    }
  
  opt.fridge[i]  <- param[which.min(value)] 
#  print(c(opt.cv[i],opt.fridge[i]))
  
  value <- c(); param <- c();
  for(l in 1:length(start)){
    optimfridge <- optim(par = start[l],tuning.fridge.ridge.sim,x0=x0,
                         cv.tuning=optm.cv[i],lower=0,upper = Inf, sigma2_hat=sigma2mml,
                         method = "L-BFGS-B",control = list(factr=1e-4))
    value[l] <- optimfridge$value   
    param[l] <- optimfridge$par  
  }
  optm.fridge[i]  <- param[which.min(value)] 

#fridge prediction, CV
  y0.fridge[i] <- t(x0)%*%svd$v%*%diag(svd$d/(svd$d^2  + opt.fridge[i]))%*%t(svd$u)%*% y
  
#fridge prediction, MML
  y0mml.fridge[i] <- t(x0)%*%svd$v%*%diag(svd$d/(svd$d^2  + optm.fridge[i]))%*%t(svd$u)%*% y
  
#ridge prediction, CV
  y0cv[i] <- t(x0)%*%svd$v%*%diag(svd$d/(svd$d^2  + opt.lambda))%*%t(svd$u)%*% y
  
#ridge prediction, MML
  y0mml[i] <- t(x0)%*%svd$v%*%diag(svd$d/(svd$d^2  + lambdamml))%*%t(svd$u)%*% y

#quadratic errors
  errmml.fridge[i] <- (y0mml.fridge[i]-y0)^2
  errcv.fridge[i] <- (y0.fridge[i]-y0)^2
  errcv[i] <- (y0cv[i]-y0)^2
  errmml[i] <- (y0mml[i]-y0)^2
  y0.vec[i] <- y0
}




