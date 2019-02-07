plotHyperpar <- function(output, includeBayes=FALSE, mixed = FALSE, onerow=TRUE, methinclude=NULL){
  
#output <- simulate10; methtau =  c("MML","MoM"); methsig=c("MML", "MoM", "PCR", "Basic"); methlambda=c("MML","MoM","GCV", "CV"); methh2 = c("MML", "MoM", "HiLMM")
nsim <- length(output)
if(mixed){
  if(is.null(methinclude)) methinclude <- c("MML", "REML", "CV")
  methtau <- intersect(methinclude,c("MML","REML"))
  methsig <- intersect(methinclude, c("MML","REML")) 
  methlambda <- intersect(methinclude,c("MML","REML", "CV"))
  methh2 <- intersect(methinclude,c("MML","REML"))
  } else {
    if(is.null(methinclude)) methinclude <- c("MML", "MoM", "PCR", "Basic","GCV", "CV","HiLMM")
    methtau <- intersect(methinclude,c("MML","MoM"))
    methsig <- intersect(methinclude, c("MML", "MoM", "PCR", "Basic")) 
    methlambda <- intersect(methinclude,c("MML","MoM", "GCV", "CV"))
    methh2 <- intersect(methinclude,c("MML","MoM", "HiLMM"))
  }
if(includeBayes && mixed) print("Bayes not implemented for mixed model. Results will only be computed for other methods.")
if(includeBayes && !mixed  && !is.element("Bayes", methtau)) {
  methtau <- c(methtau,"Bayes");methsig <- c(methsig,"Bayes");methlambda <- c(methlambda,"Bayes")
  methh2 <- c(methh2,"Bayes")
}

df_tau<- data.frame(t(sapply(1:nsim, function(i)  return(output[[i]][,1]))))
colnames(df_tau) <- rownames(output[[1]])
whtrue <- match(c("true"),colnames(df_tau))
tau2_true <- df_tau[1,whtrue]
wh <- match(methtau,colnames(df_tau))
df_tau <- df_tau[,wh]
df_tau2 <- suppressMessages(melt(df_tau))

df_sigma<- data.frame(t(sapply(1:nsim, function(i)  return(output[[i]][,2]))))
colnames(df_sigma) <- rownames(output[[1]])
whtrue <- match(c("true"),colnames(df_sigma))
sigma2_true <- df_sigma[1,whtrue]
wh <- match(methsig,colnames(df_sigma))
df_sigma <- df_sigma[,wh]
df_sigma2 <- suppressMessages(melt(df_sigma))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
color.index_tau = ifelse(is.wholenumber(as.numeric(df_tau2$variable)/2), color <- "red", color<-"blue")
color.index_sigma = ifelse(is.wholenumber(as.numeric(df_sigma2$variable)/2), color <- "red", color<-"blue")


p1 = ggplot(data = df_tau2, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = color.index_tau), outlier.shape=NA) + 
  geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
  #geom_jitter(width = 0.1) + 
  geom_hline(yintercept = tau2_true, color = "red", size = 1, linetype = "dashed") + 
  labs(y = expression(hat(tau)^2), x = " ") +  theme_classic() + 
  scale_fill_brewer(palette = "BuGn")+
  theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black") ,
        axis.text.x=element_text(angle=90, hjust=1,size=10) , 
        plot.margin=unit(c(10,2,-3,2),"mm"), #t,r,b,l
        panel.background = element_rect(colour = "black", size=0.5))
p2 = ggplot(data = df_sigma2, aes(x = variable,  y = value)) +
  geom_boxplot(aes(fill = color.index_sigma), outlier.shape=NA) + 
  geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
  geom_hline(yintercept = sigma2_true, color = "red", size = 1, linetype = "dashed") + 
  labs(y = expression(hat(sigma)^2), x = " ") + theme_classic() + 
  scale_fill_brewer(palette = "BuGn")+
  theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1,size=10),
        plot.margin=unit(c(10,2,-3,2),"mm"), #t,r,b,l
        panel.background = element_rect(colour = "black", size=0.5))
#grid.arrange(p1, p2, ncol = 2)
#grid.arrange(p1, ncol = 1)


##########FUNCTIONS OF THE VARIANCES##############
df_lambda<- data.frame(t(sapply(1:nsim, function(i)  return(output[[i]][,3]))))
colnames(df_lambda) <- rownames(output[[1]])
whtrue <- match(c("true"),colnames(df_lambda))
lambda_true <- df_lambda[1,whtrue]
df_lambda[which(df_lambda > (lambda_true*20), arr.ind = T)] <- lambda_true * 20
wh <- match(methlambda,colnames(df_lambda))
df_lambda <- df_lambda[,wh]
df_lambda2 <- suppressMessages(melt(df_lambda))
method = df_lambda2[,1] ; lambda_est = df_lambda2[,2]
index = which(lambda_est == (lambda_true*20))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
color.index_lambda = ifelse(is.wholenumber(as.numeric(method)/2), color <- "red", color<-"blue")


p3 = ggplot(data = df_lambda2, aes(x = variable, y = sqrt(value))) + 
  geom_boxplot(aes(fill = color.index_lambda), outlier.shape=NA) + 
  geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
  geom_hline(yintercept = sqrt(lambda_true), col = "red", size = 1, linetype = "dashed") +
  labs(y = expression(sqrt(hat(lambda))), x = " ") + theme_classic() + 
  scale_fill_brewer(palette = "BuGn")+
  theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1,size=10),
        plot.margin=unit(c(10,2,-3,2),"mm"), #t,r,b,l
        panel.background = element_rect(colour = "black", size=0.5))

#Plotting heritability
df_herit<- data.frame(t(sapply(1:nsim, function(i)  return(output[[i]][,4]))))
df_herit[which(df_herit > 1, arr.ind = T)] <- 1 #truncate at 1
df_herit[which(df_herit < 0, arr.ind = T)] <- 0 #truncate at 0

colnames(df_herit) <- rownames(output[[1]])
whtrue <- match(c("true"),colnames(df_herit))
herit2_true <- df_herit[1,whtrue]
wh <- match(methh2,colnames(df_herit))
df_herit <- df_herit[,wh]
df_herit2 <- suppressMessages(melt(df_herit))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
color.index_herit = ifelse(is.wholenumber(as.numeric(df_herit2$variable)/2), color <- "red", color<-"blue")

p4 = ggplot(data = df_herit2, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = color.index_herit), outlier.shape=NA) + 
  geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
  geom_hline(yintercept = herit2_true, col = "red", size = 1, linetype = "dashed") +
  geom_hline(yintercept = c(1, 0), col = "gray", linetype = "dashed", size = 0.3)+
  labs(    y = expression(hat(h)^2), x = " ") + theme_classic() + 
  scale_fill_brewer(palette = "BuGn")+
  theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1,size=10),
        plot.margin=unit(c(10,2,-3,2),"mm"), #t,r,b,l
        panel.background = element_rect(colour = "black", size=0.5))

if(onerow) grid.arrange(p1,p2,p3,p4,ncol=4,nrow=1) else grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
return(list(tauplot=p1,sigmaplot=p2, lambdaplot=p3,h2plot=p4))
}

plotHyperparlam <- function(output, methinclude=NULL, plotlog10=T){
#output <- simulatep; methinclude =  NULL; plotlog=T
  nsim <- length(output)
    if(is.null(methinclude)) methinclude <- c("MML_mgcv","GCV_mgcv", "CV_glmnet", "CV_pen")
    methlambda <- intersect(methinclude,c("MML_mgcv","GCV_mgcv", "CV_glmnet", "CV_pen"))

  df_lambda<- data.frame(t(sapply(1:nsim, function(i)  return(output[[i]][,1]))))
  colnames(df_lambda) <- rownames(output[[1]])
  whtrue <- match(c("true"),colnames(df_lambda))
  lambda_true <- df_lambda[1,whtrue]
  #df_lambda[which(df_lambda > (lambda_true*20), arr.ind = T)] <- lambda_true * 20
  wh <- match(methlambda,colnames(df_lambda))
  df_lambda <- df_lambda[,wh]
  df_lambda2 <- suppressMessages(melt(df_lambda))
  method = df_lambda2[,1] ; lambda_est = df_lambda2[,2]
  index = which(lambda_est == (lambda_true*20))
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  color.index_lambda = ifelse(is.wholenumber(as.numeric(method)/2), color <- "red", color<-"blue")
  
  if(plotlog10){
    transfy <- log10(df_lambda2$value)
    axlab <- expression(log10(hat(lambda)))
    truetr <- log10(lambda_true)
  } else {  #sqrt scale
    transfy <- sqrt(df_lambda2$value)
    axlab <- expression(sqrt(hat(lambda)))
    truetr <- sqrt(lambda_true)
  }
  
  df_lambda2 <- data.frame(df_lambda2,transfy=transfy)
  
  p3 <- ggplot(data = df_lambda2, aes(x = variable, y = transfy)) + 
    geom_boxplot(aes(fill = color.index_lambda), outlier.shape=NA) + 
    geom_jitter(width = 0.1, color = rgb(0,0,0,.2)) + 
    geom_hline(yintercept = truetr, col = "red", size = 1, linetype = "dashed") +
    labs(y = axlab, x = " ") + theme_classic() + 
    scale_fill_brewer(palette = "BuGn")+
    theme(legend.position="none",axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
          axis.text.x=element_text(angle=90, hjust=1,size=10),
          plot.margin=unit(c(3,2,-3,2),"mm"), #t,r,b,l
          panel.background = element_rect(colour = "black", size=0.5))
  show(p3)
  return(lambdaplot=p3)
}
