################################################
#MONTE_CARLO SIMULATIONS ON MULTIVARIATE T DATA
################################################

###
# Run the function file  attached along with this R file before running this code
###

setwd("D:/R")
library(mvtnorm)
require(flare)
require(glmnet)
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(lars)
library(ggplot2)
library(quantreg)


###-----generate data-----

genData <- function(choice)
{
  switch(choice,
         
         ##################################################
         #Scenario 1: FOUR REAL PREDICTORS, 16 ZERO Predictors
         ##################################################
         a={
           
           X <- matrix(rt(n*p,5),nrow=n,ncol=p)
           y <- rnorm(n,X%*%beta,2)
           
         },
         ##################################################
         #Scenario 2: ALL PREDICTORS ARE EQUAL
         ##################################################
         b={
           X <- matrix(rt(n*p,5),nrow=n,ncol=p)
           y <- rnorm(n,X%*%beta,2)
           
           
         },
         ####################################################
         #Scenario 3: ALL PREDICTORS ARE STRONGLY CORRELATED
         ####################################################
         c={
           
           x=c(1,rho^seq(1:(length(beta)-1)))
           x=toeplitz(x)
           X=matrix(rmvt(n*p,x,5),nrow=n ,ncol=p)
           y=rnorm(n,X%*%beta,2.2)
           
         },
         ##################################################
         #Scenario 4:HIGH DIMENSIONAL SETTING
         ##################################################
         d={
           X=matrix(rt(n*p,5),nrow=n,ncol=p )
           y=rnorm(n,X%*%beta,2)
           
         },
         ####################################################
         #Scenario 5: HIGH DIMENSIONAL SETTING AND STRONGLY CORRELATED
         ####################################################
         e={
           x=c(1,rho^seq(1:(length(beta)-1)))
           x=toeplitz(x)                         #create variance matrix
           X=matrix(rmvt(n*p,x,5),nrow=n,ncol=p ) # degrees of freedom=5
           y=rnorm(n,X%*%beta,2)
           
           
         })
  
  return(list(X=X,y=y))
}

####################################################
### 1).ELASTIC NET 
####################################################
elasticnet<-function(Data)
{
  
  k=10   # 10 fold cross validation
  fold <- sample(1:k, size=length(Data$y), replace=TRUE)
  lambda_values <- NULL   # list of lambda values
  cv.error <- NULL     #  error matrix
  alpha_values=seq(0:10)/10 # list of alpha values
  #10-fold Cross validation for each alpha = 0, 0.1, ... , 0.9, 1.0
  for (i in 1:10)     # i is alpha level
  {
    cvfit <- cv.glmnet(Data$X, Data$y, foldid=fold, alpha=i/10, type.measure="mse",family="gaussian")
    lambda_values[i] <- cvfit$lambda.min     # best lambda for alpha=i/10
    min.index <- which(cvfit$lambda==cvfit$lambda.min) # store location of best lambda
    cv.error[i] <- cvfit$cvm[min.index]   # store errors
  }
  min.index<-which.min(cv.error)     # locate smallest error
  selected_model<-glmnet(Data$X, Data$y, alpha=alpha_values[min.index], lambda=lambda_values[min.index])     # refit model using optimal alpha and lambda 
  
  return(as.numeric(coef(selected_model)))
}

####################################################
### 2). FLARE LQ LASSO 
####################################################

flare_lq_lasso=function(Data)
{
  result_1 = slim(Data$X,Data$y,method="lq",nlambda=10,q=1.5,lambda.min.ratio=0.35)
  Y_original= matrix(rep(Data$X%*%beta, each=10), ncol=10)
  residuals=abs((Data$X%*%result_1$beta)-Y_original)
  mse_lq_lasso<-(apply(residuals, 2,mean))^2
  min.index=which.min(mse_lq_lasso)  
  result_1$lambda[min.index] #Best performing lambda
  result_1$beta[,min.index]
  return(result_1$beta[,min.index])#Best performing lambda
}
####################################################
### 3). FLARE LAD LASSO 
####################################################
flare_lad_lasso=function(Data)
{
  result_2 = slim(Data$X,Data$y,method="lq",nlambda=10,q=1,lambda.min.ratio=0.3)
  Y_original= matrix(rep(Data$X%*%beta, each=10), ncol=10)
  residuals=abs((Data$X%*%result_2$beta)-Y_original)
  mse_lad_lasso<-(apply(residuals, 2,mean))^2
  min.index=which.min(mse_lad_lasso)  
  result_2$lambda[min.index] #Best performing lambda
  result_2$beta[,min.index]
  return(result_2$beta[,min.index])#Best performing lambda
}
####################################################
### 4). FLARE SQUARE LASSO 
####################################################
flare_sqr_lasso=function(Data)
{
  result_2 = slim(Data$X,Data$y,method="lq",nlambda=10,q=2)
  Y_original= matrix(rep(Data$X%*%beta, each=10), ncol=10)
  residuals=abs((Data$X%*%result_2$beta)-Y_original)
  mse_sqr_lasso<-(apply(residuals, 2,mean))^2
  min.index=which.min(mse_sqr_lasso)  
  result_2$lambda[min.index] #Best performing lambda
  result_2$beta[,min.index]
  return(result_2$beta[,min.index])#Best performing lambda
}
###############################
### 5). FLARE DANTZIG LASSO 
###############################
flare_dantzig_lasso=function(Data)
{
  result_3 = slim(Data$X,Data$y,method="dantzig",nlambda=10)
  #beta1=rbind(result_3$intercept,result_3$beta) # combining the intercepts and the beta
  Y_original= matrix(rep(Data$X%*%beta, each=10), ncol=10)
  residuals=abs((Data$X%*%result_3$beta)-Y_original)
  mse_dantzig<<-(apply(residuals, 2,mean))^2
  min.index=which.min(mse_dantzig)  
  result_3$lambda[min.index] #Best performing lambda
  result_3$beta[,min.index]
  return(result_3$beta[,min.index])#Best performing lambda
  
  B <- apply(results,2:3,mean)-beta
  V <- apply(results,2:3,var)
  MSE <- B^2+V
}
########################################
### 7). ADAPTIVE LAD LASSO LAMBDA=1
########################################
#Adaptive LAD Lasso for four cases of lambda
aladlasso_lambda1 <- function(Data){
  Lambda=sqrt(1.5*n*log(p))#original penalty
  weight<- coef(rq(Data$y~Data$X,method="lasso"))#Weight=LAD coefficents
  weight<-weight[-1]  
  weighted_lambda=Lambda/(abs(weight))
  weighted_lambda<-c(0,weighted_lambda)
  tempcfQR<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = weighted_lambda))
  return(tempcfQR)
}

########################################
### 8). ADAPTIVE LAD LASSO LAMBDA=2
########################################
aladlasso_lambda2 <- function(Data){
  Lambda=sqrt(2*n*log(p))
  weight<- coef(rq(Data$y~Data$X,method="lasso"))#Weight=LAD coefficents
  weight<-weight[-1]  
  weighted_lambda=Lambda/(abs(weight))
  weighted_lambda<-c(0,weighted_lambda)
  tempcfQR<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = weighted_lambda))
  return(tempcfQR)
}


########################################
### 9). ADAPTIVE LAD LASSO LAMBDA=3
########################################
aladlasso_lambda3 <- function(Data){
  Lambda=sqrt(4*n*log(p))
  weight<- coef(rq(Data$y~Data$X,method="lasso"))#Weight=LAD coefficents
  weight<-weight[-1]  
  weighted_lambda=Lambda/(abs(weight))
  weighted_lambda<-c(0,weighted_lambda)
  tempcfQR<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = weighted_lambda))
  return(tempcfQR)
}

########################################
### 10). ADAPTIVE LAD LASSO LAMBDA=4
########################################
aladlasso_lambda4 <- function(Data){
  Lambda=sqrt(10*n*log(p))
  weight<- coef(rq(Data$y~Data$X,method="lasso"))#Weight=LAD coefficents
  weight<-weight[-1]  
  weighted_lambda=Lambda/(abs(weight))
  weighted_lambda<-c(0,weighted_lambda)
  tempcfQR<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = weighted_lambda))
  return(tempcfQR)
}
########################################
### 6). ADAPTIVE LAD LASSO LARS
########################################
lsa.linear<-function(Data){
  # adaptive lasso for linear reg, tuning parameter by bic 
  # calls software from Wang and Leng (2007, JASA).
  # since regsubsets can't handle na's
  as.matrix(Data$X)->Data$X
  lm(Data$y~Data$X)->out
  lsa(out)->out.lsa
  coeff<-out.lsa$beta.bic
  coeff2<-coeff[2:(p+1)]               # get rid of intercept
  pred<-Data$X%*%coeff2+coeff[1]
  st<-sum(coeff2 !=0)                                          # number nonzero
  mse<-sum((Data$y-pred)^2)/(n-st-1)
  if(st>0) x.ind<-as.vector(which(coeff2 !=0)) else x.ind<-0
  return(coeff2)
}

cat("\014")
choice=readline("Choose from one of the 5 scenarios
                a.Uncorrelated              b.All Predictors are equal               c.Strongly Correlated Variables              d.High Dimensionsal Setting
                e.Correlated Predictors in High Dimensions     :")


switch(choice,
       a={
         n=100
         N <- 50
         beta <- c(2,-2,1,-1,rep(0,16))
         p <- length(beta)
         results <<- array(NA,dim=c(N,p,10),
                           dimnames=list(1:N,1:p,c("ELASTIC-NET","FLARE_LQ_LASSO","FLARE_LAD_LASSO","FLARE_SQR_LASSO", "FLARE_DANTZIG_LASSO","ADAPTIVE_LASSO", 
                                                   "ADAPTIVE_LAD-LASSO_LAMBDA(1)","ADAPTIVE_LAD-LASSO_LAMBDA(2)","ADAPTIVE_LAD-LASSO_LAMBDA(3)","ADAPTIVE_LAD-LASSO_LAMBDA(4)")))
         
       },b={
         n=100
         N=150
         beta=rep (0.2,20) #ALL PREDICTORS ARE EQUAL
         p=length(beta)
         results <<- array(NA,dim=c(N,p,10),
                           dimnames=list(1:N,1:p,c("ELASTIC-NET","FLARE_LQ_LASSO","FLARE_LAD_LASSO","FLARE_SQR_LASSO", "FLARE_DANTZIG_LASSO","ADAPTIVE_LASSO", 
                                                   "ADAPTIVE_LAD-LASSO_LAMBDA(1)","ADAPTIVE_LAD-LASSO_LAMBDA(2)","ADAPTIVE_LAD-LASSO_LAMBDA(3)","ADAPTIVE_LAD-LASSO_LAMBDA(4)")))
         
       },c={
         n=100;N=150
         beta=c(2,-2,1,-1,0.5,0.2,-0.3,-0.15,rep(0,12))
         p=length(beta) 
         rho=0.8
         results <<- array(NA,dim=c(N,p,10),
                           dimnames=list(1:N,1:p,c("ELASTIC-NET","FLARE_LQ_LASSO","FLARE_LAD_LASSO","FLARE_SQR_LASSO", "FLARE_DANTZIG_LASSO","ADAPTIVE_LASSO", 
                                                   "ADAPTIVE_LAD-LASSO_LAMBDA(1)","ADAPTIVE_LAD-LASSO_LAMBDA(2)","ADAPTIVE_LAD-LASSO_LAMBDA(3)","ADAPTIVE_LAD-LASSO_LAMBDA(4)")))
         
       },d={
         n=100;N=150
         beta=c(2,1,0.5,rep(0,997)) #1000 predictors
         p=length(beta)
         results <<- array(NA,dim=c(N,p,10),
                           dimnames=list(1:N,1:p,c("ELASTIC-NET","FLARE_LQ_LASSO","FLARE_LAD_LASSO","FLARE_SQR_LASSO", "FLARE_DANTZIG_LASSO","ADAPTIVE_LASSO", 
                                                   "ADAPTIVE_LAD-LASSO_LAMBDA(1)","ADAPTIVE_LAD-LASSO_LAMBDA(2)","ADAPTIVE_LAD-LASSO_LAMBDA(3)","ADAPTIVE_LAD-LASSO_LAMBDA(4)")))
       },e={
         rho=0.8
         n=100
         N=150
         beta=c(2,-2,1,-1,0.5,0.2,-0.3,-0.15,rep(0,212))
         p=length(beta)
         results <<- array(NA,dim=c(N,p,10),
                           dimnames=list(1:N,1:p,c("ELASTIC-NET","FLARE_LQ_LASSO","FLARE_LAD_LASSO","FLARE_SQR_LASSO", "FLARE_DANTZIG_LASSO","ADAPTIVE_LASSO", 
                                                   "ADAPTIVE_LAD-LASSO_LAMBDA(1)","ADAPTIVE_LAD-LASSO_LAMBDA(2)","ADAPTIVE_LAD-LASSO_LAMBDA(3)","ADAPTIVE_LAD-LASSO_LAMBDA(4)")))
       })

for (i in 1:N)
{
  
  Data <- genData(choice)
  results[i,,1]<-elasticnet(Data)[-1]
  results[i,,2]<-flare_lq_lasso(Data)
  results[i,,3]<-flare_lad_lasso(Data)
  results[i,,4]<-flare_sqr_lasso(Data)
  results[i,,5]<-flare_dantzig_lasso(Data)
  results[i,,6]<-lsa.linear(Data)
  results[i,,7]<-aladlasso_lambda1(Data)[-1]
  results[i,,8]<-aladlasso_lambda2(Data)[-1]
  results[i,,9]<-aladlasso_lambda3(Data)[-1]
  results[i,,10]<-aladlasso_lambda4(Data)[-1]
  
  #displayProgressBar(i,N)
}



B <- apply(results,2:3,mean)-beta
V <- apply(results,2:3,var)
MSE <- B^2+V
MSE_SUM=apply(MSE,2,sum)
results=results[,,-6]

ylim <- range(results)
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
  
  
}


setwd("D:/R")


pheatmap(B[1:20,],  colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(100),cluster_rows=F, cluster_cols=F, main="Bias  ")
pheatmap(V[1:20,],  colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(100),cluster_rows=F, cluster_cols=F, main="Variance   ")

png("MC2.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B[1:20,])
dev.off()


png("MC3.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(V[1:20,])
dev.off()

png("MC4.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,1],col="gray",ylim=ylim,main="Elastic net Lasso")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC5.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,2],col="bisque",ylim=ylim,main="Flare-Lq Lasso")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC6.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,3],col="cadetblue",ylim=ylim,main="Flare-lad Lasso ")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC7.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,4],col="beige",ylim=ylim,main="Flare-sqr Lasso")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC8.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,5],col="azure",ylim=ylim,main="Flare-Dantzig Lasso ")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC9.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,6],col="darkgoldenrod",ylim=ylim,main="Adaptive Lasso ")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC10.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,7],col="darkgoldenrod2",ylim=ylim,main="Ad. Lad Lasso-1")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC11.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,8],col="burlywood",ylim=ylim,main="Ad Lad Lasso-2 ")
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC12.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,9],col="burlywood",ylim=ylim,main="Ad Lad Lasso-3" )
abline(h=c(2,1,0.5), lty=2)
dev.off()

png("MC13.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,1:20,10],col="burlywood",ylim=ylim,main="Ad Lad Lasso-4 ")
abline(h=c(2,1,0.5), lty=2)
dev.off()




# ----- END plot function ----- #
