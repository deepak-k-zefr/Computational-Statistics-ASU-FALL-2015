
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment<-gl(3,3)
print (df <- data.frame(treatment,outcome,counts))

beta.1<-rep(0,5)


treatment<-as.numeric(levels(treatment))[treatment]
outcome<-as.numeric(levels(outcome))[outcome]

treatment2<-rep(0,length(counts))
treatment3<-rep(0,length(counts))

outcome2<-rep(0,length(counts))
outcome3<-rep(0,length(counts))
y=counts

for(i in 1:length(counts)){
  
  if(treatment[i]==2){
    treatment2[i]<-1}
  else if(treatment[i]==3){
    treatment3[i]<-1}
  if(outcome[i]==2){
    outcome2[i]<-1}
  else if(outcome[i]==3){
    outcome3[i]<-1}
}

df1 <- data.frame(outcome2,outcome3,treatment2,treatment3,counts)
df <- data.frame(counts,outcome2,outcome3,treatment2,treatment3)
y<-as.vector(df[,1])
X<-cbind(rep(1,nrow(X)),X)
X<-cbind(rep(1,length(counts)),df1[,1],df1[,2],df1[,3],df1[,4])
colnames(X) <-c("int","outcome2","outcome3","treatment2","treatment3")
system.time(out<-main_f(X,counts,beta.1))

beta.Est<-out$beta.MLE
beta.SE<-sqrt(diag(out$beta.cov))
beta.z<-beta.Est/beta.SE
beta.pval<-2*pnorm(-abs(beta.z))
beta.2.5<-beta.Est-1.96*beta.SE
beta.97.5<-beta.Est+1.96*beta.SE
beta.coef<-data.frame(beta.Est,beta.SE,beta.z,beta.pval,beta.2.5,beta.97.5)


library(Bhat)
lpoisson<-function(theta){
  
  #Poisson  log likelihood function
  #input:counts, sample sizes 
  #output: log-likelihood
  
  
  l<-t(y)%*%(X%*%theta) - sum(exp(X%*%theta))
  
  return(l)
  
}

d1lpoisson <- function(theta){
  d1l<-NULL
  for(i in 1:ncol(X)){
    d1l[i]=1*((sum(y*X[,i])-sum(exp(X%*%theta)*X[,i])))
  }
  return(d1l)
}


d2lpoisson <- function(theta){
  d2l<-NULL
  for(i in 1:ncol(X)){
    d2l[i]=-sum((X[,i]^2)*exp(X%*%theta))
  }
  h<-diag(d2l)
  return(h)
}

library(maxLik)
out<-maxNR(lpoisson,d1lpoisson,d2lpoisson,betas)
out<-main_f(X,counts,beta.1)
