library(NormalLaplace)
library(leaps)
M=10;
N=c(25,50,100);
intercept=0.5;
beta=c(0.3,0.25,0.2,0.15,0.1,0,0,0,0,0)



x1=matrix(rnorm(N[1]*M,mean=0,sd=1), N[1], M) 
y1=matrix(0, nrow = 25, ncol = 3)
y2=matrix(0, nrow = 25, ncol = 3)

df1=data.frame(x1,y1)
df2=data.frame(x1,y1)
for(j in 1:3)
{
  x1=matrix(rnorm(N[1]*M,mean=0,sd=1), N[1], M) 
  for(i in 1:25)
  {
beta=c(0.3,0.25,0.2,0.15,0.1,0,0,0,0,0)
error.rnorm=rnorm(1)
error.lap=rnl(1,0,1)
y1[i,j]=intercept+beta%*%x1[i,]+error.rnorm
y2[i,j]=intercept+beta%*%x1[i,]+error.lap

}
  
  regfit_25.rnorm=regsubsets(y1[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df1,nvmax=10);
  regfit_25.laplace=regsubsets(y2[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df2,nvmax=10);
}


y1=as.data.frame(y1)
names(y1)=c('response1','response2','response3')
data=data.frame(y1[,3],x1)
#fix(data)


###50 rows
x1=matrix(rnorm(N[2]*M,mean=0,sd=1), N[2], M) 
y1=matrix(0, nrow = 50, ncol = 3)
y2=matrix(0, nrow = 50, ncol = 3)

for(j in 1:3)
{
  x1=matrix(rnorm(N[2]*M,mean=0,sd=1), N[2], M) 
  for(i in 1:50)
  {
    error.rnorm=rnorm(1)
    error.lap=rnl(1,0,1)
    y1[i,j]=intercept+beta%*%x1[i,]+error.rnorm
    y2[i,j]=intercept+beta%*%x1[i,]+error.lap
         }
  
  regfit_50.rnorm=regsubsets(y2[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df1,nvmax=10);
  regfit_50.laplace=regsubsets(y2[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df2,nvmax=10);
}





###100 rows
x1=matrix(rnorm(N[3]*M,mean=0,sd=1), N[3], M) 
y1=matrix(0, nrow = 100, ncol = 3)
y2=matrix(0, nrow = 100, ncol = 3)

for(j in 1:3)
{
  x1=matrix(rnorm(N[3]*M,mean=0,sd=1), N[3], M) 
  for(i in 1:100)
  {
    error.rnorm=rnorm(1)
    error.lap=rnl(1,0,1)
    y1[i,j]=intercept+beta%*%x1[i,]+error.rnorm
    y2[i,j]=intercept+beta%*%x1[i,]+error.lap
    
  }
  
  regfit_100.rnorm=regsubsets(y2[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df1,nvmax=10);
  regfit_100.laplace=regsubsets(y2[,j]~x1[,1]+x1[,2]+x1[,3]+x1[,4]+x1[,5]+x1[,6]+x1[,7]+x1[,8]+x1[,9]+x1[,10],df2,nvmax=10);
}


















