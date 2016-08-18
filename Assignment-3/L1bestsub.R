
L1bestsub=function(y1,x1,k)
{
library(quantreg)
y1=y1[,2]
df1<-data.frame(y1,x1)
regfit = rq(y1~.,data=df1)
max.num.pred=9; 
s=1;
I=1
absolute_deviation=rep(0,10);
k1_BIC=rep(0,10);
k2_BIC=rep(0,10);
penalty=ρτ=0
position=NULL;
posit=matrix(NA, nrow = 10, ncol = 10)

BIC=penalty=0 # stores the smallest error for each number of predictors
N=length(y1);
residual=matrix(NA, nrow = 480, ncol = 12)
min_residual=rep(0,10)

for(s in 1:max.num.pred) 
{  
  BICsmall=1000;     #intial big BIC value to easily find samller BIC values
  regfits=rep(0,11)  # initial predictor matrix
  pick=NULL;    
  pick=combn(regfit$coefficients[2:11],s) # all combinations of 's' prdectors
  count=1;
  y=y1;
  while(count!=dim(pick)[2])   # Loop to find minimum BIC value for all combinations of 's' predictors
  {   
    loc=regfit$coefficients%in%pick[,count] #stores location of current predictors used
    x=x1[,loc[2:10]]
    residual[count,s]=sum(rq(y1~x)$residuals)
    regfits[loc]=regfit$coefficients[loc]   #current predictor under testing
    regfits=na.omit(regfits)   
    u=0;  # loss function
    #find sum of loss function over all the data points
     count=count+1;
  }
  min=which.min(abs(residual[,s]))
  absolute_deviation[s]=abs(min(residual[,s],na.rm=TRUE))
  loc=regfit$coefficients%in%pick[,min]
  u= abs(min(residual[,s],na.rm=TRUE)) 
  if(k==1) {
    penalty=s*log(N)/(2*N);    #penalty term for loss function
     τ=0.5;#for median regression
    ρτ=abs(u*(2*τ-2*I))  
    k1_BIC[s]=(log((ρτ))+penalty); #calculate BIC
    }
  
  if (k==2)
  { penalty=s*log(N);
    k2_BIC[s]=abs((ρτ))+penalty; #calculate BIC
    τ=0.5;#for median regression
    ρτ=u*(2*τ-2*I) 
    }
position=list(position,t(as.matrix(which(loc))))
}
 position=unlist(position)
           
     #storing the position of the smallest value
    
  #k1_BIC=k2_BIC=penalty=ρτ=0
  if(k==1)
    {plot(k1_BIC[1:9],ylab="BIC value",xlab="number of predictors",main = "BIC values vs predictors",type='l')}
 else
  { plot(k2_BIC[1:9],ylab="BIC value",xlab="number of predictors",main = "BIC values vs predictors",type='l')}
 
}