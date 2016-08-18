theta <- seq(-pi, pi,0.01)
x <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96, 2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52, 2.50)


log.lik <- function(x, theta) {-20*log(2*pi) + sum(log(1-cos(x-theta)))}
plot(theta, sapply(theta, log.lik, x=x), typ='l',main="Log-Likelihood Function",ylab=" Function ") 


log.lik2 <- function(x){ return (-20*log(2*pi) + sum(log(1-cos(x-theta))))}



newton <- function(f, f.derivative, x0, eps, maxter){
  x1 <- x0 - f(x0)/f.derivative(x0)
  ITERATION <- 1
  while(abs(x1 - x0) > eps & ITERATION < maxter){
    x0 <- x1
    x1 <- x0 - f(x0)/f.derivative(x0)
    cat("[ITERATION]", ITERATION, "  ", x1, fill=T)
    ITERATION <- ITERATION + 1
  }
  
  return (x1)
}

# First derivative of log-likelihood function
first_derivative <- function(theta){return (-sum(sin(x-theta)/(1-cos(x-theta))))}

# Second derivative
second_derivative<- function(theta){return (-sum(1/(1-cos(x-theta))))}

s=seq(-pi,pi,length.out=200)
newton(first_derivative, second_derivative, -0.0584, 1e-6, 100)
newton(first_derivative, second_derivative, -2.7, 1e-6, 100)
newton(first_derivative, second_derivative, 2.7, 1e-6, 100)

thetas=NULL
for(i in 1:200)
{
  thetas[i]=newton(first_derivative, second_derivative, s[i], 1e-3, 100)
  
}


tol=1e-1
r=rep(0,100)
for (i in 1:100)
{
  r[i]= newton(first_derivative, second_derivative, (s[11]+i*tol), 1e-6, 100)
  value=s[11]+i*tol
  if (theta[11]!=s[11])
  {
    break;
    print(i)
  }
}


