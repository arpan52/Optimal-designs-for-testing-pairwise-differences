# for pairs (1,2),(2,3)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(nloptr)
library(mvtnorm)

IUT_simple_K3 <- function(x,mu,n,sigma,q){
  
  # Calculating mean and covariance matrix of Z_ij 's
  
  t<-(1-2*x)/2
  mu12 <- (n*(mu[1]-mu[2]))/(sigma*sqrt((1/x)+(1/t)))
  mu23 <- (n*(mu[2]-mu[3]))/(sigma*sqrt((1/t)+(1/x)))
  
 
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p23<-pnorm(q,mu23)-pnorm(-q,mu23)
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  rho1223 <- -sqrt((x*x)/((x+t)*(t+x)))
  s1223<- rbind(c(1,rho1223),c(rho1223,1))
  mu1223<- c(mu12,mu23)
  p1223<- pmvnorm(mean=mu1223, sigma=s1223, lower=lw1, upper=up1)      

  p<- -(1-p12-p23+p1223)         # Power of the test
  return(p)
}

# calculating the optimal design.

del<-0.55                        # delta value
mu<- c(0,del,2*del)              # mu value
N<- 180                          # Total number of subjects
n<- sqrt(N) 
al <- 0.05                    # significance level
q <- qnorm(1-(al))              # Critical value
sigma<-1

# Lower and upper bounds
lb <- 0
ub <- 0.5

#initial values
x0 <- 0.3

local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_simple_K3,
                mu = mu,
                n = n,
                sigma = sigma,
                q = q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)                 
x2<- 1-2*res$solution
x2