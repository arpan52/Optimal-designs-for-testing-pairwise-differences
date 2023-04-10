# for pairs (1,2),(2,3),(3,4),(4,5)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(nloptr)
library(mvtnorm)

IUT_simple_K5 <- function(x,del,mu,n,sigma,q){


  # Calculating mean and covariance matrix of Z_ij 's
  
  t1<-x[1]
  t2<-x[2]
  t3<- 1-(2*t1+2*t2)
  
  mu12 <- (n*(mu[1]-mu[2]))/(sigma*sqrt((1/t1)+(1/t2)))
  mu23 <- (n*(mu[2]-mu[3]))/(sigma*sqrt((1/t2)+(1/t3)))
  mu34 <- (n*(mu[3]-mu[4]))/(sigma*sqrt((1/t3)+(1/t2)))
  mu45 <- (n*(mu[4]-mu[5]))/(sigma*sqrt((1/t2)+(1/t1)))

  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p23<- pnorm(q,mu23)-pnorm(-q,mu23)
  p34<- pnorm(q,mu34)-pnorm(-q,mu34)
  p45<- pnorm(q,mu45)-pnorm(-q,mu45)
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  rho1223 <- -sqrt((t1*t3)/((t1+t2)*(t2+t3)))
  s1223<- rbind(c(1,rho1223),c(rho1223,1))
  mu1223<- c(mu12,mu23)
  p1223<- pmvnorm(mean=mu1223, sigma=s1223, lower=lw1, upper=up1)

  rho2334 <- -sqrt((t2*t2)/((t3+t2)*(t2+t3)))
  s2334<- rbind(c(1,rho2334),c(rho2334,1))
  mu2334<- c(mu23,mu34)
  p2334<- pmvnorm(mean=mu2334, sigma=s2334, lower=lw1, upper=up1)
  
  rho1234 <- 0
  s1234<- rbind(c(1,rho1234),c(rho1234,1))
  mu1234<- c(mu12,mu34)
  p1234<-  pmvnorm(mean=mu1234, sigma=s1234, lower=lw1, upper=up1)
  
  rho1245 <- 0
  s1245<- rbind(c(1,rho1245),c(rho1245,1))
  mu1245<- c(mu12,mu45)
  p1245<-  pmvnorm(mean=mu1245, sigma=s1245, lower=lw1, upper=up1)
  
  rho2345 <- 0
  s2345<- rbind(c(1,rho2345),c(rho2345,1))
  mu2345<- c(mu23,mu45)
  p2345<-  pmvnorm(mean=mu2345, sigma=s2345, lower=lw1, upper=up1)
  
  rho3445 <- -sqrt((t3*t1)/((t3+t2)*(t2+t1)))
  s3445<- rbind(c(1,rho3445),c(rho3445,1))
  mu3445<- c(mu34,mu45)
  p3445<-  pmvnorm(mean=mu3445, sigma=s3445, lower=lw1, upper=up1)
  
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  s122334<- rbind(c(1,rho1223,rho1234),c(rho1223,1,rho2334),c(rho1234,rho2334,1))
  mu122334<- c(mu12,mu23,mu34)
  p122334<-  pmvnorm(mean=mu122334, sigma=s122334, lower=lw2, upper=up2)
  
  s122345<- rbind(c(1,rho1223,rho1245),c(rho1223,1,rho2345),c(rho1245,rho2345,1))
  mu122345<- c(mu12,mu23,mu45)
  p122345<-  pmvnorm(mean=mu122345, sigma=s122345, lower=lw2, upper=up2)
  
  s123445<- rbind(c(1,rho1234,rho1245),c(rho1234,1,rho3445),c(rho1245,rho3445,1))
  mu123445<- c(mu12,mu34,mu45)
  p123445<-  pmvnorm(mean=mu123445, sigma=s123445, lower=lw2, upper=up2)
  
  s233445<- rbind(c(1,rho2334,rho2345),c(rho2334,1,rho3445),c(rho2345,rho3445,1))
  mu233445<- c(mu23,mu34,mu45)
  p233445<-  pmvnorm(mean=mu233445, sigma=s233445, lower=lw2, upper=up2)
  
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q)
  
  s12233445<- rbind(c(1,rho1223,rho1234,rho1245),c(rho1223,1,rho2334,rho2345),c(rho1234,rho2334,1,rho3445),c(rho1245,rho2345,rho3445,1))
  mu12233445<- c(mu12,mu23,mu34,mu45)
  p12233445<-  pmvnorm(mean=mu12233445, sigma=s12233445, lower=lw3, upper=up3)
  
  
    # Power of the test  
  p<- -(1-p12-p23-p34-p45+p1223+p1234+p1245+p2334+p2345+p3445-p122334-p122345-p123445-p233445+p12233445)       
  return(p)
}

# calculating the optimal design.

del<-0.7                          # delta value
mu<- c(0,del,2*del,3*del,4*del)     # mu value
sigma<-1
N<- 180                            # Total number of subjects
n<-sqrt(N)
al<-0.05                         # significance level
q<-qnorm(1-al)                     # critical value


#initial values
x0 <-c(0.15,0.18)

# Lower and upper bounds
lb <- c(0.13,0.13)
ub <- c(0.24,0.24)


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_simple_K5,
                del = del,
                mu = mu,
                n = n,
                sigma=sigma,
                q = q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)
x3<- 1-2* res$solution[1]-2* res$solution[2]
x3