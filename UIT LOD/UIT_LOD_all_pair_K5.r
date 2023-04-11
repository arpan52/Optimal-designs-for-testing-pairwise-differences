
# for pairs (1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)

#install.packages("nloptr")
#install.packages("mvtnorm")
library('nloptr')
library(mvtnorm)

UIT_Complete_CDF_K5_power_singular <- function(x,mu,del,n,q){
  
  # Calculating mean and covariance matrix of Z_ij's

  de12 <- ((mu[1]-mu[2])*n)/sigma
  de23 <- ((mu[2]-mu[3])*n)/sigma
  de34 <- ((mu[3]-mu[4])*n)/sigma
  de45 <- ((mu[4]-mu[5])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de15 <- ((mu[1]-mu[5])*n)/sigma
  de24 <- ((mu[2]-mu[4])*n)/sigma
  de25 <- ((mu[2]-mu[5])*n)/sigma
  de35 <- ((mu[3]-mu[5])*n)/sigma
  
  
  
  a12 <- x[1]*x[3]
  b12 <- (x[1]+x[2])*(x[2]+x[3])
  rho12 <- -sqrt(a12/b12)
  
  rho13 <- 0
  rho14 <- 0
  
  a15 <- x[2]*x[3]
  b15 <- (x[1]+x[2])*(x[1]+x[3])
   rho15 <- sqrt(a15/b15)
  
  a16 <- x[2]*x[4]
  b16 <- (x[1]+x[2])*(x[1]+x[4])
  rho16 <- sqrt(a16/b16)
  
  a17 <- x[2]*x[5]
  b17 <- (x[1]+x[2])*(x[1]+x[5])
  rho17 <- sqrt(a17/b17)
  
  a18 <- x[1]*x[4]
  b18 <- (x[1]+x[2])*(x[2]+x[4])
  rho18 <- -sqrt(a18/b18)
  
  a19 <- x[1]*x[5]
  b19 <- (x[1]+x[2])*(x[2]+x[5])
  rho19 <- -sqrt(a19/b19)
  
  rho110<- 0
  
  a23 <- x[2]*x[4]
  b23 <- (x[2]+x[3])*(x[3]+x[4])
   rho23 <- -sqrt(a23/b23)
  
  rho24<- 0
  
  a25 <- x[1]*x[2]
  b25 <- (x[2]+x[3])*(x[1]+x[3])
  rho25 <- sqrt(a25/b25)
  
  rho26<- 0
  
  rho27<- 0
  
  a28 <- x[3]*x[4]
  b28 <- (x[2]+x[3])*(x[2]+x[4])
  rho28 <- sqrt(a28/b28)
  
  a29 <- x[3]*x[5]
  b29 <- (x[2]+x[3])*(x[2]+x[5])
  rho29 <- sqrt(a29/b29)
  
  a210 <- x[2]*x[5]
  b210 <- (x[2]+x[3])*(x[3]+x[5])
  rho210 <- -sqrt(a210/b210)
  
  a34 <- x[3]*x[5]
  b34 <- (x[3]+x[4])*(x[4]+x[5])
  rho34 <- -sqrt(a34/b34)
  
  a35 <- x[1]*x[4]
  b35 <- (x[3]+x[4])*(x[1]+x[3])
  rho35 <- -sqrt(a35/b35)
  
  a36 <- x[1]*x[3]
  b36 <- (x[3]+x[4])*(x[1]+x[4])
  rho36 <- sqrt(a36/b36)
  
  rho37<- 0
  
  a38 <- x[2]*x[3]
  b38 <- (x[3]+x[4])*(x[2]+x[4])
  rho38 <- sqrt(a38/b38)
  
  rho39<- 0
  
  a310 <- x[4]*x[5]
  b310 <- (x[3]+x[4])*(x[3]+x[5])
  rho310 <- sqrt(a310/b310)
  
  
  rho45<-0
  
  a46 <- x[1]*x[5]
  b46 <- (x[4]+x[5])*(x[1]+x[4])
  rho46 <- -sqrt(a46/b46)
  
  a47 <- x[1]*x[4]
  b47 <- (x[4]+x[5])*(x[1]+x[5])
  rho47 <- sqrt(a47/b47)
  
  a48 <- x[2]*x[5]
  b48 <- (x[4]+x[5])*(x[2]+x[4])
  rho48 <- -sqrt(a48/b48)
  
  a49 <- x[2]*x[4]
  b49 <- (x[4]+x[5])*(x[2]+x[5])
  rho49 <- sqrt(a49/b49)
  
  a410 <- x[3]*x[4]
  b410 <- (x[4]+x[5])*(x[3]+x[5])
  rho410 <- sqrt(a410/b410)
  
  a56 <- x[3]*x[4]
  b56 <- (x[1]+x[3])*(x[1]+x[4])
  rho56 <- sqrt(a56/b56)
  
  a57 <- x[3]*x[5]
  b57 <- (x[1]+x[3])*(x[1]+x[5])
  rho57 <- sqrt(a57/b57)
  
  rho58<- 0
  rho59<- 0
  
  a510 <- x[1]*x[5]
  b510 <- (x[1]+x[3])*(x[3]+x[5])
  rho510 <- -sqrt(a510/b510)
  
  a67 <- x[4]*x[5]
  b67 <- (x[1]+x[4])*(x[1]+x[5])
  rho67 <- sqrt(a67/b67)
  
  a68<- x[1]*x[2]
  b68 <- (x[1]+x[4])*(x[2]+x[4])
  rho68 <- sqrt(a68/b68)
  
  rho69<-0
  
  rho610<-0
  
  rho78<-0
  
  a79 <- x[1]*x[2]
  b79 <- (x[1]+x[5])*(x[2]+x[5])
  rho79 <- sqrt(a79/b79)
  
  a710 <- x[1]*x[3]
  b710 <- (x[1]+x[5])*(x[3]+x[5])
  rho710 <- sqrt(a710/b710)
  
  a89 <- x[4]*x[5]
  b89 <- (x[2]+x[4])*(x[2]+x[5])
  rho89 <- sqrt(a89/b89)
  
  rho810<- 0
  
  a910 <- x[2]*x[3]
  b910 <- (x[2]+x[5])*(x[3]+x[5])
  rho910 <- sqrt(a910/b910)
  
  
  mu1 <- (sqrt((x[1]*x[2])/(x[1] + x[2])))*de12
  mu2 <- (sqrt((x[2]*x[3])/(x[2] + x[3])))*de23
  mu3 <- (sqrt((x[3]*x[4])/(x[3] + x[4])))*de34
  mu4 <- (sqrt((x[4]*x[5])/(x[4] + x[5])))*de45
  mu5 <- (sqrt((x[1]*x[3])/(x[1] + x[3])))*de13
  mu6 <- (sqrt((x[1]*x[4])/(x[1] + x[4])))*de14
  mu7 <- (sqrt((x[1]*x[5])/(x[1] + x[5])))*de15
  mu8 <- (sqrt((x[2]*x[4])/(x[2] + x[4])))*de24
  mu9 <- (sqrt((x[2]*x[5])/(x[2] + x[5])))*de25
  mu10 <- (sqrt((x[3]*x[5])/(x[3] + x[5])))*de35
  
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10)
  S <-  rbind(c(1,rho12,rho13,rho14,rho15,rho16,rho17,rho18,rho19,rho110),c(rho12,1,rho23,rho24,rho25,rho26,rho27,rho28,rho29,rho210),c(rho13,rho23,1,rho34,rho35,rho36,rho37,rho38,rho39,rho310),c(rho14,rho24,rho34,1,rho45,rho46,rho47,rho48,rho49,rho410),c(rho15,rho25,rho35,rho45,1,rho56,rho57,rho58,rho59,rho510),c(rho16,rho26,rho36,rho46,rho56,1,rho67,rho68,rho69,rho610),c(rho17,rho27,rho37,rho47,rho57,rho67,1,rho78,rho79,rho710),c(rho18,rho28,rho38,rho48,rho58,rho68,rho78,1,rho89,rho810),c(rho19,rho29,rho39,rho49,rho59,rho69,rho79,rho89,1,rho910),c(rho110,rho210,rho310,rho410,rho510,rho610,rho710,rho810,rho910,1))
  
  lw1 <- c(-q,-q,-q,-q,-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q,q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))    #Power
  
  
}

mu<- c(0,0,-0.25,0.25,0)      # mu vector
sigma <- 1
N<- 180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/20                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6                      # delta


# equality Constraints for x
const1 <- function(x,mu,del,n,q) {
  return( (x[1]+x[2]+x[3]+x[4]+x[5]-1))
}

# Lower and upper bounds
lb <- c(0,0,0,0,0)
ub <- c(1,1,1,1,1)

#initial values
x0 <- c(0.1,0.2,0.3,0.2,0.1)


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = UIT_Complete_CDF_K5_power_singular,
                mu=mu,
                del=del,
                n=n,
                q=q,
                lb = lb,
                ub = ub,
                eval_g_eq = const1,
                opts = opts
)
print(res)

