###########################################################################
##### Simulation study to check the performance of the Max--min design for
#####  UIT Path graph versus the Balanced design for K=3
###########################################################################

# load the following libraries


library('nloptr')
library(mvtnorm)
rm(list = ls())

# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.
psi<- function(design,theta,n,sigma,q)
{
  
  t1<- design
  t2<- 1-2*t1
  q<- qnorm(1-al)
  de12 <- ((theta[1]-theta[2])*n)/sigma
  de23 <- ((theta[2]-theta[3])*n)/sigma
  mu1 <- (sqrt((t1*t2)/(t1 + t2)))*de12
  mu2 <- (sqrt((t2*t1)/(t2 + t1)))*de23
  
  a12 <- t1*t1
  b12 <- (t1+t2)*(t2+t1)
  rho12 <- -sqrt(a12/b12) 
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  S <-  rbind(c(1,rho12),c(rho12,1))
  m <- c(mu1,mu2)
  y <- -(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1))
  return(y[1])
}

# Function for calculating the average power B

B<- function(design,theta_vals,prob_vals,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],n,sigma,q)
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for calculating optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,n,sigma,q)
{
  des_0 <-0.18
  lb <- 0
  ub <- 0.4
  local_opts <- list( "algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_LN_COBYLA",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  
  res <- nloptr ( x0 = des_0,
                  theta_vals = theta_vals,
                  prob_vals = prob_vals,
                  n=n,
                  sigma=sigma,
                  q=q,
                  eval_f = B,
                  lb = lb,
                  ub = ub,
                  opts = opts
  )
  
  return(res$solution)
  
}

#Step 1

step_1<- function(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q){
  
  psi_min <- 100      
  K<- 3
  this_psi_neg<- rep(0, K-1)
  this_psi<- rep(0, K-1)
  this_psi_neg[1] <-  psi(design_l,theta_vals_l[1,],n,sigma,q)
  this_psi[1] <- -this_psi_neg[1]
  this_psi_neg[2] <-  psi(design_l,theta_vals_l[2,],n,sigma,q)
  this_psi[2] <- -this_psi_neg[2]
  
  this_psi<- as.numeric(round(this_psi, 3))
  
  this_psi_min<- min(this_psi)
  print(this_psi)
  index <- which(this_psi == this_psi_min)
  if(this_psi_min < psi_min)
  {
    psi_min <- this_psi_min
  }
  print("iteration_num")
  print(iteration_num)
  iteration_num <- iteration_num+1      
  psi_min
  B_pi_l
  
  # Stopping Criterion 
  
  if ((psi_min < -B_pi_l)&(abs(psi_min+B_pi_l)>1.0e-4)){
    stopping <- 0}
  else{
    stopping<- 1}
  
  
  if (stopping==1){
    print("Found minmax design")
    print(design_l)
    print("Least Favourable Distribution")
    print("probs")
    print(prob_vals_l)
    print("thetas")
    print(theta_vals_l)}
  else{
    step_2(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
  }
}

# Step 2 

step_2<- function(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q){
  #delta_l = unit mass to theta_l
  step_3(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
}

# Step 3

step_3<- function(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
{
  design_l1 <- design_l
  smallest_B_pi_t_l1 <- -100
  prob_vals_l1 <- prob_vals_l
  theta_vals_l1 <- theta_vals_l
  K=3
  theta_vals_l
  prob_vals_l
  
  prob_vals_t_l1<- matrix(0, length(H), K-1)
  delta_l <- rep(0, 2)
  for(i in 1:length(index)){
    delta_l[index[i]] <- 1/length(index)
  }
  
  for (t_val in 1:length(H)){
    #new T priors
    prob_vals_t_l1[t_val,] =  ((1-H[t_val])*prob_vals_l)+(H[t_val]*delta_l) 
  }
  
  # Assign l+1 values to l.
  for(i in 1:length(H)){
    design_t_l1 <- get_optimal_on_the_average_design(theta_vals_l, prob_vals_t_l1[i,],n,sigma,q)
    B_pi_t_l1 <- B(design_t_l1, theta_vals_l,prob_vals_t_l1[i,],n,sigma,q)  
    if (-smallest_B_pi_t_l1 > -B_pi_t_l1){
      smallest_B_pi_t_l1 <- B_pi_t_l1
      design_l1 <- design_t_l1
      prob_vals_l1 <- prob_vals_t_l1[i,]
      theta_vals_l1 <- theta_vals_l
    }
  }
  print("prob_vals_l1")
  print(prob_vals_l1)
  print("design_l1")
  print(design_l1)
  
  step_4(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,n,sigma,q)
}

# Step 4
step_4<- function(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,n,sigma,q){
  B_pi_l1 <- B(design_l1, theta_vals_l1, prob_vals_l1,n,sigma,q)
  print(B_pi_l1-B_pi_l)
  
  # Assigning B(l+1) values to B(l)  
  
  if (-B_pi_l1 <= -B_pi_l)
  {
    B_pi_l <- B_pi_l1
    print("B_pi_l")
    print(B_pi_l)
    design_l <- design_l1
    #"Assigned l1 to l"
    theta_vals_l <- theta_vals_l1
    prob_vals_l <- prob_vals_l1 
    # return
    step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
  }
  else{
    h_grid_space <- h_grid_space/2
    H <-seq(0, 1, by=h_grid_space)
    
    K<- 3
    this_psi_neg_1<- rep(0, K-1)
    this_psi_1<- rep(0, K-1)
    this_psi_neg_1[1] <-  psi(design_l,theta_vals_l[1,])
    this_psi_1[1] <- -this_psi_neg_1[1]
    this_psi_neg_1[2] <-  psi(design_l,theta_vals_l[2,])
    this_psi_1[2] <- -this_psi_neg_1[2]
    this_psi_1<- as.numeric(round(this_psi_1, 2))
    
    this_psi_min_1<- min(this_psi_1)
    index_1 <- which(this_psi_1 == this_psi_min_1)
    
    step_3(index_1,H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,n,sigma,q)
  }
}

N<-300                # Total sample size
n <- sqrt(N)
sigma <- 14.57        # standard deviation
al <- 0.05/4          # significance level alpha
q <- qnorm(1-al)      # critical value
del<-0.7                   # delta 
K<-3
iteration_num<-0          
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.2       # grid length
H<- seq(0, 1, by=h_grid_space)     # grid vector


print( "Initial Prior on Theta : ")
p<-0.5
# Theta values (LFC's) are their corresponding prior probabilities
theta_vals_l <-matrix(0, K-1, K)
for (i in 1:K-1){
  for (j in 1:i){
    theta_vals_l[i,j] <- del
  }}

prob_vals_l <- c(p,1-p)

# Starting optimum on the average design and corresponding power 
print( "Starting Design") 
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,n,sigma,q)
print( B_pi_l)
step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)



# values of the parameters associated


m<-1e3              # number of simulations


eta <- design_l           # max--min design for path UIT K=3
x_mm<-c(eta,1-2*eta,eta)
K<-3
true.mean<-c(39.5,35.2,39.4)   # true mean from real data
psd <- 14.57                   # pooled standard deviation of the real data


n_mm<- rep(0,K)
for (i in 1:K) {
  n_mm[i]<- ceiling(x_mm[i]*N)   # exact max--min design
}
mean_smm<-0
pwr_mm<- rep(0,m)  

n_b<- rep(0,K)
for (i in 1:K) {
  n_b[i] <- ceiling( N/K)   # exact balanced design
}

mean_b<-0
pwr_b<- rep(0,m)  

# generate from normal distribution

for (j in 1:m){
  
  for (i in 1:K) {
    
    # Simulations for Min--Max design
    sample_mm <- rep(0,n_mm[i])
    sample_mm <- rnorm(n_mm[i], mean=true.mean[i], psd)
    mean_smm[i] <- mean(sample_mm)
    
    
    # Simulations for Balanced design
    sample_b <- rep(0,n_b[i])
    sample_b <- rnorm(n_b[i], mean=true.mean[i], psd)
    mean_b[i] <- mean(sample_b)
  }
  
  pwr_mm[j] <- psi(eta,mean_smm,n,sigma,q)  # power vector for the simulated data for the min--max design
  pwr_b[j] <- psi(1/K,mean_b,n,sigma,q)     # power vector for the simulated data for the balanced design
  
  
}


print('Median of power for the max--min design')

print(median(-pwr_mm))

print('Median of power for the balanced design')

print(median(-pwr_b))