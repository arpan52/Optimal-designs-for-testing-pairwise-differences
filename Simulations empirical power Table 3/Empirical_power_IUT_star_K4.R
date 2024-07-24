###############################################################################################################
###### Simulations to compute the Empirical power for STAR graph for IUT with K=4 experimental groups  ########
###############################################################################################################

rm(list = ls())
true.mean <- c(1.2128,1.0794,1.1307,1.2763)              # Sample mean of example
K<-length(true.mean)                                     # No. of experimental groups
psd<- 0.1140                                             # Poolad standard deviation

n_m<- c(36,22,22,22)  # max-min design
n_b= c(25,25,26,26)   # balance design
n_d<- c(37,22,22,21)  # max-min design


m<-10000                                                # No. of simulations

# Intializations 

X.m_mean<-rep(0,K)
X.b_mean<-rep(0,K)
X.d_mean<-rep(0,K)

Y.m_mean<-rep(0,K)
Y.b_mean<-rep(0,K)
Y.d_mean<-rep(0,K)


t_stat_m<-rep(0,m)
t_stat_b<-rep(0,m)
t_stat_d<-rep(0,m)


ctr_m<-0
ctr_b<-0
ctr_d<-0


for (i in 1:m) {  #simulating from max-min design
  for (j in 1:K){
    X.m <- rnorm(n_m[j],true.mean[1],psd)  # null data from max--min design
    X.m_mean[j] <- mean(X.m)
    
    X.b <- rnorm(n_b[j],true.mean[1],psd)  # null data from balanced design
    X.b_mean[j] <- mean(X.b)
    
    X.d <- rnorm(n_d[j],true.mean[1],psd)  # null data from dunnett's design
    X.d_mean[j] <- mean(X.d)
  }
  
  Z_12_m = (X.m_mean[1]-X.m_mean[2])/(psd*sqrt((1/n_m[1])+(1/n_m[2])))
  Z_13_m = (X.m_mean[1]-X.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_14_m = (X.m_mean[1]-X.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  
  t_stat_m[i]<- min(c(abs(Z_12_m),abs(Z_13_m),abs(Z_14_m)))                          #test statistics under the null for max-min design
  
  Z_12_b = (X.b_mean[1]-X.b_mean[2])/(psd*sqrt((1/n_b[1])+(1/n_b[2])))
  Z_13_b = (X.b_mean[1]-X.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_14_b = (X.b_mean[1]-X.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  
  t_stat_b[i]<- min(c(abs(Z_12_b),abs(Z_13_b),abs(Z_14_b)))                          #test statistics under the null for balanced design
  
  Z_12_d = (X.d_mean[1]-X.d_mean[2])/(psd*sqrt((1/n_d[1])+(1/n_d[2])))
  Z_13_d = (X.d_mean[1]-X.d_mean[3])/(psd*sqrt((1/n_d[1])+(1/n_d[3])))
  Z_14_d = (X.d_mean[1]-X.d_mean[4])/(psd*sqrt((1/n_d[1])+(1/n_d[4])))
  
  t_stat_d[i]<- min(c(abs(Z_12_d),abs(Z_13_d),abs(Z_14_d)))                          #test statistics under the null for dunnett's design
  
}

q_m<- quantile(t_stat_m,probs = 0.95)                                                # critical value of max-min design
q_b<- quantile(t_stat_b,probs = 0.95)                                                # critical value of balanced design 
q_d<- quantile(t_stat_d,probs = 0.95)                                                # critical value of dunnett's design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
    
    Y.d <- rnorm(n_d[j],true.mean[j],psd)  # alternative data from dunett's design
    Y.d_mean[j] <- mean(Y.d)
  }
  Z_12_m.A = (Y.m_mean[1]-Y.m_mean[2])/(psd*sqrt((1/n_m[1])+(1/n_m[2])))
  Z_13_m.A = (Y.m_mean[1]-Y.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_14_m.A = (Y.m_mean[1]-Y.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  
  t_stat_m.A<- min(c(abs(Z_12_m.A),abs(Z_13_m.A),abs(Z_14_m.A)))                          #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_12_b.A = (Y.b_mean[1]-Y.b_mean[2])/(psd*sqrt((1/n_b[1])+(1/n_b[2])))
  Z_13_b.A = (Y.b_mean[1]-Y.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_14_b.A = (Y.b_mean[1]-Y.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  
  t_stat_b.A<- min(c(abs(Z_12_b.A),abs(Z_13_b.A),abs(Z_14_b.A)))                          #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
  Z_12_d.A = (Y.d_mean[1]-Y.d_mean[2])/(psd*sqrt((1/n_d[1])+(1/n_d[2])))
  Z_13_d.A = (Y.d_mean[1]-Y.d_mean[3])/(psd*sqrt((1/n_d[1])+(1/n_d[3])))
  Z_14_d.A = (Y.d_mean[1]-Y.d_mean[4])/(psd*sqrt((1/n_d[1])+(1/n_d[4])))
  
  t_stat_d.A<- min(c(abs(Z_12_d.A),abs(Z_13_d.A),abs(Z_14_d.A)))                          #test statistics under the altenative for dunnett's design
  
  if (t_stat_d.A>q_d){
    ctr_d <- ctr_d+1
  }
  
}


print('CV max-min, CV balanced, CV dunnett')
print(c(q_m,q_b,q_d))

print('power max-min, power balanced, power dunnett')
print(c(ctr_m/m,ctr_b/m,ctr_d/m))


