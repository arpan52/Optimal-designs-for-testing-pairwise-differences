###############################################################################################################
###### Simulations to compute the Empirical power for PATH graph for UIT with K=3 experimental groups  ########
###############################################################################################################

rm(list = ls())
true.mean <- c(1.3,1.6,0.9)                                # Sample mean of example
K<-length(true.mean)                                     # No. of experimental groups
psd<- 2.93                                               # Poolad standard deviation

n_m<- c(62,85,62)  # max-min design
n_b= c(70,70,69)   # balance design


m<-10000                                                  # No. of simulations

# Intializations 


X.m_mean<-rep(0,K)
X.b_mean<-rep(0,K)
Y.m_mean<-rep(0,K)
Y.b_mean<-rep(0,K)

t_stat_m<-rep(0,m)
t_stat_b<-rep(0,m)


ctr_m<-0
ctr_b<-0


for (i in 1:m) {  #simulating from max-min design
  for (j in 1:K){
    X.m <- rnorm(n_m[j],true.mean[1],psd)  # null data from max--min design
    X.m_mean[j] <- mean(X.m)
    
    X.b <- rnorm(n_b[j],true.mean[1],psd)  # null data from balanced design
    X.b_mean[j] <- mean(X.b)
  }
  Z_12_m = (X.m_mean[1]-X.m_mean[2])/(psd*sqrt((1/n_m[1])+(1/n_m[2])))
  Z_23_m = (X.m_mean[2]-X.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  t_stat_m[i]<- max(c(abs(Z_12_m),abs(Z_23_m)))                          #test statistics under the null for max-min design
  
  Z_12_b = (X.b_mean[1]-X.b_mean[2])/(psd*sqrt((1/n_b[1])+(1/n_b[2])))
  Z_23_b = (X.b_mean[2]-X.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  t_stat_b[i]<- max(c(abs(Z_12_b),abs(Z_23_b)))                          #test statistics under the null for balanced design
  
}

q_m<- quantile(t_stat_m,probs = 0.95)                                    # critical value of max-min design
q_b<- quantile(t_stat_b,probs = 0.95)                                    # critical value of balanced design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  Z_12_m.A = (Y.m_mean[1]-Y.m_mean[2])/(psd*sqrt((1/n_m[1])+(1/n_m[2])))
  Z_23_m.A = (Y.m_mean[2]-Y.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  t_stat_m.A<- max(c(abs(Z_12_m.A),abs(Z_23_m.A)))   #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_12_b.A = (Y.b_mean[1]-Y.b_mean[2])/(psd*sqrt((1/n_b[1])+(1/n_b[2])))
  Z_23_b.A = (Y.b_mean[2]-Y.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  t_stat_b.A<- max(c(abs(Z_12_b.A),abs(Z_23_b.A)))                          #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))


