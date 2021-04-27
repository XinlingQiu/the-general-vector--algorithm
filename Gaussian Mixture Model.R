library('abind')
seed=sample(100000,1)
set.seed(seed)
n <- 5000


alpha1 <- 0.4
miu1   <- 3
sigma1 <- 3


alpha2 <- 0.6
miu2   <- -4
sigma2 <- 2


kk <- 6
n1 <- floor(n*alpha1)
n2 <- n-n1

samp <-numeric(n)
samp[1:n1] <- rnorm(n1, miu1, sigma1)
samp[(n1+1):n] <- rnorm(n2, miu2, sigma2)

INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
##################### em ############################
em<-function(theta){
  alpha=theta[1:(kk/3)]
  miu=theta[(kk/3+1):(2*kk/3)]
  sigma=theta[(2*kk/3+1):kk]
  
  prob <- matrix(rep(0, kk/3*n), nrow = n)
  weight <- matrix(rep(0, kk/3*n), nrow = n)
  
  # E-step
  for (i in 1:(kk/3)) {
    prob[, i]   <- sapply(samp, dnorm, miu[i], sigma[i])
    weight[, i] <- alpha[i] * prob[, i]
  }
  row_sum <- rowSums(weight)
  prob    <- weight/row_sum
  
  
  # M-step
  for (j in 1:(kk/3)) {
    sum1     <- sum(prob[, j])
    sum2     <- sum(samp*prob[, j])
    alpha[j] <- sum1/n
    miu[j]   <- sum2/sum1
    sum3     <- sum(prob[, j]*(samp-miu[j])^2)
    sigma[j] <- sqrt(sum3/sum1)
  }
  theta<-c(alpha,miu,sigma)
  theta
}
##################### delta2 ############################
delta2<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta2+INV(INV(a)-INV(b))
}
##################### epsilon for one step ############################
eps<-function(theta1,theta2,theta3){
  theta1+INV(theta3-theta2)
}
############################### this is for figure 1 #####################################
l=500
threshold=-30
alphafirst<-runif(kk/3)
alphafirst<-alphafirst/sum(alphafirst)
miufirst<-runif(kk/3,min=-10,max=10)
sigmafirst<-runif(kk/3,min=0,max=10)
thetafirst<-c(alphafirst,miufirst,sigmafirst)

diff<-matrix(0,l,3)
limit<-matrix(0,l,3)
theta1<-thetafirst
theta2<-em(theta1)
while(metrics(theta1,theta2)>threshold){
  theta1=theta2
  theta2=em(theta1)
}
truepara=theta2
#########EM
theta<-matrix(0,l,kk)
theta[1,]<-thetafirst
for(step in 1:(l-1)){
  theta[step+1,]=em(theta[step,])
}
for(step in 2:(l-2)){
  diff[step,1]=metrics(theta[step-1,],theta[step,])
  limit[step,1]=metrics(theta[step,],truepara)
}
#########the epsilon accelerated EM algorithm
theta_1=theta
for(step in 2:(l-1)){
  theta_1[step-1,]=delta2(theta_1[step-1,],theta_1[step,],theta_1[step+1,])
}
for(step in 2:(l-2)){
  diff[step,2]=metrics(theta_1[step-1,],theta_1[step,])
  limit[step,2]=metrics(theta_1[step,],truepara)
}


###########the general epsilon accelerated EM algorithm
theta_3<-array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_3[,,2]<-theta
diff[2,3]=metrics(theta_3[2,,2],theta_3[1,,2])
limit[2,3]=metrics(theta_3[2,,2],truepara)
diff[3,3]=metrics(theta_3[2,,2],theta_3[3,,2])
limit[3,3]=metrics(theta_3[3,,2],truepara)
for(j in 3:5){
  for(step in 1:(4-(j-2))){
    theta_3[step,,j]=eps(theta_3[step+1,,j-2],theta_3[step,,j-1],theta_3[step+1,,j-1])
  }
  
}
if(metrics(theta_3[3,,2],theta_3[4,,2])<metrics(theta_3[1,,4],theta_3[2,,4])){
  diff[4,3]=metrics(theta_3[3,,2],theta_3[4,,2])
  limit[4,3]=metrics(theta_3[4,,2],truepara)
}else{
  diff[4,3]=metrics(theta_3[1,,4],theta_3[2,,4])
  limit[4,3]=metrics(theta_3[2,,4],truepara)
}
for(step in 5:l){
  diff[step,3]=metrics(theta_3[step,,2],theta_3[step-1,,2])
  limit[step,3]=metrics(theta_3[step,,2],truepara)
  for(k in 3:(step+1)){
    L=step-k+2
    theta_3[L,,k]=eps(theta_3[L+1,,k-2],theta_3[L,,k-1],theta_3[L+1,,k-1])
    if(k%%2==0&&(L>1)){
      diff_2=metrics(theta_3[L,,k],theta_3[L-1,,k])
      limit_2=metrics(theta_3[L,,k],truepara)
      if(is.finite(diff_2)&&diff_2<diff[step,3]){diff[step,3]=diff_2}
      if(is.finite(limit_2)&&limit_2<limit[step,3]){limit[step,3]=limit_2}
      
    }
  }
}



plot(diff[2:(l-4),3],type="n",ylab="m(¦È(n),¦È(n+1))",xlab="number of iterations",
     ylim=c(-30,4))
lines(diff[2:(l-4),1],col="red",lwd=3)
lines(diff[2:(l-4),2],col="green",lwd=3)
lines(diff[2:(l-4),3],col="black",lwd=3)
legend("topright",c("EM algorithm","The ¦Å-accelerated EM algorithm","the general ¦Å-accelerated EM algorithm"),
       col=c("red","green","black"),
       text.col=c("red","green","black"),
       lty=c(1,1,1),
       bty="n",cex=1.2)




plot(limit[2:(l-4),3],type="n",ylab="m(¦È(n),¦È*)",xlab="number of iterations",ylim=c(-30,4))
lines(limit[2:(l-4),1],col="red",lwd=3)
lines(limit[2:(l-4),2],col="green",lwd=3)
lines(limit[2:(l-4),3],col="black",lwd=3)
legend("topright",c("EM algorithm","The ¦Å-accelerated EM algorithm","the general ¦Å-accelerated EM algorithm"),
       col=c("red","green","black"),
       text.col=c("red","green","black"),
       lty=c(1,1,1),
       bty="n",cex=1.2)


#################################### figure 1 end #########################################

############################### this is for tabel 3 ##################################
nn=10
threshold=-10
iternum<-matrix(0,nn,3)
time<-matrix(0,nn,3)
for(i in 1:nn){
  alphafirst<-runif(kk/3)
  alphafirst<-alphafirst/sum(alphafirst)
  miufirst<-runif(kk/3,min=-10,max=10)
  sigmafirst<-runif(kk/3,min=0,max=10)
  thetafirst<-c(alphafirst,miufirst,sigmafirst)
  ############################# em #################################################
  d=proc.time()
  l=6
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,]=em(theta[step,])
  }
  diff<-metrics(theta[l,],theta[l-1,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=em(theta[l-1,])
    diff=metrics(theta[l-1,],theta[l,])
  }
  iternum[i,1]=l
  time[i,1]<-(proc.time()-d)[1]
  
  ############################# delta2 #######################################
  d=proc.time()
  l=6
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  theta[2,]<-em(theta[1,])
  for(step in 2:(l-1)){
    theta[step+1,]=em(theta[step,])
    theta[step-1,]=delta2(theta[step-1,],theta[step,],theta[step+1,])
  }
  diff=metrics(theta[l-2,],theta[l-3,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=em(theta[l-1,])
    theta[l-2,]=delta2(theta[l-2,],theta[l-1,],theta[l,])
    diff=metrics(theta[l-2,],theta[l-3,])
  }
  iternum[i,2]=l
  time[i,2]<-(proc.time()-d)[1]
  
  ############################general epsilon#######################################
  d=proc.time()
  l=4
  theta<-array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
  theta[,,1]=0
  theta[1,,2]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,,2]=em(theta[step,,2])
  }
  for(j in (l-1):(l+1)){
    for(step in 1:(l-(j-2))){
      theta[step,,j]=eps(theta[step+1,,j-2],theta[step,,j-1],theta[step+1,,j-1])
    }
  }
  if(metrics(theta[1,,4],theta[2,,4])<metrics(theta[3,,2],theta[4,,2])){
    diff=metrics(theta[1,,4],theta[2,,4])
  }else{
    diff=metrics(theta[3,,2],theta[4,,2])
  }
  while(diff>threshold){
    l=l+1
    theta=abind(theta,array(rep(0,kk*(l-1)),dim=c(l-1,kk,1)),along=3)
    theta=abind(theta,array(rep(0,(l+1)*kk),dim=c(1,kk,l+1)),along=1)
    theta[l,,1]=0
    for(k in 2:(l+1)){
      L=l-(k-2)
      if(k==2){
        theta[l,,2]=em(theta[l-1,,2])
        diff=metrics(theta[l,,2],theta[l-1,,2])
      }else{
        theta[L,,k]=eps(theta[L+1,,k-2],theta[L,,k-1],theta[L+1,,k-1])
        if(k%%2==0&&(L>1)){
          diff_2=metrics(theta[L,,k],theta[L-1,,k])
          if(is.finite(diff_2)&&diff_2<diff){diff=diff_2;
          }
        }
      }
    }
  }
  iternum[i,3]=l
  time[i,3]<-(proc.time()-d)[1]
}

summary(iternum)
summary(time)

#################################### table 3 end #########################################




