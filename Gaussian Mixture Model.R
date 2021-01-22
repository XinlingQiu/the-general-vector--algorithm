library('abind')
seed=sample(100000,1)
set.seed(seed)
n <- 5000#????????��

#??һ????˹?ֲ?
alpha1 <- 0.4
miu1   <- 3
sigma1 <- 3

# ?ڶ?????˹?ֲ?
alpha2 <- 0.6
miu2   <- -4
sigma2 <- 2

#truepara<-c(0.4,0.6,3,-4,3,2)#??ʵ????

kk <- 6#????????
n1 <- floor(n*alpha1)#??һ????˹?ֲ?????????��
n2 <- n-n1

samp <-numeric(n)
samp[1:n1] <- rnorm(n1, miu1, sigma1)
samp[(n1+1):n] <- rnorm(n2, miu2, sigma2)

#??ͼ
hist(samp, freq = FALSE,main="??˹????ģ??",xlab=sprintf("??????��:%d",n),ylab="Ƶ??")
lines(density(samp), col = 'red')

INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
#####################EM?㷨????############################
oneiter<-function(theta){
  alpha=theta[1:(kk/3)]
  miu=theta[(kk/3+1):(2*kk/3)]
  sigma=theta[(2*kk/3+1):kk]
  
  prob <- matrix(rep(0, kk/3*n), nrow = n)
  weight <- matrix(rep(0, kk/3*n), nrow = n)
  
  # E-??
  for (i in 1:(kk/3)) {
    prob[, i]   <- sapply(samp, dnorm, miu[i], sigma[i])
    weight[, i] <- alpha[i] * prob[, i]
  }
  row_sum <- rowSums(weight)
  prob    <- weight/row_sum
  
  
  # M-??
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
}#һ??????
#####################epsilon_0?????㷨############################
eps_0<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta1-INV(a-b)*sum(b^2)
}##epsilon??????��?㷨???е?һ??????
#####################epsilon_1?????㷨############################
eps_1<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta2+INV(INV(a)-INV(b))
}##epsilon??????��?㷨???е?һ??????
eps<-function(theta1,theta2,theta3){
  theta1+INV(theta3-theta2)
}
#####################һ????epsilon?㷨??ȫ????############################
geps<-function(theta){
  l=nrow(theta)
  t<-matrix(0,l,kk)
  col=kk
  l=l-1
  while(col>0){
    for(k in 1:l){
      t[k,]=t[k+1,]+INV(theta[k+1,]-theta[k,])
    }
    l=l-1
    if(l>0){
      for(k in 1:l){
        theta[k,]=theta[k+1,]+INV(t[k+1,]-t[k,])
      }
      l=l-1
    }
    col=col-1
  }
  theta[1:(l+1),]
}##һ????epsilon??????��?㷨???е???ȫ????

###############################һ??ģ???????仯ͼ#####################################
l=300
threshold=-30
alphafirst<-runif(kk/3)
alphafirst<-alphafirst/sum(alphafirst)
miufirst<-runif(kk/3,min=-10,max=10)
sigmafirst<-runif(kk/3,min=0,max=10)##????ģ????ʼ????
thetafirst<-c(alphafirst,miufirst,sigmafirst)

diff<-matrix(0,l,3)
limit<-matrix(0,l,3)
theta1<-thetafirst#??ʼ????
theta2<-oneiter(theta1)
while(metrics(theta1,theta2)>threshold){
  theta1=theta2
  theta2=oneiter(theta1)
}
truepara=theta2
#########EM(1)
theta<-matrix(0,l,kk)
theta[1,]<-thetafirst#??ʼ????
for(step in 1:(l-1)){
  theta[step+1,]=oneiter(theta[step,])
}
for(step in (2*(kk+1)):l){
  diff[step,1]=metrics(theta[step-1,],theta[step,])
  limit[step,1]=metrics(theta[step,],truepara)
}
#########(2)
theta_1=theta
for(step in 2:(l-1)){
  theta_1[step-1,]=eps_1(theta_1[step-1,],theta_1[step,],theta_1[step+1,])
}
for(step in (2*(kk+1)):l){
  diff[step,2]=metrics(theta_1[step-2,],theta_1[step-3,])
  limit[step,2]=metrics(theta_1[step-2,],truepara)
}


##########epsilon
theta_3<-geps(theta)
for(step in (2*(kk+1)):l){
  diff[step,3]=metrics(theta_3[step-2*kk,],theta_3[step-2*kk-1,])
  limit[step,3]=metrics(theta_3[step-2*kk,],truepara)
}



diff=diff[(2*(kk+1)):l,]
limit=limit[(2*(kk+1)):l,]
plot(diff[,3],type="n",ylab="distance",xlab="k",
     ylim=c(-30,4))
lines(diff[,1],col="red",lty=1,pch=15,lwd=1)
lines(diff[,2],col="green",lty=2,pch=16,lwd=2)
lines(diff[,3],col="blue",lty=3,pch=17,lwd=3)
legend("bottomleft",c("em","delta2","epsilon"),
       col=c("red","green","black"),
       text.col=c("red","green","black"),
       pch=c(15,16,17),lty=c(1,2,3),bty="n",cex=0.6)



plot(limit[,3],type="n",ylab="limit",xlab="k",ylim=c(-30,4))
lines(limit[,1],col="red",lty=1,pch=15,lwd=1)
lines(limit[,2],col="green",lty=2,pch=16,lwd=2)
lines(limit[,3],col="blue",lty=3,pch=17,lwd=3)
legend("bottomleft",c("em","delta2","epsilon"),
       col=c("red","green","black"),
       text.col=c("red","green","black"),
       pch=c(15,16,17),lty=c(1,2,3),bty="n",cex=0.6)




####################################epsilon#######################################
start=20
l=20
cl=ceiling(l/2)
theta_5=array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_5[,,1]=0
theta_5[,,2]=theta[start:(start+l-1),]
diff5<-matrix(0,l,cl)
for(m in 3:(l+1)){
  for(j in 1:(l-m+2)){
    theta_5[j,,m]=eps(theta_5[j+1,,m-2],theta_5[j,,m-1],theta_5[j+1,,m-1])
  }
}
for(j in 1:cl){
  for(i in 1:(l-2*(j-1))){
    diff5[i,j]=metrics(theta_5[i,,2*j],truepara)
  }
}

View(diff5)


#####################????ģ??Ч?ʱȽ?##################################
nn=100
threshold=-10
iternum<-matrix(0,nn,3)##????????
time<-matrix(0,nn,3)##CPUʱ??
for(i in 1:nn){
  alphafirst<-runif(kk/3)
  alphafirst<-alphafirst/sum(alphafirst)
  miufirst<-runif(kk/3,min=-10,max=10)
  sigmafirst<-runif(kk/3,min=0,max=10)##????ģ????ʼ????
  thetafirst<-c(alphafirst,miufirst,sigmafirst)
  
  #############################em?㷨#################################################
  d=proc.time()
  l=6##??ʼ????????
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,]=oneiter(theta[step,])
  }
  diff<-metrics(theta[l,],theta[l-1,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=oneiter(theta[l-1,])
    diff=metrics(theta[l-1,],theta[l,])
  }
  iternum[i,1]=l
  time[i,1]<-(proc.time()-d)[1]
  
  ##############################eps_1?????㷨#######################################
  d=proc.time()
  l=6
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  theta[2,]<-oneiter(theta[1,])
  for(step in 2:(l-1)){
    theta[step+1,]=oneiter(theta[step,])
    theta[step-1,]=eps_1(theta[step-1,],theta[step,],theta[step+1,])
  }
  diff=metrics(theta[l-2,],theta[l-3,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=oneiter(theta[l-1,])
    theta[l-2,]=eps_1(theta[l-2,],theta[l-1,],theta[l,])
    diff=metrics(theta[l-2,],theta[l-3,])
  }
  iternum[i,2]=l
  time[i,2]<-(proc.time()-d)[1]
  

  ############################general epsilon#######################################
  d=proc.time()
  l=2*kk+2
  theta<-array(rep(0,kk*l*(2*kk+2)),dim=c(l,kk,2*kk+2))
  theta[,,1]=0
  theta[1,,2]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,,2]=oneiter(theta[step,,2])
  }
  for(j in 3:(2*kk+2)){
    for(step in 1:(l+2-j)){
      theta[step,,j]=eps(theta[step+1,,j-2],theta[step,,j-1],theta[step+1,,j-1])
    }
  }
  diff=metrics(theta[1,,2*kk+2],theta[2,,2*kk+2])
  while(diff>threshold){
    l=l+1
    theta=abind(theta,array(rep(0,(2*kk+2)*kk),dim=c(1,kk,2*kk+2)),along=1)
    theta[l,,1]=0
    theta[l,,2]=oneiter(theta[l-1,,2])
    for(k in 3:(2*kk+2)){
      L=l+2-k
      theta[L,,k]=eps(theta[L+1,,k-2],theta[L,,k-1],theta[L+1,,k-1])
    }
    diff=metrics(theta[l-2*kk,,(2*kk+2)],theta[l-2*kk-1,,(2*kk+2)])
    }
  iternum[i,3]=l
  time[i,3]<-(proc.time()-d)[1]
}

par(mfrow = c(2, 3))
hist(iternum[,1],freq=TRUE,xlab="k",ylab="frequency",main="em",col="red")
hist(iternum[,2],freq=TRUE,xlab="k",ylab="frequency",main="delta2",col="green")
hist(iternum[,3],freq=TRUE,xlab="k",ylab="frequency",main="epsilon",col="black")
hist(time[,1],freq=TRUE,xlab="k",ylab="time",main="em",col="red")
hist(time[,2],freq=TRUE,xlab="k",ylab="time",main="delta2",col="green")
hist(time[,3],freq=TRUE,xlab="k",ylab="time",main="epsilon",col="black")


summary(iternum)
summary(time)









