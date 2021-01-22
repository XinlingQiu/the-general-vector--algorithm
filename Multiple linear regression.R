library("MASS")
library("abind")
#ģ??????
seed=sample(100000,1)
set.seed(seed)
N=100
p=10
kk=p+1#??????��
X<-matrix(0,N,kk)
X[,1]=1
S <- toeplitz((p:1)/p)
R <- rWishart(1, p, S)
R=R[,,1]
miu=matrix(0,p,1)
X[,2:(p+1)]=mvrnorm(N,miu,R)
inv=solve(t(X)%*%X)
truepara_1<-rnorm(kk,0,1)
Y=X%*%truepara_1
def=sample(1:N,N/2)##ȱʧ????
Y[def]=NA


INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
#####################EM?�?????############################
oneiter<-function(theta){
  theta=matrix(theta,kk,1)
  #E??
  for(i in def){
    Y[i]=X[i,]%*%theta
  }
  ##M??
  theta=inv%*%t(X)%*%Y
  theta=matrix(theta,1,kk)
  theta
}#һ??????
#####################epsilon_0?????�?############################
eps_0<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta1-INV(a-b)*sum(b^2)
}##epsilon??????��?�????е?һ??????
#####################epsilon_1?????�?############################
eps_1<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta2+INV(INV(a)-INV(b))
}##epsilon??????��?�????е?һ??????
eps<-function(theta1,theta2,theta3){
  theta1+INV(theta3-theta2)
}
#####################һ????epsilon?�???ȫ????############################
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
}##һ????epsilon??????��?�????е???ȫ????

###############################һ??ģ???????仯�?#####################################
l=100
threshold=-30
thetafirst<-rnorm(kk,0,1)

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




####################################iter_delta2#######################################
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
  thetafirst<-rnorm(kk,0,1)
  
  #############################em?�?#################################################
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
  
  ##############################eps_1?????�?#######################################
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









