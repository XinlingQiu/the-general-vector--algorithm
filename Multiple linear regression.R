library("MASS")
library("abind")
seed=sample(100000,1)
set.seed(seed)
N=100
p=10
kk=p+1
X<-matrix(0,N,kk)
X[,1]=1
S <- toeplitz((p:1)/p)
R <- rWishart(1, p, S)
R=R[,,1]
miu=matrix(0,p,1)
X[,2:(p+1)]=mvrnorm(N,miu,R)
inv=solve(t(X)%*%X)
Y=X%*%rnorm(kk,0,1)
def=sample(1:N,N/2)
Y[def]=NA


INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
##################### EM ############################
em<-function(theta){
  theta=matrix(theta,kk,1)
  #E step
  for(i in def){
    Y[i]=X[i,]%*%theta
  }
  ##M step
  theta=inv%*%t(X)%*%Y
  theta=matrix(theta,1,kk)
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
############################# this is about tabel 1 ##############################
l=100
threshold=-30
thetafirst<-rnorm(kk,0,1)


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

start=20
l=20
cl=ceiling(l/2)
theta_eps=array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_eps[,,1]=0
theta_eps[,,2]=theta[start:(start+l-1),]
diff<-matrix(0,l,cl)
for(m in 3:(l+1)){
  for(j in 1:(l-m+2)){
    theta_eps[j,,m]=eps(theta_eps[j+1,,m-2],theta_eps[j,,m-1],theta_eps[j+1,,m-1])
  }
}
for(j in 1:cl){
  for(i in 1:(l-2*(j-1))){
    diff[i,j]=metrics(theta_eps[i,,2*j],truepara)
  }
}

View(diff)
#################################### tabel 1 end #########################################

########################## this is for tabel 2 ##################################
nn=1000
threshold=-10
iternum<-matrix(0,nn,3)##
for(i in 1:nn){
  thetafirst<-rnorm(kk,0,1)
  ########################### em #################################################
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
  
  ############################## delta2 #######################################
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
  ############################ general epsilon #######################################
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
        result_2=theta[l,,2]
      }else{
        theta[L,,k]=eps(theta[L+1,,k-2],theta[L,,k-1],theta[L+1,,k-1])
        if(k%%2==0&&(L>1)){
          diff_2=metrics(theta[L,,k],theta[L-1,,k])
          if(is.finite(diff_2)&&diff_2<diff){diff=diff_2;result_2=theta[L,,k]}
        }
      }
    }
  }
  iternum[i,3]=l
}
summary(iternum)

#################################### tabel 2 end #########################################







