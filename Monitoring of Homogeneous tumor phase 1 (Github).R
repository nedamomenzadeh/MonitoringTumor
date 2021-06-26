rm(list=ls(all.names=TRUE))
#phase1
M=30
Z=matrix(0,nrow = M,ncol = 7)
adc1<-0.002
adc2<-0.001
#clustered tumor
L=read.csv("L.csv",header = T)
L=L[,-1]
library(OpenImageR)
library(qcc)
sigma<-0.01
for (i in 1:M) {
  b0=100
  s01<- exp((-adc1*b0))+rnorm(680,0,sigma)
  s02<- exp((-adc2*b0))+rnorm(151,0,sigma)
  im02<-readImage('image1.png')
  im02<-im02[,,1]
  im02[which(im02<1 & im02>0)]=s01
  im02[which(im02==1)]=s02
  SS1=im02[which(L==2 | L==3)]
  Z[i,1]=mean(SS1)
  
  b1=200
  s1<- exp((-adc1*b1))+rnorm(680,0,sigma)
  s2<- exp((-adc2*b1))+rnorm(151,0,sigma)
  im2<-readImage('image1.png')
  im2<-im2[,,1]
  im2[which(im2<1 & im2>0)]=s1
  im2[which(im2==1)]=s2
  SS2=im2[which(L==2 | L==3)]
  Z[i,2]=mean(SS2)
  
  
  
  
  
  b2=300
  im3=readImage('image1.png')
  s11<- exp((-adc1*b2))+rnorm(680,0,sigma)
  s12<- exp((-adc2*b2))+rnorm(151,0,sigma)
  im3<-im3[,,1]
  im3[which(im3<1 & im3>0)]=s11
  im3[which(im3==1)]=s12
  SS3=im3[which(L==2 | L==3)]
  Z[i,3]=mean(SS3)
  
  
  
  b3=400
  im4=readImage('image1.png')
  s21<- exp((-adc1*b3))+rnorm(680,0,sigma)
  s22<- exp((-adc2*b3))+rnorm(151,0,sigma)
  im4<-im4[,,1]
  im4[which(im4<1 & im4>0)]=s21
  im4[which(im4==1)]=s22
  SS4=im4[which(L==2 | L==3)]
  Z[i,4]=mean(SS4)
  
  
  b4=500
  im5=readImage('image1.png')
  s31<- exp((-adc1*b4))+rnorm(680,0,sigma)
  s32<- exp((-adc2*b4))+rnorm(151,0,sigma)
  im5<-im5[,,1]
  im5[which(im5<1 & im5>0)]=s31
  im5[which(im5==1)]=s32
  SS5=im5[which(L==2 | L==3)]
  Z[i,5]=mean(SS5)
  
  
  b5=600
  im6=readImage('image1.png')
  s41<- exp((-adc1*b5))+rnorm(680,0,sigma)
  s42<- exp((-adc2*b5))+rnorm(151,0,sigma)
  im6<-im6[,,1]
  im6[which(im6<1 & im6>0)]=s41
  im6[which(im6==1)]=s42
  SS6=im6[which(L==2 | L==3)]
  Z[i,6]=mean(SS6)
  
  b6=700
  im7=readImage('image1.png')
  s51<- exp((-adc1*b6))+rnorm(680,0,sigma)
  s52<- exp((-adc2*b6))+rnorm(151,0,sigma)
  im7<-im7[,,1]
  im7[which(im7<1 & im7>0)]=s51
  im7[which(im7==1)]=s52
  SS7=im7[which(L==2 | L==3)]
  Z[i,7]=mean(SS7)
  
}

A=matrix(0,nrow = M,ncol = 2)
B=c(100,200,300,400,500,600,700)
B=B-mean(B)
for (i in 1:M) {
  r=lm(log(Z[i,])~B)
  A[i,]=r$coefficients
}

MSE=NA
#w=ewma(A[,2],lambda = 0.2,nsigmas = 5.5,plot = TRUE)
for (i in 1:M) {
  MSE[i]=(summary(r)$sigma)^2
}
Sxx=sum(B^2)
MSE=mean(MSE)

lcl=mean(A[,2])+(qt(0.025,M*5)*sqrt(((M-1)*MSE)/(M*Sxx)))
ucl=mean(A[,2])+(qt(0.975,M*5)*sqrt(((M-1)*MSE)/(M*Sxx)))


