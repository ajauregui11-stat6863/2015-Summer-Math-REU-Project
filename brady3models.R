setwd("C:/Users/Adam/Downloads/04to13QBdata")
dat = read.csv("brady.tom01to13.csv", header=TRUE)
dat2014 = read.csv("bradytom2014season.csv", header=TRUE)

wk1 = subset(dat,Week==1)$Fantasy.Points
wk2 = subset(dat,Week==2)$Fantasy.Points
wk3 = subset(dat,Week==3)$Fantasy.Points
wk4 = subset(dat,Week==4)$Fantasy.Points
wk5 = subset(dat,Week==5)$Fantasy.Points
wk6 = subset(dat,Week==6)$Fantasy.Points
wk7 = subset(dat,Week==7)$Fantasy.Points
wk8 = subset(dat,Week==8)$Fantasy.Points
wk9 = subset(dat,Week==9)$Fantasy.Points
wk10 = subset(dat,Week==10)$Fantasy.Points #T.Brady didn't play in Week 10 of 2014
wk11 = subset(dat,Week==11)$Fantasy.Points
wk12 = subset(dat,Week==12)$Fantasy.Points
wk13 = subset(dat,Week==13)$Fantasy.Points
wk14 = subset(dat,Week==14)$Fantasy.Points
wk15 = subset(dat,Week==15)$Fantasy.Points
wk16 = subset(dat,Week==16)$Fantasy.Points
wk17 = subset(dat,Week==17)$Fantasy.Points

#Create histograms for each week: 
par(mfrow=c(3,3))
hist(wk1);hist(wk2);hist(wk3);hist(wk4);hist(wk5);hist(wk6);hist(wk7);hist(wk8);hist(wk9);
hist(wk11);hist(wk12);hist(wk13);hist(wk14);hist(wk15);hist(wk16);hist(wk17)

#check normality
shapiro.test(wk1);shapiro.test(wk2);shapiro.test(wk3);shapiro.test(wk4);
shapiro.test(wk5);shapiro.test(wk6);shapiro.test(wk7);shapiro.test(wk8);
shapiro.test(wk9);shapiro.test(wk10);shapiro.test(wk11);shapiro.test(wk12);
shapiro.test(wk13);shapiro.test(wk14);shapiro.test(wk15);shapiro.test(wk16);
shapiro.test(wk17)

#Combine weekly datasets: 
wkall = list(wk1=wk1, wk2=wk2, wk3=wk3, wk4=wk4, wk5=wk5, wk6=wk6, wk7=wk7, 
             wk8=wk8, wk9=wk9, wk10=wk10, wk11=wk11, wk12=wk12, wk13=wk13, 
             wk14=wk14, wk15=wk15, wk16=wk16, wk17=wk17)
attributes(wkall) = list(names = names(wkall)
       ,row.names=1:max(length(wk1),length(wk2),length(wk3),length(wk4),length(wk5),
        length(wk6),length(wk7),length(wk8),length(wk9),length(wk10),length(wk11),
        length(wk12),length(wk13),length(wk14),length(wk15),length(wk16),length(wk17)), 
       class='data.frame')

########################################################################################
#1. Frequentist's Approach: Assume yi ~ N(mu, sigma)
no.wk = dim(wkall)[2]
means = array(,no.wk) 
sds.2 = array(,no.wk)

for (i in 1:no.wk){
  means[i] = mean(wkall[,i])
  sds.2[i] = sd(wkall[,i])^2
}
sds.2

#Genenerate N observations from N(means, sds.2*(n-1)/(n)), where n is the sample size for each week
N=1000
bradymse1 = array(,no.wk)
bradymad1 = array(,no.wk)

for (j in 1:no.wk){
  if (j !=10){ #skip 10th week 
  n=length(wkall[,j])
ynew = rnorm(N,means[j],sqrt(sds.2[j]*(n-1)/n))  #used mle estimate for sigma^2
bradymse1[j] = sum(   (ynew-subset(dat2014, Week==j)$Fantasy.Points)^2   )/N
bradymad1[j] = sum(  abs(ynew-subset(dat2014, Week==j)$Fantasy.Points)   )/N
}
}

#Print MSE and MAD
bradymse1
bradymad1


###########################################################################################
#2. Bayesian Approach: Assume yi ~ N(mu, sigma^2) and prior mu~N(theta, tau^2)
overall.mean = mean(dat$Fantasy.Points)
overall.sd=sd(dat$Fantasy.Points)

N=1000
bradymse2 = array(,no.wk)
bradymad2 = array(,no.wk)

theta = overall.mean
sigma = overall.sd
tau = overall.sd

for (j in 1:no.wk){
  if (j !=10){ #skip 10th week 
    n=length(wkall[,j])
    mu.post = (tau^2)*mean(wkall[,j])/(tau^2+(sd(wkall[,j]))^2/n) + ((sd(wkall[,j]))^2/n)*theta/(tau^2 + (sd(wkall[,j]))^2/n) 
    rho = (tau^2*sigma^2)/(n*tau^2 + sigma^2)
    ynew = rnorm(N, mu.post, sqrt(sigma^2+rho))
    bradymse2[j] = sum((ynew-subset(dat2014, Week==j)$Fantasy.Points)^2)/N
    bradymad2[j] = sum(abs(ynew-subset(dat2014, Week==j)$Fantasy.Points))/N
  }
}
bradymse2
bradymad2
##################################################################################################
#3. Bayesian Approach: Assume yi ~ N(mu, sigma^2) 
# prior mu~N(theta, tau^2). theta = career mean and tau^2 = career variance # prior sigma^2 ~ gamma(1,1). This is exponential.  

bradymse3 = array(,no.wk)
bradymad3 = array(,no.wk)

for (j in 1:no.wk){
  if (j !=10){ #skip 10th week 
    n=length(wkall[,j])
    n.1=n+1
    wk.dat=list(Y=c(wkall[,j], NA), W=n.1)
    cat("model{
    for(w in 1:W){
    Y[w]~dnorm(mu, 1/sigma2)
    }
    mu~dnorm(17.06, 1/(7.98^2))
    sigma2~dgamma(1,1)
    
    }", file="wkly.jag")
    
    inits<-function(){
      list(mu=rnorm(1, 17.064, 1/(7.98^2)), sigma2=rgamma(1, 1, 1))}
    wk = jags.model(file = "wkly.jag", data=wk.dat, n.chains=1,  n.adapt=1000)
    result= coda.samples(wk, var=c("mu", "sigma2"), n.iter=1000)
    mu.post = result[[1]][,1] #sampled posterior mean
    sigma2.post = result[[1]][,2] #sample posterior variance
    ynew = rnorm(N, mu.post, sqrt(sigma2.post))
    bradymse3[j] = sum((ynew-subset(dat2014, Week==j)$Fantasy.Points)^2)/N
    bradymad3[j] = sum(abs(ynew-subset(dat2014, Week==j)$Fantasy.Points))/N
  }
} #Print MSE and MAD
bradymse3 
bradymad3

###########################################################################################


bradymse = rbind(bradymse1, bradymse2, bradymse3)
bradymad = rbind(bradymad1, bradymad2, bradymad3)
rowMeans(bradymse,na.rm=T)
rowMeans(bradymad,na.rm = T)

par(mfrow=c(2,1))
plot(bradymse[1,], type="o", pch=16, col="black", ylim=c(0,max(bradymse[1,]+10,na.rm=TRUE)), main="Tom Brady", ylab="MSE", xlab="Week")
par(new=T)
plot(bradymse[2,], type="o", pch=16, col="red", ylim=c(0,max(bradymse[1,]+10,na.rm=TRUE)), main="Tom Brady", ylab=" ", xlab=" ")
par(new=T)
plot(bradymse[3,], type="o", pch=16, col="green", ylim=c(0,max(bradymse[1,]+10,na.rm=TRUE)), main="Tom Brady", ylab=" ", xlab=" ")

plot(bradymad[1,], type="o", pch=16, col="black", ylim=c(0,max(bradymad[1,]+2,na.rm=TRUE)), main="Tom Brady", ylab="MAD", xlab="Week")
par(new=T)
plot(bradymad[2,], type="o", pch=16, col="red", ylim=c(0,max(bradymad[1,]+2,na.rm=TRUE)), main="Tom Brady", ylab=" ", xlab=" ")
par(new=T)
plot(bradymad[3,], type="o", pch=16, col="green", ylim=c(0,max(bradymad[1,]+2,na.rm=TRUE)), main="Tom Brady", ylab=" ", xlab=" ")
par(mfrow=c(1,1))  
  