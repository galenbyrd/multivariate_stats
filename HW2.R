# HOMEWORK 2
# Galen Byrd

# CHAPTER 4
#10a
mu<-c(3,1,4)
SIGMA <-rbind(c(6,1,-2),c(1,13,4),c(-2,4,4))
c<-c(2,-1,3)
zmu<-c%*%mu
zSIGMA<-t(c)%*%SIGMA%*%c
# z~(17,21)
#10b
cb<-rbind(c(1,1,1),c(1,-1,2))
zbmu<-cb%*%mu
zbSIGMA<-cb%*%SIGMA%*%t(cb)
# z1,z2~N2((8,10),[(29,-1),(-1,9)])
#10c
# y2~N(1,13)
#10d
# y1,y3~N2((3,4),[(6,-2),(-2,4)])
#10e
C<-rbind(c(1,0,0),c(0,0,1),c(.5,.5,0))
cmu<-C%*%mu
cSIGMA<-C%*%SIGMA%*%t(C)
# y1,y3,.5(y1+y2)~N3((3,4,2),[(6,-2,3.5),(-2,4,1),(3.5,1,5.25)])
  
#11b
SIGMARoot <- eigen(SIGMA)$vectors%*%diag(sqrt(eigen(SIGMA)$values))%*%t(eigen(SIGMA)$vectors)
SIGMARootInverse<-solve(SIGMARoot)
# z=SIGMARootInverse%*%(y-3,y-1,y-4)
#11c
# Distributed chi-square(3) because statistical distance is distributed chi-square(p)

#14a
mu<-c(2,-3,4)
SIGMA<-rbind(c(4,-3,0),c(-3,6,0),c(0,0,5))
# Dependent due to covariance not equal to 0
#14b
# Independent due to covariance equal to 0
#14c
# Independent due to covariance equal to 0
#14d
c<-rbind(c(1,1,0),c(0,0,1))
dSIGMA<-c%*%SIGMA%*%t(c)
# Independent due to covariance equal to 0
#14e
c<-rbind(c(1,0,1),c(0,1,0))
eSIGMA<-c%*%SIGMA%*%t(c)
# Dependent due to covariance not equal to 0

#16a
mu<-c(2,-1,3,1)
SIGMA<-rbind(c(7,3,-3,2),c(3,6,0,4),c(-3,0,5,-2),c(2,4,-2,4))
(n<-mu[1:2])
(z<-SIGMA[1:2,3:4]%*%solve(SIGMA[3:4,3:4]))
# n+z%*%c(x1-3,x2-1)
#16b
SIGMA[1:2,1:2]-SIGMA[1:2,3:4]%*%solve(SIGMA[3:4,3:4])%*%SIGMA[3:4,1:2]


# CHAPTER 5
#11
nullvec<-c(6,11)
Y<-cbind(c(3,6,5,10),c(10,12,14,9))
meanvec<-colMeans(Y)
varmat<-var(Y)
n=nrow(Y)
p=ncol(Y)
T2 <- n*t(meanvec-nullvec)%*%solve(varmat)%*%(meanvec-nullvec)
(critval <- (n-1)*p/(n-p)*qf(.95,p,n-p))
(T2>critval)
# Fail to reject H0

#12
probe <- read.table("T3_6_PROBE.dat")
probe<-probe[,-1]
#a
nullvec<-c(30,25,40,25,30)
meanvec <- colMeans(probe)
varmat <- var(probe)
n=nrow(probe)
p=ncol(probe)
(T2 <- n*t(meanvec-nullvec)%*%solve(varmat)%*%(meanvec-nullvec))
(critval <- (n-1)*p/(n-p)*qf(.95,p,n-p))
(T2>critval)
# reject H0

#b
t.test(probe$V2,mu=nullvec[1])
# reject H0, true mean not equal to 30
t.test(probe$V3,mu=nullvec[2])
# Fail to reject H0, true mean equal to 25
t.test(probe$V4,mu=nullvec[3])
# reject H0, true mean not equal to 40
t.test(probe$V5,mu=nullvec[4])
# Fail to reject H0, true mean equal to 25
t.test(probe$V6,mu=nullvec[5])
# Fail to reject H0, true mean equal to 30


fBeetles <- read.table("T5_5_FBEETLES.dat")
fBeetles<-fBeetles[,-1]
#16a
g1<- subset(fBeetles,V2==1)[,-1]
g2<- subset(fBeetles,V2==2)[,-1]
p <- ncol(g1)
n1 <- nrow(g1)
n2 <- nrow(g2)
mean1 <- colMeans(g1)
mean2 <- colMeans(g2)
S1 <- var(g1)
S2 <- var(g2)
Sg <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
T2 <- n1*n2/(n1+n2)*t(mean1-mean2)%*%solve(Sg)%*%(mean1-mean2)
a <- p*(n1+n2-2)/(n1+n2-p-1)
crit.val <- a*qf(.95,p,n1+n2-p-1)
T2>crit.val
# Reject H0

#16b
#####Calculating a t-test for each variable
p.vals <- rep(1,p)

#Test on y1
t.stat <- (mean1[1]-mean2[1])/sqrt((1/n1+1/n2)*Sg[1,1])
p.vals[1] <- 2*pt(abs(t.stat), n1+n2-2, lower.tail = F)
#Test on y2
t.stat <- (mean1[2]-mean2[2])/sqrt((1/n1+1/n2)*Sg[2,2])
p.vals[2] <- 2*pt(abs(t.stat), n1+n2-2, lower.tail = F)
#Test on y3
t.stat <- (mean1[3]-mean2[3])/sqrt((1/n1+1/n2)*Sg[3,3])
p.vals[3] <- 2*pt(abs(t.stat), n1+n2-2, lower.tail = F)
#Test on y4
t.stat <- (mean1[4]-mean2[4])/sqrt((1/n1+1/n2)*Sg[4,4])
p.vals[4] <- 2*pt(abs(t.stat), n1+n2-2, lower.tail = F)
# Reject H0 for all

#c
solve(Sg)%*%(mean1-mean2)

#20a
goods <- read.table("T5_8_GOODS.dat")
goods<-goods[,-1]

g1<- subset(goods,V2==1)[,-1]
g2<- subset(goods,V2==2)[,-1]
p <- ncol(g1)
n1 <- nrow(g1)
n2 <- nrow(g2)
mean1 <- colMeans(g1)
mean2 <- colMeans(g2)
S1 <- var(g1)
S2 <- var(g2)
Sg <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
T2 <- n1*n2/(n1+n2)*t(mean1-mean2)%*%solve(Sg)%*%(mean1-mean2)
a <- p*(n1+n2-2)/(n1+n2-p-1)
crit.val <- a*qf(.95,p,n1+n2-p-1)
T2>crit.val
# Reject H0

#20b
solve(Sg)%*%(mean1-mean2)


essay <- read.table("T5_9_ESSAY.dat")
essay <- essay[,-1]
#21a
dword<-as.matrix(essay[,1]-essay[,3])
dverb<-as.matrix(essay[,2]-essay[,4])
diff <-cbind(dword,dverb)
p<-ncol(diff)
n<-nrow(diff)
meandiff <- colMeans(diff)
varmat <- var(diff)
(T2 <- n*t(meandiff)%*%solve(varmat)%*%meandiff)
(critval <- p*(n-1)/(n-p)*qf(.95,p,n-p))
T2>critval
#reject H0

#21b
solve(varmat)%*%meandiff

#21c
t.test(dword)
# Reject H0
t.test(dverb)
# Reject H0

