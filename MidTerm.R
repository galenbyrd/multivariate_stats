cork <- read.table("cork.txt", header=T)
cork <- cork[,-1]
mean.vec<- colMeans(cork)

# 1 B
C<-rbind(c(1,-1,0,0),c(1,0,-1,0),c(0,1,0,-1))
mean.vec<-C%*%mean.vec

# 1 C
null.vec <- c(0,0,0)
n <- nrow(cork)
p <- ncol(cork)
varmat <- var(cork)
(varmat2 <- C%*%varmat%*%t(C))
T2 <- n*t(mean.vec-null.vec)%*%solve(varmat2)%*%(mean.vec-null.vec)

#Calculating the critical value
crit.val <- p*(n-1)/(n-p)*qf(.95,p,n-p)
# reject H0 when bigger than critical value
# so reject H0


#Finding the p-value
T2.p <- (n-p)/((n-1)*p)*T2
(p.val <- pf(T2.p, p, n-p, lower.tail=F))
# reject H0 at 5% significance level

# 1 D
p.vals <- rep(1,3)

p.vals[1] <- t.test(cork$N-cork$E,mu=null.vec[1])$p.value
p.vals[2] <- t.test(cork$N-cork$S,mu=null.vec[2])$p.value
p.vals[3] <- t.test(cork$E-cork$W,mu=null.vec[3])$p.value

# Fail to reject H0 for first variable (N-E)
.05/3



# 2 ################################################################################################
sparrow <- read.table("sparrow.txt", header = T)
sparrow <- sparrow[,-1]

sparrow[,6]<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

SA <-sparrow[1:21,]
SD <-sparrow[22:49,]
colMeans(SA)
colMeans(SD)

m1 <- manova(cbind(totlen,alarext,beaklen,humlen,sternlen)~as.factor(V6), data=sparrow)
summary(m1, test="Wilks")


# 2 D
m11 <- aov(totlen~as.factor(V6),data=sparrow)
summary(m11)
#there is the biggest difference among groups

m12 <- aov(alarext~as.factor(V6),data=sparrow)
summary(m12)

m13 <- aov(beaklen~as.factor(V6),data=sparrow)
summary(m13)

m14 <- aov(humlen~as.factor(V6),data=sparrow)
summary(m14)

m15 <- aov(sternlen~as.factor(V6),data=sparrow)
summary(m15)



# 3 ################################################################################################
us <- readxl::read_xlsx("USStates.xlsx")

m1 <- manova(cbind(Smokers,PhysicalActivity,Obese,College,NonWhite)~as.factor(Region), data=us)
summary(m1, test="Wilks")
# Reject H0

summary(m1,test="Wilks")$SS

#Calculating the E and H matrix
H <- summary(m1,test="Wilks")$SS$`as.factor(Region)`
E <- summary(m1,test="Wilks")$SS$Residuals

#Calculating Wilks' test statistic
wilks <- det(E)/det(E+H)

# 3 D
summary(aov(Smokers~Region,data=us))
summary(aov(PhysicalActivity~Region,data=us))
summary(aov(Obese~Region,data=us))
summary(aov(College~Region,data=us))
summary(aov(NonWhite~Region,data=us))

# Small p values mean significant.


# 4 ################################################################################################
#Calculating the basic statistics
N <- nrow(us)
k <- length(unique(us$Region))
p <- 5

W <- subset(us, Region=="W", select=c(Smokers,PhysicalActivity,Obese,College,NonWhite))
n1 <- nrow(W)
mean1 <- colMeans(W)
S1 <- var(W)

S <- subset(us, Region=="S", select=c(Smokers,PhysicalActivity,Obese,College,NonWhite))
n2 <- nrow(S)
mean2 <- colMeans(S)
S2 <- var(S)

MW <- subset(us, Region=="MW", select=c(Smokers,PhysicalActivity,Obese,College,NonWhite))
n3 <- nrow(MW)
mean3 <- colMeans(MW)
S3 <- var(MW)

NE <- subset(us, Region=="NE", select=c(Smokers,PhysicalActivity,Obese,College,NonWhite))
n4 <- nrow(NE)
mean4 <- colMeans(NE)
S4 <- var(NE)


n <- c(n1,n2,n3,n4)
#Calculating the pooled covariance matrix assuming the covariance across countries are the same
Spl <- ((n1-1)*S1+(n2-1)*S2+(n3-1)*S3+(n4-1)*S4)/(N-k)

u <- det(Spl)/prod(diag(Spl))

u.1 <- -((N-1)-1/6*(2*p+5))*log(u)

pchisq(u.1,.5*p*(p-1),lower.tail = F)
# Reject H0, the variables are NOT independent
