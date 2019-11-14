### 1 ##############################################################
library(readxl)
firms <- read_excel("firms.xlsx")
# A
#Separating data
bankrupt <- subset(firms, Status=="B")[,-5]
sound <- subset(firms, Status=="S")[,-5]
#Calculating the basic statistics
p <- ncol(bankrupt)
n1 <- nrow(bankrupt)
n2 <- nrow(sound)
#Calculating the mean vectors and covariance matrices
mean.b <- colMeans(bankrupt)
mean.s <- colMeans(sound)
S.b <- var(bankrupt)
S.s <- var(sound)
S.pl <- ((n1-1)*S.b+(n2-1)*S.s)/(n1+n2-2)
#Calculating Hotelling's T2
T2 <- n1*n2/(n1+n2)*t(mean.b-mean.s)%*%solve(S.pl)%*%(mean.b-mean.s)
#Calculating the critical value
a <- p*(n1+n2-2)/(n1+n2-p-1)
crit.val <- a*qf(.95,p,n1+n2-p-1)
T2>crit.val
p.val <- 1-pf(1/a*T2,p,n1+n2-p-1)
# Reject null. There is evidence (p-val=1.36e-05) to suggest a difference among bankrupt/financially sound banks
# in at least one of the variables x1,x2,x3,x4.

# B
# We perform a hypothesis test to make sure the groups are different enough
# in at least one of the given variables.

# C
#Finding the E and H matrix using MANOVA
m1 <- manova(cbind(x1,x2,x3,x4)~as.factor(Status),data=firms)
H <- summary(m1)$SS[[1]]
E <- summary(m1)$SS[[2]]
#Calculating the eigenvalues and vectors for the discriminant analysis
e.vals <- Re(round(eigen(solve(E)%*%H)$values,digits=4))
e.vecs <- Re(round(eigen(solve(E)%*%H)$vectors,digits=4))
a1 <- e.vecs[,1]
# So the variables that contribute to the separation from most to least are: x2,x4,x3,x1

# D
t(a1)%*%mean.b
t(a1)%*%mean.s
(zc <- .5*t(a1)%*%(mean.b+mean.s))
# for new data point, if its mean is >-.26 it's Bankrupt, <-.26 it's Financially stable

# E
library(MASS)
k <- 2
LDA <- lda(Status~x1+x2+x3+x4 , data=firms, prior=rep(1,k)/k)
(error <- mean(firms$Status != predict(LDA)$class) )
#Apparent error rate = .087
LDA.CV <- lda(Status~x1+x2+x3+x4 , data=firms, prior=rep(1,k)/k, CV=T)
(error <- mean(firms$Status != LDA.CV$class) )
#Error rate using cross validation = .109

# Apparent error rate underestimates actual error rate. Using cross validation we remove each observation 
# individially and recalculate the classification rules, which should pull the apparent error rate towards
# the actual error rate, which it does.

### 2 ##############################################################
cereal <- read_excel("cereal.xlsx")
# A
D <- dist(cereal[,-1],diag=T, upper=T)
m.sl <-hclust(d = D, method="single")
plot(as.dendrogram(m.sl), main="Dendrogram for Single Linkage")
rect.hclust(m.sl,k=2,border="red")

m.cl <-hclust(d = D, method="complete") 
plot(as.dendrogram(m.cl), main="Dendrogram for Complete Linkage")
rect.hclust(m.cl,k=2,border="red")

m.al <-hclust(d = D, method="average") 
plot(as.dendrogram(m.al), main="Dendrogram for Average Linkage")
rect.hclust(m.al,k=2,border="red")

# B
# I prefer the average method because it does not stretch/shrink the data like single and complete linkage do

# C
# I would use 2 clusters

# D
# Based on that answer, Product/Total are in a group, and the rest are in one big group.
# This makes sense as Product/Total both have a lot more vitamins than the others.

# E
c1<-rbind(cereal[7,-1],cereal[9,-1])
c2<-rbind(cereal[1:6,-1],cereal[8,-1],cereal[10:12,-1])

#Calculating the basic statistics
p <- ncol(c1)
n1 <- nrow(c1)
n2 <- nrow(c2)
#Calculating the mean vectors and covariance matrices
mean1 <- colMeans(c1)
mean2 <- colMeans(c2)
S.pl <- var(cereal[,-1])
S1 <- var(c1)
S2 <- var(c2)
S.pl <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
#Calculating Hotelling's T2
T2 <- n1*n2/(n1+n2)*t(mean1-mean2)%*%solve(S.pl)%*%(mean1-mean2)
#Calculating the critical value
a <- p*(n1+n2-2)/(n1+n2-p-1)
crit.val <- a*qf(.95,p,n1+n2-p-1)
T2>crit.val
(p.val <- 1-pf(1/a*T2,p,n1+n2-p-1))
# Reject null. There is evidence (p-val=5.84e-05) to suggest a difference in means among cereal groups
# meaning we should not combine the clusters.

### 3 ##############################################################
house.med <- read.table("housdat.txt",header=T)
house <- house.med[,-14]
# A
n <- nrow(house)
p <- ncol(house)
diag(cov(house))
# Use correlation matrix because the variances in the data are not similar.
# Some have small variances and others have huge variances, which would effect the PC's.
# The variance for FTAX is the one that is too large.

# B
R <- cor(house)
e.vec <- eigen(R)$vectors
e.val <- eigen(R)$values

plot(1:p,e.val, xlab="Number",ylab="Eigenvalue",main="Scree Plot for House", type="l")
points(1:p,e.val)

percentage <- rep(0,p)
for (i in 1:p){
  percentage[i] <- sum(e.val[1:i])/sum(e.val)   
}
(percentage)

# I would use the first 5 PC's so that we retain 80% of the variability in the data.
# This is also where the scree plot really levels off.

# C
# 80%

# D
z<-as.matrix(house)%*%e.vec
plot(z[,1],z[,2])
# Yes there are a few groups 

# E
plot(house.med[,14],z[,1])
plot(house.med[,14],z[,2])
# It does not appear these two PC's are very useful in predicting median home prices. There are big groups horizontally,
# which means that for one value of our PC, there is a big range of potential house prices.

### 4 ##############################################################
house <- house[,-4]
p <- 12
# A
# y = mu + L %*% f + error
# Where L is a p by m matrix of loadings and f is the vector of common factors
# Assume: E(f)=0, E(error)=0, var(f)=I, var(error)=Psi, and f and error are independent.

# B
library(psych)
R <- cor(house)

FA.ML4 <- factanal(covmat=R, factors=4, rotation="varimax")
Psi.ml4 <- diag(diag(R-FA.ML4$loadings%*%t(FA.ML4$loadings)))

FA.ML5 <- factanal(covmat=R, factors=5, rotation="varimax")
Psi.ml5 <- diag(diag(R-FA.ML5$loadings%*%t(FA.ML5$loadings)))

FA.ML6 <- factanal(covmat=R, factors=6, rotation="varimax")
Psi.ml6 <- diag(diag(R-FA.ML6$loadings%*%t(FA.ML6$loadings)))

FA.ML7 <- factanal(covmat=R, factors=7, rotation="varimax")
Psi.ml7 <- diag(diag(R-FA.ML7$loadings%*%t(FA.ML7$loadings)))

#Loglikelihood
m=4
FA4 <- (FA.ML4$loadings%*%t(FA.ML4$loadings)+Psi.ml4)
ll4 <- -n/2*(log(det(FA4))+sum(diag(solve(FA4)%*%R)))
AIC4 <- -2*ll4+2*(p*(m+1)-m*(m-1))

m=5
FA5 <- (FA.ML5$loadings%*%t(FA.ML5$loadings)+Psi.ml5)
ll5 <- -n/2*(log(det(FA5))+sum(diag(solve(FA5)%*%R)))
AIC5 <- -2*ll5+2*(p*(m+1)-m*(m-1))

m=6
FA6 <- (FA.ML6$loadings%*%t(FA.ML6$loadings)+Psi.ml6)
ll6 <- -n/2*(log(det(FA6))+sum(diag(solve(FA6)%*%R)))
AIC6 <- -2*ll6+2*(p*(m+1)-m*(m-1))

m=7
FA7 <- (FA.ML7$loadings%*%t(FA.ML7$loadings)+Psi.ml7)
ll7 <- -n/2*(log(det(FA7))+sum(diag(solve(FA7)%*%R)))
AIC7 <- -2*ll7+2*(p*(m+1)-m*(m-1))

# We should use 7 factors as this minimizes the value of AIC

e.val <- eigen(R)$values
percentage <- rep(0,p)
for (i in 1:p){
  percentage[i] <- sum(e.val[1:i])/sum(e.val)   
}
(percentage)

# We retain 90% of the variability in the data

# C
FA.ML7$loadings
# considering loading >.6 significant.
# Factor 1: CRIM,INDEX,FTAX. Adding crime rates with closeness to highway and property tax rate
# Factor 2: PLAND,WDIS. Adding closeness to jobs and average size of lots.
# Factor 3: ARM,LSP. Contrasting average rooms with % lower status.
# Factor 4: PAGE
# Factor 5: PTR
# Factor 6: Trivial factor
# Factor 7: Trivial Factor

