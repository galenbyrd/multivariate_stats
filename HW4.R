library(readxl)
track <- read_excel("track.xlsx")

track3<-cbind(track[,1],100/track[,2],200/track[,3],400/track[,4],800/track[,5]*60,
              1500/(track[,6]*60),5000/(track[,7]*60),10000/(track[,8]*60),42195/(track[,9]*60))

# 1 #################################################
# A
#Getting the sample size and number of variables
n <- nrow(track)
p <- ncol(track)-1

###Finding the eigenvalue and vectors
R <- cov(track[,-1]) #covariance matrix
e.vec <- eigen(R)$vectors
(e.val <- eigen(R)$values)

# B
mean(e.val)

plot(1:p,e.val, xlab="Number",ylab="Eigenvalue",main="Scree Plot for Track", type="l")
points(1:p,e.val)

percentage <- rep(0,p)
for (i in 1:p){
  percentage[i] <- sum(e.val[1:i])/sum(e.val)   
}
(percentage)
# Use just first as it accounts for 98% of variability

# C
(e.vec)
# First one is -.97 marathon, -.18 10,000m and -.11 400m. We can call this long distance running
# Second one is not super interpretable, it is .82 400m, .35 200m, .29 10,000m, .21 100m.

# D
z<-as.matrix(track[,-1])%*%e.vec[,1]
order(z)
# Top 5: usa,australia,japan,portugal,netherla
# Bottom 5: cookis,wsamoa,singapor,domrep,malaysia

# 2 #################################################
# A
R <- cor(track[,-1]) #correlation matrix
e.vec <- eigen(R)$vectors
(e.val <- eigen(R)$values)

# B
mean(e.val)

plot(1:p,e.val, xlab="Number",ylab="Eigenvalue",main="Scree Plot for Track", type="l")
points(1:p,e.val)

percentage <- rep(0,p)
for (i in 1:p){
  percentage[i] <- sum(e.val[1:i])/sum(e.val)   
}
(percentage)
# Use just first two as they account for 93% of variability and the second eigenvalue is >=avg(e.val)

# C
(e.vec)
# First one is .3+ for all variables. We can call this general running ability
# Second one is not super interpretable, it is .56 100m, .46 200m, .43 marathon, .31 5,000m. We can call this short/long distance running

# D
track.sc <- scale(track[,-1],center=T,scale=apply(track[,-1],2,sd))

z<-as.matrix(track.sc)%*%e.vec[,1:2]
order(z[,1])
# Top 5: usa,gbni,italy,ussr,gdr
# Bottom 5: cookis,wsamoa,maritiu,png,singapor

# E
# The two results are different. They both agree that usa/cookis and wsamoa are the
# best/worst respectively, but after that they give different answers.

# 3 #################################################
# A
track3<-cbind(track[,1],100/track[,2],200/track[,3],400/track[,4],800/(track[,5]*60),
              1500/(track[,6]*60),5000/(track[,7]*60),10000/(track[,8]*60),42195/(track[,9]*60))

# B
R <- cov(track3[,-1]) #covariance matrix
e.vec <- eigen(R)$vectors
(e.val <- eigen(R)$values)

# C
(mean(e.val))

plot(1:p,e.val, xlab="Number",ylab="Eigenvalue",main="Scree Plot for Track", type="l")
points(1:p,e.val)

percentage <- rep(0,p)
for (i in 1:p){
  percentage[i] <- sum(e.val[1:i])/sum(e.val)   
}
(percentage)
# Use first two as they account for 94% of variability

# D
(e.vec)
# First one is .3+ for all variables. We can call this general running
# Second one is not super interpretable, it is .60 100m, .47 200m, .42 marathon, .29 5,000/10,000m. We can call this short/long distance running

# E
track.sc <- scale(track3[,-1],center=T,scale=apply(track3[,-1],2,sd))
z<-as.matrix(track.sc)%*%e.vec
order(z[,1])
# Top 5: usa,gbni,italy,ussr,gdr
# Bottom 5: cookis,wsamoa,mauritiu,png,singapor

# F
# I prefer the method from __________________________

# 4 #################################################
# A
track<-cbind(track[,1],track[,2:4],track[,5:9]*60)
track <- read_excel("track.xlsx")

R <- cov(track[,-1])
meanvec<-colMeans(track[,-1])
d2 <- NULL
for (i in 2:9){
  d2 <- c(d2,(as.matrix(track[i,-1])-meanvec)%*%solve(R)%*%(t(as.matrix(track[i,-1])-meanvec)))
}
j=2:9
chisqpoints <- qchisq((j-1/2)/n,df=p)
qqplot(chisqpoints,d2, main="Chisq Q-Q Plot")


# B
# Data is MVN because Chi-Square Q-Q plot is linear so use Maximum likelihood method

# C
library(psych)
R <- cor(track[,-1])
FA.ML1 <- factanal(covmat=R, factors=1, rotation="none")
Psi.ml1 <- diag(diag(R-FA.ML1$loadings%*%t(FA.ML1$loadings)))
(FA.ML1.res <- round(R-(FA.ML1$loadings%*%t(FA.ML1$loadings)+Psi.ml1),2))
FA.ML1$loadings

FA.ML2 <- factanal(covmat=R, factors=2, rotation="none")
Psi.ml2 <- diag(diag(R-FA.ML2$loadings%*%t(FA.ML2$loadings)))
(FA.ML2.res <-round(R-(FA.ML2$loadings%*%t(FA.ML2$loadings)+Psi.ml2),2))
FA.ML2$loadings

# We should use m=2

# D
FA.ML2 <- factanal(covmat=R, factors=2, rotation="varimax")
Psi.ml2 <- diag(diag(R-FA.ML2$loadings%*%t(FA.ML2$loadings)))
FA.ML2.res <-round(R-(FA.ML2$loadings%*%t(FA.ML2$loadings)+Psi.ml2),2)
(L <-FA.ML2$loadings)
FA.ML2$scores
# Factor 1 is long distance running ability
# Factor 2 is short distance running ability

# E
#D <- sqrt(solve(diag(diag(cov(track[,-1])))))
v<-track[53,-1]

(f <- t(as.matrix(L))%*%solve(R)%*%t(as.matrix(v-meanvec)))


