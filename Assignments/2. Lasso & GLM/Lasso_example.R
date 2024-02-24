# prostate.Lasso.R Implementation of Lasso estimation for 
#          linear multiple regression 
#
# Pedro Delicado, april 2016
#

SoftThr <- function(x,lambda){sign(x)*max(0,abs(x)-lambda)}

# Prostate data for illustrating Lasso regression 
# 
# Data downloaded from
# http://statweb.stanford.edu/~tibs/ElemStatLearn/
# 10-04-2016
#
# Goal: To examine the correlation between the level of 
# log of prostate-specific antigen (lpsa)  
# and a number of clinical measures 
# in 97 men who were about to receive a radical prostatectomy.
#
# 
prostate <- read.table("prostate_data.txt", header=TRUE, row.names = 1)
#plot(prostate)
train.sample <- which(prostate$train==TRUE)
val.sample <- which(prostate$train==FALSE)

n.train <- length(train.sample)
n.val <- length(val.sample)

# centering the data and scaling the predcitors
mean.training <- apply(prostate[train.sample,],2,mean)
sd.training <- apply(prostate[train.sample,],2,sd)
  
Y <- scale( prostate$lpsa[train.sample], center=TRUE, scale=FALSE)
X <- scale( as.matrix(prostate[train.sample,1:8]), center=TRUE, scale=TRUE)

n <- dim(X)[1]
p <- dim(X)[2]

m.Y.tr <- mean.training[9]
m.X.tr <- mean.training[1:8]
s.Y.tr <- sd.training[9]
s.X.tr <- sd.training[1:8]
  
Y.val <- prostate$lpsa[val.sample]-m.Y.tr
X.val <- (as.matrix(prostate[val.sample,1:8])- matrix(m.X.tr,nrow=n.val,ncol=p,byrow = TRUE))%*%diag(1/s.X.tr)

# Prostate data: Lasso estimation for a given lambda
# 
# lambda <- 0 # OLS estimation
# lambda <- max(apply(X,2,function(x,y){abs(mean(x*y))},y=Y))
lambda <- .5*max(apply(X,2,function(x,y){abs(mean(x*y))},y=Y))
max.iter <- 50
eps.beta <- 1e-15
eps.r <- 1e-15

beta.0 <- rep(0,length(X[1,]))
beta.1 <- beta.0 
stop.rule <- FALSE
r.0 <- Y - X%*%beta.0 
beta.00 <- beta.0
r.00 <- r.0
n.iter <- 0
while(!stop.rule){
  n.iter <- n.iter + 1
  for (j in 1:p){
    beta.1[j] <- SoftThr(beta.0[j] + mean(X[,j]*r.0),lambda)
    r.1 <- r.0 - (beta.1[j]-beta.0[j])*X[,j]
    # Updating beta.0[j] and r.0
    beta.0[j] <- beta.1[j]
    r.0 <- r.1
  }
  if ((n.iter==max.iter) | 
      (mean((beta.00-beta.0)^2)<=eps.beta) | 
      (abs(1-mean(r.0^2)/mean(r.00^2))<=eps.r)) stop.rule <- TRUE
  beta.00 <- beta.0
  r.00 <- r.0
}
beta.lasso.lambda <- beta.0
# MSPE in the validation sample
Y.hat.val <- X.val %*% beta.0
plot(Y.hat.val,Y.val,asp=1)
abline(a=0,b=1,col=2)

SE.val <- (Y.val-Y.hat.val)^2
MSPE.val.lambda <- mean(SE.val)
# With no model: prediction of Y.val=(0,...,0)
# having a MSPE equal to mean(Y.val^2).
# Proportion at which the MSPE has been reduced:
MSPE.val.lambda/mean(Y.val^2)


# Prostate data: Lasso estimation and 
#   coefficient path.
# 
# Lambda.v vector
lambda.max <- max(apply(X,2,function(x,y){abs(mean(x*y))},y=Y))
n.lambdas <- 100
lambda.v <- seq(lambda.max,0,length=n.lambdas)

beta.path <- matrix(0,nrow=n.lambdas, ncol=p)
names(beta.path) <- names(prostate[1:8])
MSPE.val <- numeric(n.lambdas)
sd.SE.val <- numeric(n.lambdas)

max.iter <- 10
eps.beta <- 1e-15
eps.r <- 1e-15

Y.hat.val <- X.val %*% beta.path[1,]
SE.val <- (Y.val-Y.hat.val)^2
MSPE.val[1] <- mean(SE.val)
sd.SE.val[1] <- sd(SE.val)/sqrt(n.val)

for (l in 2:n.lambdas){
  beta.0 <- beta.path[l-1,]
  beta.1 <- beta.0 
  stop.rule <- FALSE
  lambda <- lambda.v[l]
  r.0 <- Y - X%*%beta.0 
  beta.00 <- beta.0
  r.00 <- r.0
  n.iter <- 0
  while(!stop.rule){
    n.iter <- n.iter + 1
    for (j in 1:p){
      beta.1[j] <- SoftThr(beta.0[j] + mean(X[,j]*r.0),lambda)
      r.1 <- r.0 - (beta.1[j]-beta.0[j])*X[,j]
      # Updating beta.0[j] and r.0
      beta.0[j] <- beta.1[j]
      r.0 <- r.1
    }
    if ((n.iter==max.iter) | 
        (mean((beta.00-beta.0)^2)<=eps.beta) | 
        (abs(1-mean(r.0^2)/mean(r.00^2))<=eps.r)) stop.rule <- TRUE
    beta.00 <- beta.0
    r.00 <- r.0
  }
  beta.path[l,] <- beta.0
# MSPE in the validation sample
  Y.hat.val <- X.val %*% beta.0
  SE.val <- (Y.val-Y.hat.val)^2
  MSPE.val[l] <- mean(SE.val)
  sd.SE.val[l] <- sd(SE.val)/sqrt(n.val)
}

# shrinkage factor s
s <- apply(beta.path,1,function(beta){sum(abs(beta))})
s <- s/max(s)

l.min <- which.min(MSPE.val)
lambda.min <- lambda.v[l.min]
s.min <- s[l.min]

print(cbind(names(prostate)[1:8],beta.path[l.min,]))

op <- par(mfrow=c(2,1))
# MSPE.val 
plot(lambda.v,MSPE.val, xlab="lambda",
     ylim=c(min(MSPE.val-sd.SE.val),max(MSPE.val+sd.SE.val)))
segments(lambda.v,MSPE.val-sd.SE.val, 
         lambda.v,MSPE.val+sd.SE.val, 
         col=8,lend=2)

abline(v=lambda.min,col=2,lty=2)

# coefficients path
plot(range(lambda.v), range(beta.path),type="n",
     xlab="lambda",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(lambda.v,beta.path[,j],col=4)
  #  points(lambda.v,beta.path[,j],pch=19,cex=.7,col=4)
}
text(0*(1:p), beta.path[n.lambdas,],names(prostate)[1:p],pos=4)
abline(v=lambda.min,col=2,lty=2)
par(op)

op <- par(mfrow=c(2,1))
# MSPE.val against the log(lambda)
lg.lambda.v <- log(lambda.v[-n.lambdas])
plot(lg.lambda.v,MSPE.val[-n.lambdas], xlab="log lambda", 
     ylim=c(min(MSPE.val-sd.SE.val),max(MSPE.val+sd.SE.val)))
segments(lg.lambda.v,MSPE.val[-n.lambdas]-sd.SE.val[-n.lambdas], 
         lg.lambda.v,MSPE.val[-n.lambdas]+sd.SE.val[-n.lambdas], 
         col=8,lend=2)

abline(v=log(lambda.min),col=2,lty=2)

# coefficients path aginst the log(lambda)
plot(range(lg.lambda.v), range(beta.path),type="n",
     xlab="log lambda",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(lg.lambda.v,beta.path[-n.lambdas,j],col=4)
  #  points(lambda.v,beta.path[,j],pch=19,cex=.7,col=4)
}
text(0*(1:p)+min(lg.lambda.v), beta.path[n.lambdas,],names(prostate)[1:p],pos=4)
abline(v=log(lambda.min),col=2,lty=2)
par(op)

op <- par(mfrow=c(2,1))
# MSPE.val aginst the shrinkage factor s
plot(s,MSPE.val, xlab="s, shrinkage factor", 
     ylim=c(min(MSPE.val-sd.SE.val), max(MSPE.val+sd.SE.val)))
segments(s,MSPE.val-sd.SE.val, 
         s,MSPE.val+sd.SE.val, 
         col=8,lend=2)

abline(v=s.min,col=2,lty=2)

# coefficients path aginst the shrinkage factor s
plot(c(0,1.1), range(beta.path),type="n",
     xlab="s, shrinkage factor",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(s,beta.path[,j],col=4)
  #  points(lambda.v,beta.path[,j],pch=19,cex=.7,col=4)
}
text(1+0*(1:p), beta.path[n.lambdas,],names(prostate)[1:p],pos=4)
abline(v=s.min,col=2,lty=2)

par(op)

# Prostate data: Lasso with glmnet
#
# Comparing with glmnet results
require(glmnet)

lasso.1 <- glmnet(X,Y, standardize=FALSE, intercept=FALSE)
cv.lasso.1 <- cv.glmnet(X,Y, standardize=FALSE, intercept=FALSE,nfolds=n)

op <- par(mfrow=c(2,1))
plot(cv.lasso.1)
plot(lasso.1,xvar="lambda")
abline(v=log(cv.lasso.1$lambda.min),col=2,lty=2)
abline(v=log(cv.lasso.1$lambda.1se),col=2,lty=2)
print(coef(lasso.1,s=cv.lasso.1$lambda.min))
print(coef(lasso.1,s=cv.lasso.1$lambda.1se))
par(op)

# To scale or not to scale?
#
# results of glmnet without previous centering and scaling

Y.no.scl <- prostate$lpsa[train.sample]
X.no.scl <- as.matrix(prostate[train.sample,1:8])

lasso.2 <- glmnet(X.no.scl,Y.no.scl, standardize=TRUE, intercept=TRUE)
cv.lasso.2 <- cv.glmnet(X.no.scl,Y.no.scl, standardize=TRUE, intercept=TRUE,nfolds=n)
print(coef(lasso.2,s=cv.lasso.2$lambda.1se))
mean.Y <- mean(Y.no.scl)
mean.X <- apply(prostate[train.sample,1:8],2,mean)
sd.X <- apply(prostate[train.sample,1:8],2,sd)

# intercept:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[1]
# b) from the fitted model centering and scaling
mean.Y - sum((mean.X/sd.X) * coef(lasso.1,s=cv.lasso.2$lambda.1se)[-1])

# coefficients:
# a) No centering, no scaling
coef(lasso.2,s=cv.lasso.2$lambda.1se)[-1]
# b) from the fitted model centering and scaling
coef(lasso.1,s=cv.lasso.2$lambda.1se)[-1] / sd.X

# prediction in the validation set:
# a) No centering, no scaling
Y.val.no.scl <- prostate$lpsa[val.sample]
Y.val.no.scl.hat <- predict(lasso.2,newx=as.matrix(prostate[val.sample,1:8]),s=cv.lasso.2$lambda.1se)
plot(Y.val.no.scl.hat,Y.val.no.scl,asp=1)
abline(a=0,b=1,col=2)
# b) from the fitted model centering and scaling
Y.val.hat <- predict(lasso.1,newx=X.val,s=cv.lasso.1$lambda.1se)
Y.val.no.scl.hat.1 <- mean.Y + Y.val.hat 
plot(Y.val.no.scl.hat.1,Y.val.no.scl,asp=1)
abline(a=0,b=1,col=2)
plot(Y.val.no.scl.hat,Y.val.no.scl.hat.1,asp=1)
abline(a=0,b=1,col=2)
