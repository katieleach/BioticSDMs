## ## CAR model fitted with INLA ## ##

## Load packages
require(spdep)
require(RANN)
require(INLA)

## Read in data
data <- read.table("sample.txt", header=T)

## Function from Beale et al. 2014
sortBase <- function(vec, n.knots = 2) {
  ## Function to calculate bases for regression splines. Modified from code
  ## provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
  ## penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1–24(2005).
  ## Parameter vec is a vector defining the raw data vector, n.knots defines the
  ## number of knots in the GAM.
  N              <- length(vec)
  x.time         <- c(vec)
  ZFE            <- cbind(rep(1,N), x.time)
  x.knots        <- quantile(unique(x.time), seq(0, 1, 
                    length = (n.knots+2))[-c(1, (n.knots+2))], na.rm = TRUE)
  Z_K            <- (abs(outer(x.time,x.knots,"-")))^3
  OMEGA.all      <- (abs(outer(x.knots,x.knots,"-")))^3
  svd.OMEGA.all  <- svd(OMEGA.all)
  sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u)*
                    sqrt(svd.OMEGA.all$d)))
  Z.out          <- t(solve(sqrt.OMEGA.all, t(Z_K)))
  return(Z.out)
}

## Define regression splines
Z.var1 <- matrix(c(scale(sortBase(data[,4]))), ncol = 2)
Z.var2 <- matrix(c(scale(sortBase(data[,5]))), ncol = 2)

## Calculate distance between cells for CAR component
coords <- SpatialPoints(as.matrix(data[,1:2]))
k1 <- knn2nb(knearneigh(coords, longlat=TRUE))
all.linked <- max(unlist(nbdists(k1, coords)))
neb <- dnearneigh(coords, 0, all.linked)
neb2 <- nb2mat(neb)
neb2[neb2>0] <- 1

## Prepare model
y=1000 #Number of grid cells
node_CAR <- rep(1:y, times=1, length.out=y, each=1)
z=10 #Shape parameter
hyperpar_CAR <- list(param=c(z, 0.001))
formula_CAR_mod1 <- spp1 ~ Z.var1[,1] + Z.var1[,2] + Z.var2[,1] + Z.var2[,2] +
                    f(node_CAR, model="besag", graph=neb2, 
                    hyper=list(theta=hyperpar_CAR))

## INLA model
model_CAR_mod1   <- inla(formula=formula_CAR_mod1, data=data, family="binomial", 
                    control.compute=list(dic=TRUE, waic=TRUE), 
                    control.fixed=list(prec.intercept=0.001, prec=0.001), 
                    verbose=FALSE)
summary(model_CAR_mod1)

## Recombining regression splines from Beale et al. 2014
inv.logit.new     <- function(x) {p <- 1 / (1+exp(-x)); p}
var1Mean <- inv.logit.new(outer(Z.var1[,1],  model_CAR_mod1$summary.fixed[2,1], "*") + 
            matrix(scale(model_CAR_mod1$summary.fixed[1,1], scale = F), 
            nrow=y, ncol=1, byrow=T) + 
            outer(Z.var1[,2], model_CAR_mod1$summary.fixed[3,1], "*"))
var2Mean <- inv.logit.new(outer(Z.var2[,1],  model_CAR_mod1$summary.fixed[4,1], "*") + 
            matrix(scale(model_CAR_mod1$summary.fixed[1,1], scale = F), 
            nrow=y, ncol=1, byrow=T) + 
            outer(Z.var2[,2], model_CAR_mod1$summary.fixed[5,1], "*"))
