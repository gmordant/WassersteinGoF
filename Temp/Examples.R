#--------------------------------------------------------------------------------------------
#            Example: Distance computation between a point cloud and a Reference measure
#--------------------------------------------------------------------------------------------

points<-matrix(rnorm(200),ncol=2)
Ref<-matrix(rt(20000,6),ncol=2) # Ref must contain way more observations than the samplesize

WassersteinDist(points,Ref,p=2)



#--------------------------------------------------------------------------------------------
#            Example: Distance computation between a point cloud and a Gaussian reference
#--------------------------------------------------------------------------------------------

points<-matrix(rnorm(200),ncol=2)

WassersteinDistGaus(points, 100000, p=2)

#--------------------------------------------------------------------------------------------
#            Example: Generate quantiles for isotropic Gaussian reference
#--------------------------------------------------------------------------------------------


library(pdist)
library(doParallel)

begin<-Sys.time()
u<-generateGausWassersteinQuantile(nRep=90,samplesize=200,dim=2, p=1,B=100000, nCores=4)
Sys.time()-begin

plot(ecdf(u))


#--------------------------------------------------------------------------------------------
#            Example: Produce Gaussian residuals
#--------------------------------------------------------------------------------------------
# Residuals are empirically centred and (Cholesky-) rescaled Gaussian observations
# Produce residuals of a Gaussian sample of size 200 in dimension 3
x<-rGausResid(200,3)
plot(x[,1:2])


#--------------------------------------------------------------------------------------------
#            Example: Generate quantiles for Gaussian family
#--------------------------------------------------------------------------------------------

library(pdist)
library(doParallel)

begin<-Sys.time()
u<-generateGausFamWassersteinQuantile(quantiles=c(0.9,0.95),nRep=30,samplesize=200,dim=2, p=1,B=100000, nCores=4)
Sys.time()-begin

plot(ecdf(u$obs))


#--------------------------------------------------------------------------------------------
#            Example: Hypothesis test for the Gaussian family
#--------------------------------------------------------------------------------------------

# Set p
p=1

# Generate a sample
points<-matrix(rt(400,7),ncol=2)
points<-scale(points, scale = FALSE) # center
points<-points%*%solve(chol(cov(points))) # rescale

# Generate the data to produce quantiles
begin<-Sys.time()
u<-generateGausFamWassersteinQuantile(nRep=1000,samplesize=200,dim=2, p=p,B=200000, nCores=4)
Sys.time()-begin
#---------------------Results----------------------
#        90%              95%         99%
#   0.2139332         0.2179838    0.2257982

# Compute the test statistic

Tm<-WassersteinDistGaus(points, B=200000,p)

# Compute the p-value
pVal= mean(Tm>u)
