#' Computes the Wassertein distance between a target and a reference
#'
#' @param points input points
#' @param Ref points from the reference,continuous measure
#' @param p 1= W1, 2=W2
#' @return Wasserstein distance
#'
#' @import pdist
#' @export
WassersteinDist<-function(points,Ref, p=2)
{
  if (ncol(points)!=ncol(Ref))
  {
    stop("Dimensions of points and Ref do not agree.")
  }
  if(p==2){
    Mat<-matrix(pdist::pdist(Ref, points)@dist^2, nrow=nrow(points))
  }
  else{
    Mat<-matrix(pdist::pdist(Ref, points)@dist, nrow=nrow(points))
  }
  psi<- WassersteinGoF::SGD_OT(1, Mat)
  #min(psi)
  cost<--WassersteinGoF::Objective(Mat,psi)
  return(cost)
}

#' Computes the Wassertein distance to an isotropic Gaussian
#'
#' @param points input points
#' @param B number of Monte-Carlo samples
#' @param p 1= W1, 2=W2
#' @return Wasserstein distance
#'
#' @import pdist
#' @export
WassersteinDistGaus<-function(points,B, p=2)
{
  dim=ncol(points)
  Ref<-matrix(rnorm(B*dim), ncol=dim)
  if(p==2){
    Mat<-matrix(pdist::pdist(Ref, points)@dist^2, nrow=nrow(points))
  }
  else{
    Mat<-matrix(pdist::pdist(Ref, points)@dist, nrow=nrow(points))
  }
  psi<- WassersteinGoF::SGD_OT(1, Mat)
  #min(psi)
  cost<--WassersteinGoF::Objective(Mat,psi)
  return(cost)
}

#' Generate the quantiles for an isotropic Gaussian.
#'
#' @param quantiles a list of quantiles you want to return. By default, the function returns only a sample vector from the distribution under the null.
#' @param nRep takes the number of replications
#' @param samplesize sample size to be used
#' @param p 1=W1, 2=W2
#' @param B the number of points to approximate the continuous distribution
#' @param dim the dimension
#' @param nCores number of cores for the parellelisation. By default, the number is the number of cores of your machine minus 1.
#' @return list with sample vector from the distribution under the null and specified quantiles (if any)
#'
#' @import pdist
#' @import foreach
#' @import doParallel
#' @export
generateGausWassersteinQuantile<-function(quantiles=1.2, nRep,samplesize, p=2,B, dim, nCores=0)
{
  if( nCores ==0){nCores<-parallel::detectCores()-1} # I never use all cores to ensure that other processes can run.
  cl <-parallel::makeCluster(nCores[1])
  parallel::clusterCall(cl=cl,function() library(WassersteinGoF))
  doParallel::registerDoParallel(cl)
  store<-foreach(i = 1:nRep,.combine='c') %dopar% {
    WassersteinDistGaus(points=matrix(rnorm(samplesize*dim),ncol=dim),B=B,p=p)
  }
  parallel::stopCluster(cl)
  if(quantiles[1]==1.2)
  {return(store)}
  else
  {
    return(list("obs"=store, "quantiles"=quantile(store,quantiles)))
  }
}
#' Generate the quantiles for the Gaussian family
#'
#' @param quantiles a list of quantiles you want to return. By default, the function returns only a sample vector from the distribution under the null.
#' @param nRep takes the number of replications
#' @param samplesize sample size to be used
#' @param p 1=W1, 2=W2
#' @param B the number of points to approximate the continuous distribution
#' @param dim the dimension
#' @param nCores number of cores for the parellelisation. By default, the number is the number of cores of your machine minus 1.
#' @return list with sample vector from the distribution under the null and specified quantiles (if any)
#'
#' @import pdist
#' @import foreach
#' @import doParallel
#' @export
generateGausFamWassersteinQuantile<-function(quantiles=1.2, nRep,samplesize, p=2,B, dim, nCores=0 )
{

  if(nCores ==0){nCores<-parallel::detectCores()-1} # I never use all cores to ensure that other processes can run.
  cl <-parallel::makeCluster(nCores[1])
  parallel::clusterCall(cl=cl,function() library(WassersteinGoF))
  doParallel::registerDoParallel(cl)
  store<-foreach(i = 1:nRep,.combine='c') %dopar% {
    WassersteinDistGaus(points=rGausResid(samplesize,dim),B=B,p=p)
  }
  parallel::stopCluster(cl)
  if(quantiles[1]==1.2)
  {return(store)}
  else
  {
    return(list("obs"=store, "quantiles"=quantile(store,quantiles)))
  }
}

