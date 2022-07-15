#' get_toy_example
#'
#' Function to simulate datasets following **Toy Example** structure.
#'
#' @param n integer, number of observations.
#' @param sqrt_1_minus_sig2 real, standard deviation of informative part of model.
#' @param p integer, number of covariates.
#'
#' @return List(X,Y) of dataset
#' @export
get_toy_example <- function(n=50,sqrt_1_minus_sig2=0.95,p=1000){
  # Structure
  repX <- 50
  A <- sqrt_1_minus_sig2*matrix(c(rep(1,repX),rep(0,p-repX)),nrow = 1,byrow = T)
  D <- sqrt_1_minus_sig2*matrix(1,nrow = 1,byrow = T)
  # Observations
  d <- ncol(A)+nrow(A)+ncol(D)
  psi <- MASS::mvrnorm(n = n,mu = rep(0,d),Sigma = diag(d))
  phi <- psi[,1:nrow(A),drop=F]
  epsilonX_info <- psi[,nrow(A)+1:repX]*sqrt(1-sqrt_1_minus_sig2^2)
  epsilonX_noise <- psi[,nrow(A)+repX+1:(ncol(A)-repX)]
  epsilonY_info <- psi[,nrow(A)+ncol(A)+1,drop=F]*sqrt(1-sqrt_1_minus_sig2^2)
  # X and Y
  X <- phi%*%A + cbind(epsilonX_info,epsilonX_noise)
  Y <- phi%*%D + epsilonY_info
  list(X=X,Y=Y)
}

#' get_design_1
#'
#' Function to simulate datasets following **Design 1** structure.
#'
#' @param n integer, number of observations.
#' @param sqrt_1_minus_sig2 real, standard deviation of informative part of model.
#' @param p integer, number of covariates.
#' @param q integer, number of response variables.
#'
#' @return List(X,Y) of dataset
#' @export
get_design_1 <- function(n=50,sqrt_1_minus_sig2=0.99,p=1000,p1=50,q=3){
  # Structure
  alpha3 <- 1/sqrt(3)
  alpha2 <- 1/sqrt(2)
  repX <- p1
  A1 <- c(rep(alpha3,repX),rep(0,p-repX))
  A2 <- c(rep(0,repX),rep(alpha2,repX),rep(0,p-2*repX))
  A <- matrix(c(rep(A1,3),rep(A2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2
  D1 <- c(rep(alpha3,1),rep(0,q-1))
  D2 <- c(rep(0,1),rep(alpha2,1),rep(0,q-2))
  D <- matrix(c(rep(D1,3),rep(D2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2
  # Observations
  d <- ncol(A)+nrow(A)+ncol(D)
  psi <- MASS::mvrnorm(n = n,mu = rep(0,d),Sigma = diag(d))
  phi <- psi[,1:nrow(A)]
  epsilonX_info <- psi[,nrow(A)+1:(2*repX)]*sqrt(1-sqrt_1_minus_sig2^2)
  epsilonX_noise <- psi[,nrow(A)+(2*repX)+1:(ncol(A)-2*repX)]
  epsilonY_info <- psi[,nrow(A)+ncol(A)+1:2]*sqrt(1-sqrt_1_minus_sig2^2)
  epsilonY_noise <- psi[,d]
  # X and Y
  X <- phi%*%A + cbind(epsilonX_info,epsilonX_noise)
  Y <- phi%*%D + cbind(epsilonY_info,epsilonY_noise)
  list(X=X,Y=Y)
}

#' get_design_2
#'
#' Function to simulate datasets following **Design 2** structure.
#'
#' @param seed real, the seed to be used.
#' @param n integer, number of observations.
#' @param q integer, number of response variables.
#' @param p1 integer, number of covariates in **X1**.
#' @param p2 integer, number of covariates in **X2**.
#' @param sigma1 real, standard deviation of additive noise in **X1**.
#' @param sigma2 real, standard deviation of additive noise in **X2**.
#' @param sigmaY real, standard deviation of additive noise in **Y**.
#' @param ncpX integer, number of latent variables in **X1** and **X2**.
#' @param ncpXCom integer, number of latent variables in common between **X1** and **X2**.
#' @param ncpXYCom integer, number of latent variables in common between (**X1**/**X2**) and **Y**.
#' @param plot boolean, whether or not to plot the dataset.
#'
#' @return list(Xs,Y,S) Xs is a list of 2 matrices, Y is a mtrix, S is the list of the two spectra **S1** and **S2** and
#'  the list of the common spectra in **S1** and **S2**.
#' @export
get_design_2 <- function(seed=1,n=50,q=5,p1=500,p2=5000,
                         sigma1=0.05,sigma2=0.05,sigmaY=0.1,
                         ncpX=10,ncpXCom=5,ncpXYCom=3,plot=F){
  # Internal Function
  simulateData <- function(n, ncp, p, sigma, sigmaNoise=0.1, ConcVarFact=0.8, n_min_peaks=3){
    meanPeakSigma <- sigma
    sigPeakSigma <- sigma / 4
    axis <- 1:p
    S <- matrix(0,ncp,p)
    C <- matrix(0,n,ncp)
    for(i in 1:ncp){
      npeaks <- 3+ceiling(10*runif(1))
      peakheights <- runif(npeaks)
      sigmas <- runif(npeaks) * sigPeakSigma + meanPeakSigma
      position <- runif(npeaks) * p
      for(j in 1:npeaks){
        S[i,] <- S[i,] + peakheights[j] * exp(-0.5 * ((axis - position[j]) / sigmas[j])^2)
      }
    }
    meanC <- sort(10^runif(ncp),decreasing = T)
    varC <- ConcVarFact * meanC * runif(ncp)
    for(i in 1:ncp){
      C[,i] <- rnorm(n = n,mean = meanC[i], sd = varC[i]/2)
    }
    X <- C%*%S;X <- X/max(abs(X))
    E <- matrix(rnorm(n*p,sd = sigmaNoise),nrow = n,ncol = p)
    X <- X*sqrt(1-sigmaNoise^2) + E
    list(X=X, C=C, S=S, E=E)
  }

  set.seed(seed)
  # ncpX for each X separately
  Data_1 <- simulateData(n=n, ncp=ncpX, p=p1, sigma=20, sigmaNoise = sigma1, ConcVarFact=0.8, n_min_peaks=5)
  Data_2 <- simulateData(n=n, ncp=ncpX, p=p2, sigma=30, sigmaNoise = sigma2, ConcVarFact=0.8, n_min_peaks=5)
  S1 <- Data_1$S;C1 <- Data_1$C;X1 <- Data_1$X
  S2 <- Data_2$S;C2 <- Data_2$C
  # ncpXCom in common
  C2[,1:ncpXCom] <- C1[,1:ncpXCom,drop=F]
  X2 <- C2%*%S2;X2 <- X2/max(abs(X2))
  E2 <- matrix(rnorm(n*p2,sd = sigma2),nrow = n,ncol = p2)
  X2 <- X2*sqrt(1-sigma2^2) + E2
  # Build Y on ncpXYCom components
  Y <- scale(Data_1$C[,1:ncpXYCom,drop=F])
  # Add extra variables and noise
  E_y <- matrix(rnorm(n*q,sd = sigmaY),nrow = n,ncol = q)
  if(q>ncpXYCom){
    Y <- cbind(Y,matrix(rnorm(n*(q-ncpXYCom)),nrow = n))*sqrt(1-sigmaY^2) + E_y
  }else{
    Y <- Y*sqrt(1-sigmaY^2) + E_y
  }
  if(plot){
    layout(matrix(c(1,2,2,3,3,3),nrow = 2,byrow = T))
    matplot(t(X1),lty=1,type="l")
    matplot(t(X2),lty=1,type="l")
    corXY <- cor(cbind(X1,X2),Y)
    matplot(corXY,type="l")
  }
  list(Xs=list(X1=X1,X2=X2),Y=Y,S=list(S1=S1,S2=S2,SXY=list(S1=S1[1:ncpXYCom,],S2=S2[1:ncpXYCom,])))
}
