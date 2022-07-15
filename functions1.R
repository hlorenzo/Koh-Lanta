#' get_design_1_na
#'
#' Function to simulate datasets following **Design 1** structure.
#'
#' @param n integer, number of observations.
#' @param sqrt_1_minus_sig2 real, standard deviation of informative part of model.
#' @param p integer, number of covariates.
#' @param q integer, number of response variables.
#' @param nNA integer, number of response variables.
#'
#' @return List(X,Y) of dataset
#' @export
get_design_1_na <- function(n=50,sqrt_1_minus_sig2=0.99,p=1000,p1=50,q=3,nNA=0){
  # Structure
  out <- get_design_1(n,sqrt_1_minus_sig2,p = p,p1 = p1,q = q)
  if(nNA!=0){
    id_na <- sample(1:(n*p),nNA)
    out$X[id_na] <- NA
  }
  out
}


#' get_design_1_na_MNAR
#'
#' Function to simulate datasets following **Design 1** structure.
#'
#' @param n integer, number of observations.
#' @param sqrt_1_minus_sig2 real, standard deviation of informative part of model.
#' @param p integer, number of covariates.
#' @param q integer, number of response variables.
#' @param nNA integer, number of response variables.
#'
#' @return List(X,Y) of dataset
#' @export
get_design_1_na_MNAR <- function(n=50,sqrt_1_minus_sig2=0.99,p=1000,q=3,
                                 threshold=NA){
  # Structure
  out <- get_design_1(n,sqrt_1_minus_sig2,p,q)
  if(!is.na(threshold)){
    id_na <- which(out$X<as.numeric(threshold))
    out$X[id_na] <- NA
  }
  out
}
