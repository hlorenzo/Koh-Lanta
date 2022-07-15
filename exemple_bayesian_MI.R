doImputationDensity <- function(x,y,maxIter=100,errorMin=1e-4)
{
  test <- T
  iter <- 1
  q <- ncol(y)
  p <- ncol(x)
  sigmays <- NULL
  id_na_x <- which(is.na(x))
  id_na_y <- which(is.na(y))
  id_na_vect <- which(is.na(x),arr.ind = T)
  N_na_x <- length(id_na_x)
  N_na_y <- length(id_na_y)

  ## Initialisation du modèle de régression de y
  theta0 <- list(b=t(t(rnorm(p+1))),var=var(y,na.rm = T))
  errs <- NULL
  y_hypo <- y
  x_hypo <- x
  x_hypo[id_na_x] <- colMeans(x,na.rm = T)[id_na_vect[,2]]

  while(test){
    ## y, equation (7), Imputation step for y
    for(id in 1:N_na_y){
      i <- id_na_y[id]
      x_i <- t(x_hypo[i,,drop=F])
      mu_i <- crossprod(theta0$b,rbind(1,x_i))
      # distri_y[,id] <- rnorm(n = 1,mean = mu_i,sd = sqrt(theta0$var))
      y_hypo[id_na_y[id]] <- rnorm(n = 1,mean = mu_i,sd = sqrt(theta0$var))
    }

    ## xi paramètre de régression , equation (5), Posterior step for x parameter
    y_m <- y
    y_m[id_na_y] <- y_hypo[id_na_y]
    y_m <- cbind(1,y_m)
    b_x <- solve(crossprod(y_m))%*%crossprod(y_m,x_hypo)
    residuals <- x_hypo-y_m%*%b_x
    sigma_x <- apply(residuals,2,var)
    xi <- list(b=b_x,var=sigma_x)

    ## x, equation (6)
    # distri_x <- matrix(NA,1,N_na_x)
    for(id in 1:N_na_x){
      i_j <- id_na_vect[id,]
      i <- i_j[1]
      j <- i_j[2]
      y_i <- t(y_hypo[i,,drop=F])
      mu_i_j <- t(xi$b[,j,drop=F])%*%rbind(1,y_i)
      x_hypo[id_na_x[id]] <- rnorm(1,mu_i_j,sd = sqrt(xi$var[j]))#distri_x[,id]
    }
    # x_hypo[id_na_x] <- colMeans(distri_x)

    ## theta, equation (4)
    x_m <- x
    d_m <- x_hypo[id_na_x]#distri_x[m,]
    for(id in 1:N_na_x){
      i_j <- id_na_vect[id,]
      i <- i_j[1]
      j <- i_j[2]
      x_m[i,j] <- d_m[id]
    }
    x_m <- cbind(1,x_m)
    b_y <- solve(crossprod(x_m))%*%crossprod(x_m,y_hypo)
    residuals <- y_hypo-x_m%*%b_y
    sigma_y <- var(residuals)

    theta <- list(b=b_y,var=sigma_y)
    sigmays <- c(sigmays,sigma_y)
    # ide0 <- ide
    theta0 <- theta

    if(iter==maxIter) test <- F
    iter <- iter + 1
  }
  out <- list(x=x_hypo,y=y_hypo,
              b=list(x_knowd_y =b_x,y_knowd_x =b_y),
              sigma2=list(x=sigma_x,y=sigmays))
  out
}



errMin <- 1e-6
sdX <- 3
p <- 10
n <- 50
M <- 100
x = x_init <- matrix(rnorm(n*p),n,p)
sdY_th <- 5
y = y_init <- x[,1,drop=F]*20 + 10 + rnorm(n,sd = sdY_th)
p2 <- 3
xp <- x[,p]
x[,(p-p2+1):(p)] <- matrix(rep(xp,p2),nrow = n,byrow = F)#+matrix(rnorm(n*p2,sd = 0.5),n,p2)
z_init <- cbind(y_init,x_init)

prop <- 0.1
n_na_x <- ceiling(n*p*prop)
n_na_y <- ceiling(n*prop)
id_i_na_x <- sample(1:(n*p),n_na_x)
id_i_na_y <- sample(1:(n),n_na_y)
x[id_i_na_x] <- NA
y[id_i_na_y] <- NA
x_miss <- rnorm(n_na_x)
x_hypo <- x
x_hypo[id_i_na_x] <- x_miss
y_miss <- rnorm(n_na_y)
y_hypo <- y
y_hypo[id_i_na_y] <- y_miss



res <- doImputationDensity(x,y,maxIter=50,errorMin=1e-4)
plot(sqrt(res$sigma2$y),ylim=c(0,max(c(sqrt(res$sigma2$y),sdY_th))));abline(h=sdY_th,col=2)
sigmas_y <- rep(NA,M)
for(m in 1:M)
{
  idm <- sample(1:n,size = n,replace = T)
  resm <- doImputationDensity(x[idm,,drop=F],y[idm,,drop=F],maxIter=50,errorMin=1e-4)
  points(sqrt(resm$sigma2$y),ylim=c(0,max(c(sqrt(resm$sigma2$y),sdY_th))),col=m,cex=0.2,pch=16)
  sigmas_y[m] <- sqrt(resm$sigma2$y)[length(sqrt(resm$sigma2$y))]
}
hist(sigmas_y,20)

