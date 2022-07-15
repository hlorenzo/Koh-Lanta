## y de dim 1
## theta de dimension 3
## x de dim 1
## eta de dimension 3
# n <- 300
# x <- t(t(rnorm(n)))
# y <- 10 + 5*x + t(t(rnorm(n)))
x <- scale(1:10)
y <- x + rnorm(10,sd = 0.1)
n <- length(x)

nNa <- 1
naX <- 1
x[naX] <- NA
# y[naY] <- NA

## A priori non informatif sur les coeff de régression de y par x
prec <- 1e3
thetas <- list(
  a0 = seq(-3,3,length.out = prec),
  a  = seq(-3,3,length.out = prec),
  s2 = seq(0,5,length.out = prec)
)
poste_thetas <- list(
  a0 = dnorm(thetas$a0,mean = 0,sd = 1),
  a  = dnorm(thetas$a,mean = 0,sd = 1),
  s2 = dchisq(thetas$s2,df = 2)
)


M <- 1000

iter <- 1
test <- T

x_cur <- x
y_cur <- y

while(test){
  a0_out = a_out = s2_out <- rep(0,prec)
  for(m in 1:M){
    cat(m);cat(" ")
    a0 <- sample(thetas$a0,1,prob = poste_thetas$a0)
    a <- sample(thetas$a,1,prob = poste_thetas$a)
    s2 <- sample(thetas$s2,1,prob = poste_thetas$s2)
    for(i in 1:n){
      if(i %in% naX){
        mux_m <- a0+a*y_cur[i]
        varx_m <- s2
        x_cur[i] <- rnorm(1,mux_m,sqrt(varx_m))
      }
    }
    momo <- lm(x~y,data = data.frame(list(x=x_cur,y=y_cur)))
    a0_out <- a0_out + momo$coefficients[1]/M
    a_out <- a_out + momo$coefficients[2]/M
    s2_out <- s2_out + mean((momo$residuals)^2)/M
  }
  ## Draw to get theta
  id_samp_theta <- sample(1:length(newDistriSigm),size = M,replace = T,prob = newDistriSigm)
  density_rho <- unlist(lapply(id_samp_theta,function(ii){
    covo <- Sigs_list[[ii]]
    covo[1,2]/sqrt(covo[1,1]*covo[2,2])
  }))

  density_sigma <- newDistriSigm
  iter <- iter + 1
  if(iter==12) test <- F
  plot(density(density_rho),main=iter)
}

iter <- iter + 1



