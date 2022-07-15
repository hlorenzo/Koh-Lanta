## y de dim 1
## theta de dimension 3
## x de dim 1
## eta de dimension 3
# n <- 300
# x <- t(t(rnorm(n)))
# y <- 10 + 5*x + t(t(rnorm(n)))
x <- t(t(c(1,1,-1,-1,2,2,-2,-2,NA,NA,NA,NA)))
y <- t(t(c(1,-1,1,-1,NA,NA,NA,NA,2,2,-2,-2)))
n <- length(x)

nNa <- 4
naX <- which(is.na(x))#sample(1:n,size = nNa)
naY <- which(is.na(y))#sample(1:n,size = nNa)
# x[naX] <- NA
# y[naY] <- NA

## A priori non informatif sur les coeff de régression de y par x
prec <- 2^8
GO <- c(-2,2)
valsPoss <- do.call(cbind,lapply(1:8,function(i){rep(c(rep(GO[1],2^(i-1)),rep(GO[2],2^(i-1))),2^(8-i))}))


Sigs_list <- lapply(1:prec,function(i){
  x_cur <- x
  y_cur <- y
  x_cur[naX] <- valsPoss[i,1:4]#sample(c(-2,2),nNa,replace = T)#rnorm(nNa,mean(na.omit(x)),sd(na.omit(x)))
  y_cur[naY] <- valsPoss[i,1:4+4]#rnorm(nNa,mean(na.omit(y)),sd(na.omit(y)))
  var(cbind(x_cur,y_cur))
})
density_rho <- unlist(lapply(1:prec,function(ii){
  covo <- Sigs_list[[ii]]
  covo[1,2]/sqrt(covo[1,1]*covo[2,2])
}))
density_sigma <- rep(1,prec)/prec

M <- 1000

iter <- 1
test <- T

x_cur <- x
y_cur <- y
x_cur[naX] <- rnorm(nNa,mean(na.omit(x)),sd(na.omit(x)))
y_cur[naY] <- rnorm(nNa,mean(na.omit(y)),sd(na.omit(y)))

density_rho0 <- density_rho

while(test){
  newDistriSigm <- rep(0,length(density_sigma))
  for(m in 1:M){
    cat(m);cat(" ")
    sig_id <- sample(1:prec,1,prob = density_sigma)
    Sig_m <- Sigs_list[[sig_id]]
    sig12_m <- Sig_m[1,1]
    sig22_m <- Sig_m[2,2]
    rho_m <- Sig_m[1,2]/sqrt(sig12_m*sig22_m)
    for(i in 1:n){
      if(i %in% naX){
        mux_m <- rho_m*sqrt(sig12_m/sig22_m)*y_cur[i]
        varx_m <- sig12_m*(1-rho_m^2)
        x_cur[i] <- rnorm(1,mux_m,sqrt(varx_m))
      }
      if(i %in% naY){
        muy_m <- rho_m*sqrt(sig22_m/sig12_m)*x_cur[i]
        vary_m <- sig22_m*(1-rho_m^2)
        y_cur[i] <- rnorm(1,muy_m,sqrt(vary_m))
      }
    }
    matou_m <- var(cbind(x_cur,y_cur))
    for(j in 1:length(density_sigma)){
      newDistriSigm[j] <-  newDistriSigm[j] + diwish(W=matou_m,n,S=Sigs_list[[j]])/M
    }
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



