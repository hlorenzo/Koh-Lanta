prodNA <- function(A,b){
  n <- nrow(A)
  if(!any(is.na(A))){
    out <- A%*%b
  }else{
    out <- rep(0,n)
    for(j in 1:nrow(A)){
      out[j] <- sum(na.omit(A[j,]*b))
    }
    out <- t(t(out))
  }
  out
}

plpls <- function(x,y,errMin=1e-3){
  p <- ncol(x)
  w0 <- t(t(rep(0,p)))
  u <- rowMeans(y,na.rm = T)
  test <- T
  it <- 0
  while(test){
    w <- prodNA(t(x),u)
    w <- w/sqrt(sum(w^2))
    t <- prodNA(x,w)
    c <- prodNA(t(y),t)/sum(t^2)
    u <- prodNA(y,c)
    err <- sum((w-w0)^2)
    w0 <- w
    if(err<errMin) test <- F
    it <- it + 1
  }
  p <- prodNA(t(x),t)/sum(t^2)
  if(which.max(abs(p))!=which.max((p))){
    w <- -w
    p <- -p
    t <- -t
    c <- -c
    u <- -u
  }
  list(w=w,t=t,c=c,u=u,p=p)
}

pls_NA <- function(x,y,R,errMin=1e-3,doCV=F,B=NULL){
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  t = matrix(NA,n,R)
  U = matrix(NA,n,R)
  P = matrix(NA,p,R)
  W = matrix(NA,p,R)
  C = matrix(NA,q,R)
  if(doCV){
    idB <- sample(c(unlist(lapply(1:B,function(b){rep(b,(n-n%%B)/B)})),1:(n%%B)),n,replace=F)
    MSES <- rep(NA,R)
    for(r in 1:R){
      mse_r <- 0
      for(b in 1:B){
        idb <- which(idB==b)
        x_train <- x[-idb,,drop=F]
        y_train <- y[-idb,,drop=F]
        x_test <- x[idb,,drop=F]
        y_test <- y[idb,,drop=F]
        mo <- pls_NA(x_train,y_train,r,errMin=errMin)
        x_test_imp <- apply(x_test,2,function(vv){
          idNA <- which(is.na(vv))
          if(length(idNA)>0){
            vv[idNA] <- mo$mux[idNA]
          }
          vv
        })
        y_est <- predict_pls_NA(mo,x_test_imp)
        mse_r <- mse_r + sum((y_est-y_test)^2)
      }
      MSES[r] <- mse_r/(q*B)
    }
    R <- which.min(MSES)
  }
  mux <- colMeans(x,na.rm = T)
  muy <- colMeans(y,na.rm = T)
  na_x <- which(is.na(mux))
  na_y <- which(is.na(muy))
  if(length(na_x)>0){
    mux[na_x] <- 0
  }
  if(length(na_y)>0){
    muy[na_y] <- 0
  }
  xi <- scale(x,scale = F)
  yi <- scale(y,scale = F)
  for(r in 1:R){
    momo <- plpls(xi,yi,errMin=errMin)
    t[,r] <- momo$t
    U[,r] <- momo$u
    P[,r] <- momo$p
    W[,r] <- momo$w
    C[,r] <- momo$c
    xi <- xi - tcrossprod(t[,r],P[,r])
    yi <- yi - tcrossprod(t[,r],C[,r])
  }
  B <- tcrossprod(W%*%solve(crossprod(P,W)),C)
  list(t=t,W=W,U=U,P=P,C=C,B=B,mux=mux,muy=muy,R=R)
}

predict_pls_NA <- function(mo,x_test){
  n_test <- nrow(x_test)
  MUX <- matrix(rep(mo$mux,n_test),nrow = n_test,byrow = T)
  MUY <- matrix(rep(mo$muy,n_test),nrow = n_test,byrow = T)
  Y_est <- MUY + (x_test-MUX)%*%mo$B
}


MI_PLS <- function(x,y,R,M=5,proper=T,errMin=1e-9){
  datas_imputed0 <- list()
  models = Rs <- list()
  stat_simu <- matrix(NA,M,2)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  idNA <- which(is.na(x),arr.ind = T)
  for(m in 1:M){
    ## Hyper-parameter and model estimations
    test <- T
    while(test){
      id_m <- sample(1:n,n,replace = T)
      if(length(c(1:n)[-id_m])>0) test <- F
    }
    idNA_m <- which(is.na(x[id_m,]))
    idNA_m_id <- which(is.na(x[id_m,]),arr.ind = T)
    it <- 0
    mo_m <- pls_NA(x[id_m,],y[id_m,,drop=F],R = R,doCV=F,B=B)
    mu_X <- matrix(rep(mo_m$mux,n),nrow = n,byrow = T)
    mu_Y <- matrix(rep(mo_m$muy,n),nrow = n,byrow = T)
    x_id_m <- mo_m$t%*%t(mo_m$P)+mu_X
    W <- mo_m$W
    t <- apply(W,2,function(w){prodNA(x,w)})
    P <- apply(t,2,function(tt){prodNA(t(x),tt)/sum(tt^2)})
    C <- apply(t,2,function(tt){prodNA(t(y),tt)/sum(tt^2)})
    B <- W%*%solve(t(P)%*%W)%*%t(C)
    B_R <- lapply(1:R,function(r){W[,1:r,drop=F]%*%solve(t(P[,1:r,drop=F])%*%W[,1:r,drop=F])%*%t(C[,1:r,drop=F])})

    Esti <- x
    Esti_0 <- t%*%t(P)+mu_X
    if(nrow(idNA)>0){
      for(k in 1:nrow(idNA)){
        Esti[idNA[k,1],idNA[k,2]] <- Esti_0[idNA[k,1],idNA[k,2]]
      }
    }
    Residual <- x-Esti
    esti_m <- x[id_m,]
    esti_m_0 <- mo_m$t%*%t(mo_m$P)+mu_X
    if(nrow(idNA_m_id)>0){
      for(k in 1:nrow(idNA_m_id)){
        esti_m[idNA_m_id[k,1],idNA_m_id[k,2]] <- esti_m_0[idNA_m_id[k,1],idNA_m_id[k,2]]
      }
    }
    residual_m <- x[id_m,]-esti_m
    ## Variances of noise and imputation
    Noise <- matrix(0,n,p)
    sds <- rep(0,p)
    for(jj in 1:p){
      idNA_jj_pos <- which(idNA[,2]==jj)
      sds[jj] <- sd(residual_m[,jj],na.rm = T)
      if(is.na(sds[jj])) sds[jj] <- 0
      if(length(idNA_jj_pos)>0){
        idNA_jj <- idNA[idNA_jj_pos,1]
        Esti[idNA_jj,jj] <- Esti[idNA_jj,jj] + rnorm(length(idNA_jj),sd = sds[jj])
      }
    }

    ## Estimation of OOB error
    id_OOB_m <- c(1:n)[-id_m]
    if(length(id_OOB_m)>0){
      x_test <- x[id_OOB_m,,drop=F]
      y_test_oob <- y[id_OOB_m,,drop=F]
      t_oob <- apply(W,2,function(w){prodNA(x_test,w)})
      P_oob <- apply(t,2,function(t_oob){prodNA(t(x),t_oob)/sum(t_oob^2)})
      X_oob <- lapply(1:R,function(r){
        t_oob[,1:r,drop=F]%*%t(P_oob[,1:r,drop=F])
      })
      id_na_OOB <- which(is.na(x_test),arr.ind = T)
      mses <- rep(0,R)
      for(r in 1:R){
        for(jj in 1:p){
          idNA_jj_pos <- which(id_na_OOB[,2]==jj)
          if(length(idNA_jj_pos)>0 & sds[jj]>0){
            idNA_jj <- id_na_OOB[idNA_jj_pos,1]
            X_oob[[r]][idNA_jj,jj] <- X_oob[[r]][idNA_jj,jj] +
              rnorm(length(idNA_jj),sd = sds[jj])
          }
        }
        mu_x <- matrix(rep(mo_m$mux,length(id_OOB_m)),ncol=p,byrow = T)
        mu_y <- matrix(rep(mo_m$muy,length(id_OOB_m)),ncol=q,byrow = T)
        y_est_oob <- (X_oob[[r]]-mu_x)%*%B_R[[r]]+mu_y
        mses[r] <- mean((y_est_oob-y_test_oob)^2)
      }
    }
    Rs[[m]] <- which.min(mses)
    datas_imputed0[[m]] <- list(x=Esti)
    ## Prediction model estimation
    mo_m <- pls_NA(Esti,y,R = Rs[[m]])
    models[[m]] <- mo_m
  }
  models
}

predict_MI_PLS <- function(mos,x_test){
  ress <- lapply(models,function(mo){predict_pls_NA(mo,x_test)})
  y_est <- ress[[1]]
  M <- length(ress)
  if(M>1){
    for(m in 2:M){
      y_est <- y_est + ress[[m]]
    }
    y_est <- y_est/M
  }
  y_est
}
