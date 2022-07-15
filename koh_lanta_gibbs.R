softPCA <- function(x,lambda,minValues=1e-2)
{
  covo <- cov(x)
  diag(covo) <- 0
  co <- abs(covo)-lambda
  co[which(co<0)] <- 0
  co <- co*sign(covo)
  eigen_co <- eigen(co)
  R <- length(which(eigen_co$values>minValues))
  if(R>0)
  {
    w <- eigen_co$vectors[,1:R,drop=F]
    t <- scale(x,scale=F)%*%w
    values <- eigen_co$values[1:R]
  }
  else
  {
    values = t = w <- NULL
  }
  list(w=w,values=values,t=t)
}

getResiduals <- function(x,y,momo)
{
  n <- nrow(x)
  sdy <- momo$model$sdY
  sdx <- momo$model$sdX
  muy <- momo$model$muY
  mux <- momo$model$muX
  resY <- y
  resX <- x
  for(i in 1:n)
  {
    difyi <-y[i,,drop=F]-muy
    difxi <-x[i,,drop=F]-mux
    # yi_pred <- difxi%*%momo$model$B
    resY[i,] <- difyi/sdy#(difyi - pipi%*%difyi)/sdy
    resX[i,] <- difxi/sdx#(difyi - pipi%*%difyi)/sdy
  }
  if(is.null(dim(momo$model$U)))
  {
    pipi <- matrix(0,n,n)
  }
  else
  {
    pipi <- resX%*%momo$model$U
    pipi <- pipi%*%solve(crossprod(pipi))%*%t(pipi)
  }
  projection <- pipi%*%resY
  resY <- resY - projection

  list(residuals=resY,projection=projection)
}

do_SPCA_boot <- function(x_init,y_init,momoXGivenY,n_B=1e2,n_lambdas=20,toPlot=F)
{
  getLambdasHIHI <- function(xSC,ySC,n,p,q){
    getLambda0 <- function(xSC,ySC,n,p,q){
      coco <- cov(xSC,ySC)
      diag(coco) <- 0
      Sig_est <- matrix(rep(coco,n),nrow = n,byrow = T)
      theta <- colMeans((xSC*matrix(rep(ySC,p),n,byrow = F)-Sig_est)^2)
      mean(sqrt(log(max(p,q))*theta/n))
    }
    mean(unlist(lapply(1:q,function(j){getLambda0(xSC,ySC[,j],n,p,q)})))
  }
  n <- nrow(y_init)
  p <- ncol(x_init)
  momoXGivenY <- ddsPLS(y_init,x_init)
  R2 = Q2 = R_est_all <- rep(NA,n_lambdas)
  R_optim = lam_optim <- NULL
  residualisation <- getResiduals(y_init,x_init,momoXGivenY)
  resX_init <- residualisation$residuals
  x_hat_supervized <- residualisation$projection
  if(p>1)
  {
    y_scale <- scale(y_init)
    lambda_0 <- getLambdasHIHI(resX_init,y_init,n,p,q)
    coco <- cov(resX_init,y_scale)
    diag(coco) <- 0
    maxcoco <- max(abs(coco))
    if(maxcoco<lambda_0) lambda_0 <- 0
    lams <- seq(lambda_0,maxcoco,length.out = n_lambdas)
    if(toPlot) par(mfrow=c(1,2))
    flag_On_Continue = T
    i_lam <- 1
    while(flag_On_Continue)
    {
      la <- lams[i_lam]
      q2 = r2 <- rep(NA,n_B)
      q2_r_s = r2_r_s <- list()
      for(ii in 1:n_B)
      {
        id <- sample(1:n,size = n,replace = T)
        id_oob <- which(!(1:n %in% unique(id)))
        x_id <- x_init[id,,drop=F]
        y_id <- y_init[id,,drop=F]
        resX_id <- getResiduals(y_id,x_id,momoXGivenY)$residuals
        mo <- softPCA(resX_id,la)
        if(length(mo$values)>0)
        {
          R_est <- length(mo$values)
          pi_train <- tcrossprod(mo$w)
          proj_train <- resX_id%*%pi_train
          R2_r <- rep(0,R_est)
          for(r in 1:R_est)
          {
            # pi_train_r <- tcrossprod(mo$w[,1:r,drop=F])
            # proj_train_r <- resX_id%*%pi_train_r
            res_unscaled <- resX_id
            t_a <- res_unscaled%*%(mo$w[,1:r,drop=F])
            pipi <- t_a%*%solve(crossprod(t_a))%*%t(t_a)
            proj_train_r <- pipi%*%resX_id#resX_init%*%pipi#
            R2_r[r] <-  1-mean((resX_id-proj_train_r)^2)/mean((resX_id)^2)
          }
          r2_r_s[[ii]] <- R2_r
        }
        else
        {
          proj_train <- matrix(0,nrow(resX_id),ncol(resX_id))
        }
        r2[ii] <- 1-mean((resX_id-proj_train)^2)/mean((resX_id)^2)
        x_T <- x_init[id_oob,,drop=F]
        y_T <- y_init[id_oob,,drop=F]
        resX_oob <- getResiduals(y_T,x_T,momoXGivenY)$residuals
        if(length(mo$values)>0)
        {
          pi_test <- pi_train
          proj <- resX_oob%*%pi_test
          Q2_r <- rep(0,R_est)
          for(r in 1:R_est)
          {
            res_unscaled <- resX_oob
            t_b <- res_unscaled%*%(mo$w[,1:r,drop=F])
            pipi <- t_b%*%solve(crossprod(t_b))%*%t(t_b)
            proj_test_r <- pipi%*%resX_oob#resX_init%*%Pi_pca#
            # pi_test_r <- tcrossprod(mo$w[,1:r,drop=F])
            # proj_test_r <- resX_oob%*%pi_test_r
            Q2_r[r] <-  1-mean((resX_oob-proj_test_r)^2)/mean((resX_oob)^2)
          }
          q2_r_s[[ii]] <- Q2_r
        }
        else
        {
          proj <- matrix(0,nrow(resX_oob),ncol(resX_oob))
        }
        q2[ii] <- 1-mean((resX_oob-proj)^2)/mean((resX_oob)^2)
      }
      R_all <- unlist(lapply(q2_r_s,length))
      tab_R_est <- table(R_all)
      if(dim(tab_R_est)>0){
        R_est = R_est_all[i_lam] <- as.numeric(names(tab_R_est[which.max(tab_R_est)]))
        if(R_est>0)
        {
          id_model_ok <- which(R_all>=R_est)
          R2[i_lam] <- mean(unlist(lapply(r2_r_s[id_model_ok],function(kiki){kiki[R_est]})))
          Q2[i_lam] <- mean(unlist(lapply(q2_r_s[id_model_ok],function(kiki){kiki[R_est]})))
        }
      }
      diffR2Q2 <- abs(R2-Q2)
      if(toPlot)
      {
        matplot(lams,cbind(R2,Q2),col=c(2,4),pch=16,type="b")
        text(lams,(R2+Q2)/2,R_est_all)
        plot(lams,diffR2Q2,col=3,type="b")
        id_min <- which.min(diffR2Q2)
        points(lams[id_min],diffR2Q2[id_min],col=1,pch=16,cex=1.5)
      }
      if(is.na(R2[i_lam]) || i_lam==n_lambdas)
      {
        flag_On_Continue <- F
      }
      else
      {
        i_lam <- i_lam + 1
      }
    }
    res_unscaled <- resX_init
    if(length(na.omit(diffR2Q2))==0)
    {
      pipi <- matrix(0,n,n)
      mo <- list(w=NULL,t=NULL,values=NULL)
    }
    else
    {
      i_lam_optim <- which.min(diffR2Q2)
      lam_optim <- lams[i_lam_optim]
      R_optim <- R_est_all[i_lam_optim]
      mo <- softPCA(resX_init,lam_optim)
      R_optim <- min(R_optim,length(mo$values))
      if(R_optim>0)
      {
        mo$w <- mo$w[,1:R_optim,drop=F]
        mo$t <- mo$t[,1:R_optim,drop=F]
        mo$values <- mo$values[1:R_optim]
        ## Look at final variances
        if(is.null(dim(mo$w)))
        {
          pipi <- matrix(0,p,p)#n,n)
        }
        else
        {
          t_c <- res_unscaled%*%mo$w
          pipi <- t_c%*%solve(crossprod(t_c))%*%t(t_c)
        }
      }
      else
      {
        pipi <- matrix(0,n,n)
        mo <- list(w=NULL,t=NULL,values=NULL)
      }
    }
    x_hat <- x_hat_supervized + pipi%*%res_unscaled#resX_init%*%pipi#
    for(i in 1:n)
    {
      x_hat[i,] <- x_hat[i,]* momoXGivenY$model$sdY + momoXGivenY$model$muY
      x_hat_supervized[i,] <- x_hat_supervized[i,]* momoXGivenY$model$sdY + momoXGivenY$model$muY
    }
    x_residual_final <- x_init - x_hat
    sd_addi_noise <- apply(x_residual_final,2,sd)
  }
  else
  {
    x_hat <- resX_init
    for(i in 1:n)
    {
      x_hat[i,] <- x_hat[i,]* momoXGivenY$model$sdY + momoXGivenY$model$muY
      x_hat_supervized[i,] <- x_hat_supervized[i,]* momoXGivenY$model$sdY + momoXGivenY$model$muY
    }
    x_residual_final <- x_init - x_hat
    lams=R2=Q2=diffR2Q2 <- NULL
    R_optim=lam_optim <- NULL
    mo <- NULL
    sd_addi_noise <- apply(x_residual_final,2,sd)
  }
  list(
    x_hat = x_hat,
    x_hat_supervized = x_hat_supervized,
    x_residual=x_residual_final,
    results=list(
      lambda=lams,
      R2=R2,
      Q2=Q2,
      diffR2Q2=diffR2Q2
    ),
    optimum=list(
      R=R_optim,
      lambda=lam_optim
    ),
    model=mo,
    sd_addi_noise = sd_addi_noise
  )
}

one_cycle_koh_lanta <- function(x,y,id_na_x_id,id_na_y_id,n_B=1e2,n_lambdas=20,toPlot=F)
{
  momoXGivenY <- ddsPLS(y,x,n_B = n_B,n_lambdas = n_lambdas)
  momo_hat_XGivenY <- do_SPCA_boot(x,y,momoXGivenY,n_B=n_B,n_lambdas=n_lambdas,
                                   toPlot=toPlot)
  n_na_x <- nrow(id_na_x_id)
  if(n_na_x>0)
  {
    x_all <- momo_hat_XGivenY$x_hat
    for(i in 1:n_na_x)
    {
      x[id_na_x_id[i,1],id_na_x_id[i,2]] <- x_all[id_na_x_id[i,1],id_na_x_id[i,2]] +
        rnorm(1,sd = momo_hat_XGivenY$sd_addi_noise[id_na_x_id[i,2]])
    }
  }

  momoYGivenX <- ddsPLS(x,y,n_B = n_B,n_lambdas = n_lambdas)
  momo_hat_YGivenX <- do_SPCA_boot(y,x,momoYGivenX,n_B=n_B,n_lambdas=n_lambdas,
                                   toPlot=toPlot)
  n_na_y <- nrow(id_na_y_id)
  if(n_na_y>0)
  {
    y_all <- momo_hat_YGivenX$x_hat
    for(i in 1:n_na_y)
    {
      y[id_na_y_id[i,1],id_na_y_id[i,2]] <- y_all[id_na_y_id[i,1],id_na_y_id[i,2]] +
        rnorm(1,sd = momo_hat_YGivenX$sd_addi_noise[id_na_y_id[i,2]])
    }
  }

  list(data_proper=list(x=x,y=y),
       data_sup_and_unsup=list(x=momo_hat_XGivenY$x_hat,y=momo_hat_YGivenX$x_hat),
       data_sup=list(x=momoXGivenY$Y_est,y=momoYGivenX$Y_est),
       para_sup=list(x=momoXGivenY,y=momoYGivenX),
       para_sup_and_unsup=list(x=momo_hat_XGivenY$model,y=momo_hat_YGivenX$model),
       optimal_para = list(sup=list(x=momoXGivenY$lambda,y=momoYGivenX$lambda),
                           unsup=list(x=momo_hat_XGivenY$optimum$lambda,y=momo_hat_YGivenX$optimum$lambda)
       ),
       sd_residual = list(x=momo_hat_XGivenY$sd_addi_noise,y=momo_hat_YGivenX$sd_addi_noise)
  )
}

nouveau_koh_lanta_m <- function(x,y,n_B,n_lambdas,maxIter,errorMin)
{
  n <- nrow(x)
  id_na_x <- which(is.na(x))
  n_na_x <- length(id_na_x)
  id_na_x_id <- which(is.na(x),arr.ind = T)
  id_na_y <- which(is.na(y))
  n_na_y <- length(id_na_y)
  id_na_y_id <- which(is.na(y),arr.ind = T)

  mux <- colMeans(x,na.rm=T)
  muy <- colMeans(y,na.rm=T)

  x0 <- x
  y0 <- y
  for(i in 1:n_na_x)
  {
    x0[id_na_x_id[i,1],id_na_x_id[i,2]] <- mux[id_na_x_id[i,2]]
  }

  for(i in 1:n_na_y)
  {
    y0[id_na_y_id[i,1],id_na_y_id[i,2]] <- muy[id_na_y_id[i,2]]
  }

  ## Gibbs
  test <- T
  cycle <- 1
  while(test)
  {
    one_cycle <- one_cycle_koh_lanta(x0,y0,id_na_x_id,id_na_y_id,n_B=n_B,n_lambdas=n_lambdas,toPlot=F)
    x0 <- one_cycle$data_proper$x
    y0 <- one_cycle$data_proper$y
    if(cycle==maxIter ){#| errors[cycle] < errorMin){
      test <- F
    }else{
      cycle <- cycle + 1
    }
  }
  one_cycle
}

nouveau_koh_lanta_all_M <- function(x,y,M=10,n_B=20,n_lambdas=10,maxIter=5,errorMin=1e-2)
{
  n <- nrow(x)
  id_na_x <- which(is.na(x))
  n_na_x <- length(id_na_x)
  id_na_x_id <- which(is.na(x),arr.ind = T)
  id_na_y <- which(is.na(y))
  n_na_y <- length(id_na_y)
  id_na_y_id <- which(is.na(y),arr.ind = T)

  mux <- colMeans(x,na.rm=T)
  muy <- colMeans(y,na.rm=T)

  results <- list()
  data_imputed <- list()
  n <- nrow(y)
  for(m in 1:M)
  {
    cat("\n--------------- ");cat(m);cat("\n")
    id_m <- sample(1:n,n,replace=T)
    ## Estimate new model from bootstrap sample
    results[[m]] <- nouveau_koh_lanta_m(x[id_m,,drop=F],y[id_m,,drop=F],n_B,
                                        n_lambdas,maxIter,errorMin)
    ## Gibbs
    test <- T
    cycle <- 1
    xm <- x
    ym <- y
    if(n_na_x>0)
    {
      for(i in 1:n_na_x)
      {
        xm[id_na_x_id[i,1],id_na_x_id[i,2]] <- mux[id_na_x_id[i,2]]
      }
    }
    if(n_na_y>0)
    {
      for(i in 1:n_na_y)
      {
        ym[id_na_y_id[i,1],id_na_y_id[i,2]] <- muy[id_na_y_id[i,2]]
      }
    }
    ## Estimate new data from chosen model
    while(test)
    {
      ## Side X
      xm_sup <- predict(results[[m]]$para_sup$x,ym)$y_est
      residual_xm <- xm-xm_sup
      u_x <- results[[m]]$para_sup_and_unsup$x$w
      t_res <- residual_xm%*%u_x
      pipi <- t_res%*%solve(crossprod(t_res))%*%t(t_res)
      xm_sup_and_unsup <- xm_sup + pipi%*%residual_xm#residual_xm%*%pipi#
      xm0 <- xm
      if(n_na_x>0)
      {
        for(i in 1:n_na_x)
        {
          xm[id_na_x_id[i,1],id_na_x_id[i,2]] <-
            xm_sup_and_unsup[id_na_x_id[i,1],id_na_x_id[i,2]] +
            rnorm(1,sd=results[[m]]$sd_residual$x[id_na_x_id[i,2]])
        }
      }
      ## Side Y
      ym_sup <- predict(results[[m]]$para_sup$y,xm)$y_est
      residual_ym <- ym-ym_sup
      u_y <- results[[m]]$para_sup_and_unsup$y$w
      if(is.null(dim(u_y)))
      {
        pipi <- matrix(0,n,n)
      }
      else
      {
        t_res <- residual_ym%*%u_y
        pipi <- t_res%*%solve(crossprod(t_res))%*%t(t_res)
        # pipi <- tcrossprod(residual_ym%*%u_y)
      }
      ym_sup_and_unsup <- ym_sup + pipi%*%residual_ym#residual_ym%*%pipi#
      if(n_na_y>0)
      {
        for(i in 1:n_na_y)
        {
          ym[id_na_y_id[i,1],id_na_y_id[i,2]] <-
            ym_sup_and_unsup[id_na_y_id[i,1],id_na_y_id[i,2]] +
            rnorm(1,sd=results[[m]]$sd_residual$y[id_na_y_id[i,2]])
        }
      }

      if(cycle==maxIter ){#| errors[cycle] < errorMin){
        test <- F
      }else{
        cycle <- cycle + 1
      }
    }
    data_imputed[[m]] <- list(x=xm,y=ym)
  }
  data_imputed
}
