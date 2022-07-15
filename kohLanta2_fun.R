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

nipals_pca_full <- function(x0,R,n,p){
  x <- scale(x0,scale = F)
  ts <- matrix(NA,n,R)
  W = P <- matrix(NA,p,R)
  for(r in 1:R){
    w0 <- rnorm(p)
    w0 <- w0/sqrt(sum(w0^2))
    test <- T
    while(test){
      t <- apply(x,1,function(xi){sum(na.omit(xi*w0))})
      Pr <- apply(x,2,function(xj){sum(na.omit(xj*t))})/sum(t^2)
      w <- Pr/sqrt(sum(Pr^2))
      test <- sum((w-w0)^2)>1e-5
      w0 <- w
    }
    W[,r] <- w
    P[,r] <- Pr
    ts[,r] <- t
  }
  xout <- x0
  xReconst <- tcrossprod(ts,P)
  id_na <- which(is.na(x0),arr.ind = T)
  mus <- colMeans(x0,na.rm = T)
  for(i in 1:nrow(id_na)){
    id <- id_na[i,]
    xout[id[1],id[2]] <- mus[id[2]] + xReconst[id[1],id[2]]
  }
  xout
}

#' Title
#'
#' @param X
#' @param y
#' @param B
#' @param M Number of multiple imputations
#' @param LD Wether use low dimensional approach
#' @param init How to initialize imputation: "mean" or "nipals". "nipals" corresponds to full rank Nipals PCA reconstruction
#' @param n_lambdas
#' @param symmetric
#' @param errorMinEM
#' @param maxIterEm
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
KohLanta <- function(X,y,B=500,M=20,LD=F,init=c("mean","nipals"),n_lambdas=30,symmetric=F,
                     errorMinEM=1e-2,maxIterEm=20,verbose=F){
  n <- nrow(y)
  p <- ncol(X)
  q <- ncol(y)
  id_IB_s <- matrix(sample(1:n,size = n*M,replace = T),nrow = M)
  x_m_s <- list()
  Rs <- list()
  R_imputation = R_prediction <- rep(0,M)
  lambdas <- list()
  sigma2_k_s <- list()
  var_x <- apply(X,2,var,na.rm=T)
  sd_residuals <- list()
  models <- list()
  model_imputation <- list()

  x_m0 <- X
  y_m0 <- y
  id_NA <- which(is.na(X),arr.ind = T)
  id_NA_nono <- which(is.na(X),arr.ind = F)
  MUy <- colMeans(y_m0,na.rm = T)
  denomY <- solve(crossprod(scale(y_m0,scale = F)))
  sigma2 <- apply(x_m0,2,var,na.rm = T)
  MU <- colMeans(x_m0,na.rm = T)
  if(nrow(id_NA)>0){
    if(init=="mean"){
      for( jj in unique(id_NA[,2])){
        idNA <- id_NA[which(id_NA[,2]==jj),1]
        x_m0[idNA,jj] <- MU[jj]
      }
    }else{
      x_m0 <- nipals_pca_full(X,R = min(n-1,p-1),n = n,p = p)
    }
  }
  defaultW <- getOption("warn")
  options(warn = -1)
  cor_max_total <- max(abs(cor(x_m0,y_m0)))
  options(warn = defaultW)
  idNACorM <- which(is.na(cor_max_total))
  if(length(which(is.na(cor_max_total)))>0){
    cor_max_total[idNACorM] <- 0
  }
  lambdas_imput <- seq(0,cor_max_total,length.out = n_lambdas+1)[-n_lambdas-1]
  for(m in 1:M){
    if(verbose==T){cat("m=");cat(m);cat(" ")}
    x_m <- X[id_IB_s[m,],,drop=F]
    y_m <- y[id_IB_s[m,],,drop=F]
    id_na_m <- which(is.na(x_m),arr.ind = T)
    id_na_m_nono <- which(is.na(x_m),arr.ind = F)
    Rs[[m]] <- 0
    lambdas[[m]] <- NULL
    sigma2_k_s[[m]] <- apply(x_m,2,var,na.rm = T)#list()
    id_na_s2 <- which(is.na(sigma2_k_s[[m]]))
    if(length(id_na_s2)>0){
      sigma2_k_s[[m]][id_na_s2] <- 0
    }
    MU <- colMeans(x_m,na.rm = T)
    ## Impute to mean or nipalsPCA
    MUy <- colMeans(y_m,na.rm = T)
    denomY <- solve(crossprod(scale(y_m,scale = F)))
    if(nrow(id_na_m)>0){
      if(init=="mean"){
        for( jj in unique(id_na_m[,2])){
          idNA <- id_na_m[which(id_na_m[,2]==jj),1]
          x_m[idNA,jj] <- MU[jj]

        }
      }else{
        x_m <- nipals_pca_full(x_m,R = min(n-1,p-1),n = n,p = p)
      }
      defaultW <- getOption("warn")
      options(warn = -1)
      cocoM <- cor(x_m,y_m)
      options(warn = defaultW)
      idNACocoM <- which(is.na(cocoM))
      if(length(which(is.na(cocoM)))>0){
        cocoM[idNACocoM] <- 0
      }

      lambda_max <- min(max(abs(cocoM)),cor_max_total)

      ########################
      # Parameter estimation #
      ########################

      test_em <- TRUE
      it_em <- 0
      R0 <- 0
      errors = R_m <- NULL
      lambdas_m <- list()
      lambda_m_est <- NULL
      while(test_em & it_em<maxIterEm){
        if(verbose==T) cat("-")
        it_em <- it_em + 1
        model_inverse <- ddsPLS(y_m,x_m,n_B=B,lambdas = lambdas_imput,#NULL,
                                lambda_roof=lambda_max,LD=LD,
                                criterion = "diffR2Q2",#"Q2",#
                                n_lambdas = n_lambdas)
        x_m_est <- model_inverse$Y_est[id_na_m]
        sd2_na_m <- sigma2_k_s[[m]][id_na_m[,2]]
        errs <- (x_m_est-x_m[id_na_m])^2
        error_em <- rep(0,length(errs))
        for(iii in 1:length(errs)){
          if(is.na(sd2_na_m[iii])) browser()
          if(sd2_na_m[iii]>0){
            error_em[iii] <- errs[iii]/sd2_na_m[iii]
          }
        }
        error_em <- mean(error_em)
        # error_em <- mean(/
        #                    sigma2_k_s[[m]][id_na_m[,2]])#/n_selectedX
        errors <- c(errors,error_em)
        R0 <- model_inverse$R
        lambdas_m[[it_em]] <- model_inverse$lambda
        R_m <- c(R_m,R0)
        model0 <- model_inverse
        x_m[id_na_m] <- x_m_est
        if(error_em<errorMinEM){
          test_em<-FALSE
        }
        lambda_m_est <- model_inverse$lambda
      }
      lambdas[[m]] <- lambda_m_est
      model_imputation[[m]] <- model_inverse


      ##############
      # IMPUTATION #
      ##############

      x_0 <- x_m0
      y_0 <- y_m0

      x_0[id_NA_nono] <- predict(model_inverse,y_0)$y_est[id_NA_nono]

      #############################
      ## Estimate additive noise ##
      #############################

      Residual_m <- X[id_IB_s[m,],,drop=F]-model_inverse$Y_est#x_m#_est
      x_test <- X-x_0
      sd_residuals_m <- rep(0,p)
      for(jj in 1:p){
        sd_residuals_m[jj] <- sd(Residual_m[,jj],na.rm = T)
        if(is.na(sd_residuals_m[jj])) sd_residuals_m[jj] <- 0
        idNA_jj <- which(is.na(X[,jj]))
        if(length(idNA_jj)>0){
          val <- rnorm(length(idNA_jj),sd=sd_residuals_m[jj])
          x_0[idNA_jj,jj] <- x_0[idNA_jj,jj] + val
          x_test[idNA_jj,jj] <- val
        }
      }

      if(m==1) x_ref <- x_0

      if(verbose==T)cat(paste("   Imputation (R=",model_inverse$R,", Q2=",round(model_inverse$Q2[model_inverse$R],2),")\n",sep=""))
    }else{
      x_ref <- x_m
      x_0 <- x_m0
      y_0 <- y_m0
      sd_residuals_m <- NULL
    }
    if(!symmetric){
      model <- ddsPLS(x_0,y_0,doBoot = T,lambdas = NULL,LD=LD,
                      n_B = B,n_lambdas = n_lambdas)
      if(m==1){
        model_ref <- ddsPLS(x_ref,y_0,doBoot = T,lambdas = NULL,LD=LD,
                            n_B = B,n_lambdas = n_lambdas)
      }
    }else{
      model <- ddsPLS(x_0,y_0,doBoot = F,lambdas = lambdas[[m]],LD=LD,
                      n_B = B,n_lambdas = n_lambdas)
      if(m==1){
        model_ref <- ddsPLS(x_ref,y_0,doBoot = F,lambdas = lambdas[[m]],LD=LD,
                            n_B = B,n_lambdas = n_lambdas)
      }
    }
    if(verbose==T) cat(paste("   Prediction (R=",model$R,", Q2=",round(model$Q2[model$R],2),")\n",sep=""))
    x_m_s[[m]] <- x_0
    models[[m]] <- model
    R_prediction[m] <- model$R
    sd_residuals[[m]] <- sd_residuals_m
  }
  R_out <- list(Imputation=R_imputation,Prediction=R_prediction)
  out <- list(original=X,imputed=x_m_s,
              models=models,model_imputation=model_imputation,
              R=R_out)
  out$lambda <- lapply(1:M,function(m){out$models[[m]]$lambda})
  q <- ncol(y)
  out$Statistics <- build_statistics(out,n,p,q,M,B)
  out$Q2 <- lapply(1:M,function(m){out$models[[m]]$Q2})
  out$Q2h <- lapply(1:M,function(m){out$models[[m]]$Q2h})
  out$R2 <- lapply(1:M,function(m){out$models[[m]]$R2})
  out$R2h <- lapply(1:M,function(m){out$models[[m]]$R2h})
  out$varExplained <- lapply(1:M,function(m){out$models[[m]]$varExplained})
  out$Selection <- lapply(1:M,function(m){out$models[[m]]$Selection})
  out$Reference <- list(model=model_ref,imputed=x_ref)
  class(out) <- "KohLanta"
  return(out)
}

#' predict.KohLanta
#'
#' @param x model
#' @param X_test new data set
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.KohLanta <- function(x,X_test=NULL,...){
  M <- length(x$models)
  n_test <- nrow(X_test)
  Q2_M <- do.call(rbind,lapply(1:M,function(m){
    mo <- x$models[[m]]
    Rm <- x$R$Prediction[m]
    if(Rm>0){
      out <- mo$Q2[Rm]
    }else{
      out <- 0
    }
    c(Rm,out)
  }))
  meanQ2 <- do.call(rbind,lapply(sort(unique(x$R$Prediction)),function(r){
    c(r,mean(Q2_M[which(x$R$Prediction==r),2]))
  }))
  R_ok <- meanQ2[which(meanQ2[,2]>0),1]
  if(length(data_test$x)>0){
    id <- 1:M
    if(length(R_ok)>0){
      id <- which(Q2_M[,1] %in% R_ok)
      out <- list()
      y_est <- list()
      for(m in 1:M){
        if(is.null(X_test)){
          y_est[[m]] <- x$models[[m]]$Y_est
        }else{
          y_est[[m]] <- predict(x$models[[m]],X_test)$y_est
        }
      }
      denom <- sum(Q2_M[id,2])
      m <- id[1]
      out_sum <- y_est[[m]]*Q2_M[m,2]/denom
      for(i_m in 1:length(id)){
        m <- id[i_m]
        out[[i_m]] <- y_est[[m]]
        if(i_m>1) out_sum <- out_sum + out[[i_m]]*Q2_M[i_m,2]/denom
      }
    }else{
      out = y_est <- list()
      out_sum <- matrix(0,n_test,length(x$models[[m]]$model$muY))
      for(m in 1:M){
        y_est[[m]] <- matrix(rep(x$models[[m]]$model$muY,n_test),nrow = n_test,byrow = T)
        out[[m]] <- y_est[[m]]
        out_sum <- out_sum + out[[m]]/m
      }
    }
  }
  list(est=out_sum,distri=out)
}

#' Function to plot bootstrap performance results of the KohLanta algorithm
#'
#' @param x A KohLanta object.
#' @param R The number of components of the models desired, default to majority.
#' @param type The type of graphics. One of `B`, `U`, `V`, `P`, `C`, `U_star`,
#'  `T`, `selectionComp` or `selectionGlobal`.
#' @param whether or not to link the dots in the plot. Default to FALSE.
#' @param digits double. Rounding of the written explained variance.
#' @param legend.position character. Where to put the legend.
#' @param horiz boolean. Whether to plot horizontally.
#' @param col vector. Mainly to modify bars in weight plots.
#' @param mar vector. The margins for the plot.
#' @param cex.names double. Size factor for variable names.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @importFrom graphics layout
#'
#' @export
#'
#' @seealso \code{\link{KohLanta}}, \code{\link{ddsPLS}}
#'
#' @useDynLib KohLanta
plot.KohLanta <- function(x,R=NULL,stat="criterion",
                          bivariate=FALSE,
                          procruste=FALSE,
                          components=c(1,2),
                          addLine=FALSE,
                          digits=1,
                          legend.position="topright",
                          horiz=TRUE,
                          col=NULL,
                          cex.names=1,mar=c(5, 4, 4, 2) + 0.1,
                          ...){
  ploplot <- function(im,stat,addLine=FALSE){

    ## -----------------------------------
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
    Rs <- sort(unique(im$R$Prediction))
    if(any(Rs==0)) Rs <- Rs[-1]
    q <- ncol(im$models[[1]]$Y_obs )
    p <- ncol(im$original)
    if(stat=="selectionGlobal"){
      par(mfrow=c(1,length(Rs)))
      par(mar=c(3,3,2,1))
      countsR <- lapply(Rs,function(R){
        idR <- which(im$R$Prediction==R)
        moR <- im$models[idR]
        MR <- length(idR)
        stats <- table(unlist(lapply(1:MR,function(m){moR[[m]]$Selection$X})))
      })
      ylim <- c(0,max(unlist(countsR)))
      for(iR in 1:length(Rs)){
        R <- Rs[iR]
        counts <- rep(0,p)
        stats <- countsR[[iR]]
        counts[as.numeric(names(stats))] <- stats
        bb <- barplot(counts,ylim=ylim,xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
        title(main = paste("Variable selected (X) on",R,"cp(s)."),line = 1)
        title(xlab = "Id",line = 2)
        title(ylab = "Count",line = 2)
        labY2 <- unique(round(seq(0,ylim[2],length.out = 7)))
        axis(2,labY2,tick = T,at = labY2,line = 0,las=2)
        axis(1,1:p,tick = F,at = bb,line = 0)
      }
    }else if (stat=="selectionComp"){
      par(mfrow=c(length(Rs),max(Rs)))
      par(mar=c(3,3,3,1))
      countsR <- lapply(Rs,function(R){
        idR <- which(im$R$Prediction==R)
        moR <- im$models[idR]
        MR <- length(idR)
        stats <- lapply(1:MR,function(m){
          U <- moR[[1]]$model$U
          apply(abs(U),2,function(u){which(u>1e-5)})
        })
        lapply(1:R,function(r){
          table(unlist(lapply(1:MR,function(m){stats[[m]][[r]]})))
        })
      })
      ylim <- c(0,max(unlist(countsR)))
      for(iR in 1:length(Rs)){
        R <- Rs[iR]
        for(r in 1:R){
          counts <- rep(0,p)
          stats <- countsR[[iR]][[r]]
          counts[as.numeric(names(stats))] <- stats
          bb <- barplot(counts,ylim=ylim,xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
          title(main = paste("Variable selected (X) on cp",r,"\n among",R,"cp(s)."),line = 1)
          title(xlab = "Id",line = 2)
          title(ylab = "Count",line = 2)
          labY2 <- unique(round(seq(0,ylim[2],length.out = 7)))
          axis(2,labY2,tick = T,at = labY2,line = 0,las=2)
          axis(1,1:p,tick = F,at = bb,line = 0)
        }
        if(R<max(Rs)){
          for(r in (R+1):max(Rs)){
            plot(1,1,col="white",xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
          }
        }
      }
    }else if (stat!="B"){
      par(mfrow=c(length(Rs),max(Rs)))
      par(mar=c(2,2,2,1))
      for(R in Rs){
        id <- which(names(im$Statistics[[R]])==stat)
        limi <- 1.15#max(abs(do.call(cbind,im$Statistics[[R]][[id]])))*1.15
        for(r in 1:R){
          mu <- im$Statistics[[R]][[id]]$mean[,r]
          sd <- sqrt(im$Statistics[[R]][[id]]$var)[,r]
          ddff <- cbind(mu,mu+sd,mu-sd)
          matplot(ddff,type="p",pch=c(16,NA,NA),lty=c(1,1,1),
                  col=r,ylim=c(-1,1)*limi,
                  xlim=c(0,nrow(ddff)+1),
                  xlab="",ylab="",
                  xaxt="n",yaxt="n",
                  main="",bty="n")
          title(main = paste(stat,", cp",r,"on",R,"cp(s)."),line = 0)
          title(xlab = "Id",line = 1)
          # title(ylab = "Value",line = 1)
          bb <- round(range(c(-1,1)*limi))
          axis(1,unique(round(seq(1,nrow(ddff),length.out = 5))),tick = T,
               at = unique(round(seq(1,nrow(ddff),length.out = 5))),line = -3/4)
          axis(2,seq(bb[1],bb[2],length.out = 3),tick = T,
               at = seq(bb[1],bb[2],length.out = 3),line = -3/4)
          for(i in 1:length(mu)){
            points(c(i-1/4,i+1/4),(mu+sd)[i]*c(1,1),type="l",col=r)
            points(c(i-1/4,i+1/4),(mu-sd)[i]*c(1,1),type="l",col=r)
            points(c(i,i),c((mu-sd)[i],(mu+sd)[i]),type="l",col=r)
          }
          if(addLine){
            points(mu,type="l",col=r,lwd=2)
          }
          abline(h=0,lty=3,lwd=1/2)
          abline(h=c(-1,1),lty=1,lwd=1/2)
        }
        if(R<max(Rs)){
          for(r in (R+1):max(Rs)){
            plot(0,0,xlab="",xaxt="n",ylab="",yaxt="n",
                 bty="n",col="white")
          }
        }
      }
    }else{
      par(mfrow=c(length(Rs),q))
      par(mar=c(2,2,2,1))
      for(R in Rs){
        id <- which(names(im$Statistics[[R]])==stat)
        q <- ncol(im$Statistics[[R]][[id]]$mean)
        limi <- max(abs(do.call(cbind,im$Statistics[[R]][[id]])))
        for(j in 1:q){
          mu <- im$Statistics[[R]][[id]]$mean[,j]
          sd <- sqrt(im$Statistics[[R]][[id]]$var)[,j]
          ddff <- cbind(mu,mu+sd,mu-sd)
          matplot(ddff,type="p",pch=c(16,NA,NA),lty=c(1,1,1),
                  col=j,ylim=c(-1,1)*limi,xlab="",ylab="",
                  main="",bty="n",xaxt="n",yaxt="n")
          title(main = paste(stat,"for variable",j,"on",R,"cps."),line = 0)
          title(xlab = "Id",line = 1)
          # title(ylab = "Value",line = 1)
          bb <- round(range(c(-1,1)*limi))
          axis(1,unique(round(seq(1,nrow(ddff),length.out = 5))),tick = T,
               at = unique(round(seq(1,nrow(ddff),length.out = 5))),line = -3/4)
          axis(2,seq(bb[1],bb[2],length.out = 3),tick = T,
               at = seq(bb[1],bb[2],length.out = 3),line = -3/4)
          for(i in 1:length(mu)){
            points(c(i-1/4,i+1/4),(mu+sd)[i]*c(1,1),type="l",col=j)
            points(c(i-1/4,i+1/4),(mu-sd)[i]*c(1,1),type="l",col=j)
            points(c(i,i),c((mu-sd)[i],(mu+sd)[i]),type="l",col=j)
          }
          if(addLine){
            points(mu,type="l",col=j,lwd=2)
          }
          abline(h=0,lty=3)
        }
      }
    }
  }
  bivarPlot <- function(im,stat,cps=c(1,2),procruste=F,
                        addLine=FALSE){
    ## -----------------------------------
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
    Rs <- sort(unique(im$R$Prediction))
    if(stat!="B"){
      K <- 5
      m0 <- matrix(rep(c(1:length(Rs)),K),nrow=K,byrow = T)
      if(stat!="T"){
        layout(mat = rbind(m0,rep(length(Rs)+1,length(Rs)) ))
      }else{
        layout(mat = m0)
      }
      for(R in Rs){
        if(R<max(cps)){
          plot(0,0,xaxt="n",xlab="",yaxt="n",
               ylab="",main="",bty="n",col="white")
          text(0,0,labels = "Selected components not built in this model")
        }else{
          id <- which(names(im$Statistics[[R]])==stat)
          idR <- which(im$R$Prediction==R)
          varExpl <- do.call(rbind,lapply(im$models[idR],function(ll){
            oo <- ll$varExplained$Comp[cps]
            if(length(oo)<2){oo <- c(oo,0)}
            else if(length(oo)>2){oo <- oo[1:2]}
            oo}))
          muVarexpl <- colMeans(varExpl)
          sdVarexpl <- apply(varExpl,2,sd)
          if(stat=="T"){
            limi <- max(abs(cbind(im$Statistics[[R]][[id]]$mean,
                                  sqrt(im$Statistics[[R]][[id]]$var))))*c(-1,1)*1.2
          }else{
            limi <- 1.2*c(-1,1)
          }
          mu1 <- im$Statistics[[R]][[id]]$mean[,cps[1]]
          sd1 <- sqrt(im$Statistics[[R]][[id]]$var)[,cps[1]]
          mu2 <- im$Statistics[[R]][[id]]$mean[,cps[2]]
          sd2 <- sqrt(im$Statistics[[R]][[id]]$var)[,cps[2]]
          colo <- color[c(1:p)%%(length(color)+1)]
          if(stat=="T"){
            colo <- "white"
            coloS <- color[c(1:n)%%(length(color)+1)]
          }
          plot(mu1,mu2,
               xlab=paste("Comp ",cps[1]," (",round(muVarexpl)[1],
                          "% +-",round(sdVarexpl,2)[1]," )",sep=""),
               ylab=paste("Comp ",cps[2]," (",round(muVarexpl)[2],
                          "% +-",round(sdVarexpl,2)[2]," )",sep=""),
               xlim=limi,ylim=limi,col=colo,pch=1,
               main=paste(stat,", model on ",R," components",sep=""))
          if(stat=="T"){
            text(mu1,mu2,labels = 1:length(mu1),col=coloS)
            colo <- coloS
          }
          W <- diff(limi)/100
          for(i in 1:length(mu1)){
            xi <- mu1[i];sdxi <- sd1[i];yi <- mu2[i];sdyi <- sd2[i]
            points(c(xi-W,xi+W),(yi+sdyi)*c(1,1),type="l",col=colo[i])
            points(c(xi-W,xi+W),(yi-sdyi)*c(1,1),type="l",col=colo[i])
            points(c(xi,xi),c(yi-sdyi,yi+sdyi),type="l",col=colo[i],lwd=1/2)
            yi <- mu1[i];sdyi <- sd1[i];xi <- mu2[i];sdxi <- sd2[i]
            points((yi+sdyi)*c(1,1),c(xi-W,xi+W),type="l",col=colo[i])
            points((yi-sdyi)*c(1,1),c(xi-W,xi+W),type="l",col=colo[i])
            points(c(yi-sdyi,yi+sdyi),c(xi,xi),type="l",col=colo[i],lwd=1/2)
          }
          abline(h=0,v=0,lty=3,col="gray")
          if(stat!="T"){
            th <- seq(-pi,pi,length.out = 100)
            points(cos(th),sin(th),col="gray",type="l")
            points(cos(th)/2,sin(th)/2,col="gray",type="l",lty=2)
          }
        }
      }
      if(R>=max(cps)){
        if(stat!="T"){
          colo <- color[c(1:p)%%(length(color)+1)]
          par(mar=c(0,0,0,0))
          if(!(stat %in% c("C","V"))){
            p <- ncol(im$original)
            plot(1:p,rep(0,p),col=colo,bty="n",pch=16,
                 xaxt="n",yaxt="n",xlab="",ylab="")
            title(main = "Variable legend",line = -1 )
            text(1:p,0,pos = 1,labels = 1:p)
          }else{
            colo <- color[c(1:n)%%(length(color)+1)]
            p <- ncol(im$Reference$model$Y_obs)
            plot(1:p,rep(0,p),col=colo,bty="n",pch=16,
                 xaxt="n",yaxt="n",xlab="",ylab="")
            title(main = "Variable legend",line = -1 )
            text(1:p,0,pos = 1,labels = 1:p)
          }
        }
      }
    }
  }

  plot_ellipse <- function(dd,main,col=1,get=FALSE,plot=T){
    coco <- var(dd,na.rm = T)
    eiei <- eigen(coco)
    angle <- atan2(eiei$vectors[2,1],eiei$vectors[1,1])
    if(angle < 0) angle <- angle + 2*pi
    avg = colMeans(dd)
    chisquare_val = 2.4477
    theta_grid = seq(0,2*pi,length.out = 180)
    phi = angle
    X0=avg[1]
    Y0=avg[2]
    a <- chisquare_val*sqrt(eiei$values[1])
    b <- chisquare_val*sqrt(eiei$values[2])
    ellipse_x_r  <- a*cos( theta_grid )
    ellipse_y_r  <- b*sin( theta_grid )
    R <- cbind(c(cos(phi),sin(phi)),c(-sin(phi), cos(phi)))
    r_ellipse = cbind(ellipse_x_r,ellipse_y_r)%*% R
    out <- cbind(X0+r_ellipse[,1],Y0+r_ellipse[,2])
    if(plot){
      points(X0+r_ellipse[,1],Y0+r_ellipse[,2],type="l",col=col)
      text(X0,Y0,main,col=col)
    }
    if(get){
      return(out)
    }
  }

  color = c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon")

  procruste_pages <- function(ref,new){
    decomp <- svd(crossprod(new,ref))
    alpha <- sum(new^2)/sum(decomp$d)
    new%*%tcrossprod(decomp$u,decomp$v)/alpha
  }

  plot.im.indiv <- function(im){
    ## Reset personnal plot par() settings
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
    id_na_ix <- which(is.na(im$original),arr.ind = T)
    people_ix <- unique(id_na_ix[,1])
    t_ref <- im$Reference$model$model$t
    varExpl <- do.call(rbind,lapply(im$models,function(ll){
      oo <- ll$varExplained$Comp
      if(length(oo)<2){oo <- c(oo,0)}
      else if(length(oo)>2){oo <- oo[1:2]}
      oo}))
    muVarexpl <- colMeans(varExpl)
    sdVarexpl <- apply(varExpl,2,sd)
    M <- length(im$imputed)
    scores_1 <- list()
    propro <- list()
    for(m in 1:M){
      t_m <- im$models[[m]]$model$t
      t_m_rotate <- procruste_pages(t_ref,t_m)
      scores_1[[m]] <- t_m_rotate
    }
    for(i in 1:n){
      dd <- na.omit(do.call(rbind,lapply(scores_1,function(sc){sc[i,]})))
      propro[[i]] <- plot_ellipse(dd,paste(i),get=T,plot=F)
    }
    plot(t_ref,col="white",
         xlim=c(-1,1)*max(abs(do.call(rbind,propro)[,1])),
         ylim=c(-1,1)*max(abs(do.call(rbind,propro)[,2])),
         xlab=paste("Component 1 (",round(muVarexpl)[1],"% +-",round(sdVarexpl,2)[1]," )"),
         ylab=paste("Component 2 (",round(muVarexpl)[2],"% +-",round(sdVarexpl,2)[2]," )"),
         main="Scores")
    for(i in 1:n){
      if(i %in% people_ix){
        col <- color[which(people_ix==i)]
        pm <- propro[[i]]
        # plot_ellipse(dd,paste(i),col=col)
        points(pm[,1],pm[,2],type="l",col=col)
        text(mean(pm[,1]),mean(pm[,2]),i,col=col)
      }else{
        dd <- na.omit(do.call(rbind,lapply(scores_1,function(sc){sc[i,]})))
        mu_i <- colMeans(dd)
        text(mu_i[1],mu_i[2],paste(i),cex = 1/2)
      }
    }
  }

  plot.im.variables <- function(im,type="weight"){
    ## Reset personnal plot par() settings
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
    p <- ncol(im$original)
    if(type%in%c("weight","loading")){
      varExpl <- do.call(rbind,lapply(im$models,function(ll){
        oo <- ll$varExplained$Comp
        if(length(oo)<2){oo <- c(oo,0)}
        else if(length(oo)>2){oo <- oo[1:2]}
        oo}))
      muVarexpl <- colMeans(varExpl)
      sdVarexpl <- apply(varExpl,2,sd)
      if(type=="loading"){
        P_ref <- im$Reference$model$model$P
      }else if(type=="weight"){
        P_ref <- im$Reference$model$model$U
      }
      M <- length(im$imputed)
      loadings_1 <- list()
      if(type=="loading"){
        main <- "Loadings"
      }else if(type=="weight"){
        main <- "Weights"
      }
      for(m in 1:M){
        if(type=="loading"){
          P_m <- im$models[[m]]$model$P
        }else if(type=="weight"){
          P_m <- im$models[[m]]$model$U
        }
        P_m_rotate <- procruste_pages(P_ref,P_m)
        loadings_1[[m]] <- P_m_rotate
      }
      theta <- seq(-pi,pi,length.out = 200)
      layout(t(t(c(rep(1,7),2))))
      xlim <- max(1,max(abs(unlist(lapply(loadings_1,function(ll){range(ll[,1])})))))*c(-1,1)
      ylim <- max(1,max(abs(unlist(lapply(loadings_1,function(ll){range(ll[,2])})))))*c(-1,1)
      plot(cos(theta),sin(theta),col="gray",xlim=xlim,ylim=ylim,
           xlab=paste("Component 1 (",round(muVarexpl)[1],"% +-",round(sdVarexpl,2)[1]," )"),
           ylab=paste("Component 2 (",round(muVarexpl)[2],"% +-",round(sdVarexpl,2)[2]," )"),
           main=main,type="l")
      points(cos(theta)/2,sin(theta)/2,col="gray",lty=2,type="l")
      abline(h=0,v=0,col="gray")
      for(m in 1:M){
        col <- color
        points(loadings_1[[m]],pch=16,col=color)
      }
      par(mar=c(0,0,0,0))
      p <- ncol(im$original)
      plot(1:p,rep(0,p),col=color[1:p],bty="n",pch=16,
           xaxt="n",yaxt="n",xlab="",ylab="")
      title(main = "Variable legend",line = -1 )
      text(1:p,0,pos = 1,labels = 1:p)
    }else{
      cat(
        "Please select a correct stat among `weight` or `loading`.")
    }
  }
  ## -----------------------------------
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  tryCatch(
    expr = {
      if(!bivariate){
        ploplot(x,stat=stat,addLine = addLine)

      }else{
        bivarPlot(x,stat,cps=components,procruste=procruste,
                  addLine=addLine)
      }},
    error = function(e){
      cat(
        "Please select a correct stat among `B`, `U`, `V`, `P`, `C`,
          `U_star`, `T`, `selectionComp` or `selectionGlobal`")
    }
  )
}

build_statistics <- function(im,n,p,q,M,B){
  unique_R_prediction <- sort(unique(im$R$Prediction))
  out_total <- list()
  for(r in unique_R_prediction){
    if(r!=0){
      id_r <- which(im$R$Prediction==r)
      if(length(id_r)>0){
        momo <- im$models[id_r]
        p <- ncol(momo[[1]]$X)
        q <- ncol(momo[[1]]$Y_est)
        n <- nrow(momo[[1]]$Y_est)
        M_momo <- length(id_r)
        BB <- lapply(momo,function(mm){mm$model$B})
        TT <- lapply(momo,function(mm){mm$model$t})
        VV <- lapply(momo,function(mm){mm$model$V})
        PP <- lapply(momo,function(mm){mm$model$P})
        CC <- lapply(momo,function(mm){mm$model$C})
        UU <- lapply(momo,function(mm){mm$model$U})
        UUSSTTAARR <- lapply(momo,function(mm){mm$model$U_star})
        TT_stat <- list(mean=matrix(0,n,r),var=matrix(0,n,r))
        VV_stat <- list(mean=matrix(0,q,r),var=matrix(0,q,r))
        CC_stat <- list(mean=matrix(0,q,r),var=matrix(0,q,r))
        PP_stat <- list(mean=matrix(0,p,r),var=matrix(0,p,r))
        UU_stat <- list(mean=matrix(0,p,r),var=matrix(0,p,r))
        UUSSTTAARR_stat <- list(mean=matrix(0,p,r),var=matrix(0,p,r))
        for(s in 1:r){
          for(i_m in 1:M_momo){
            sign_i_B <- sign(VV[[i_m]][which.max(abs(VV[[i_m]][,s])),s])
            if(sign_i_B<0){
              VV[[i_m]][,s] <- VV[[i_m]][,s]*(-1)
              PP[[i_m]][,s] <- PP[[i_m]][,s]*(-1)
              CC[[i_m]][,s] <- CC[[i_m]][,s]*(-1)
              UUSSTTAARR[[i_m]][,s] <- UUSSTTAARR[[i_m]][,s]*(-1)
              UU[[i_m]][,s] <- UU[[i_m]][,s]*(-1)
              TT[[i_m]][,s] <- TT[[i_m]][,s]*(-1)
            }
          }
        }
        for(i_m in 1:M_momo){
          TT_stat$mean <- TT_stat$mean + TT[[i_m]]/M_momo
          VV_stat$mean <- VV_stat$mean + VV[[i_m]]/M_momo
          PP_stat$mean <- PP_stat$mean + PP[[i_m]]/M_momo
          CC_stat$mean <- CC_stat$mean + CC[[i_m]]/M_momo
          UUSSTTAARR_stat$mean <- UUSSTTAARR_stat$mean + UUSSTTAARR[[i_m]]/M_momo
          UU_stat$mean <- UU_stat$mean + UU[[i_m]]/M_momo
        }
        for(i_m in 1:M_momo){
          TT_stat$var <-  TT_stat$var +
            (TT[[i_m]]-TT_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$t,function(mp){mp$var}))/
            M_momo
          VV_stat$var <-  VV_stat$var +
            (VV[[i_m]]-VV_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$V,function(mp){mp$var}))/
            M_momo
          PP_stat$var <-  PP_stat$var +
            (PP[[i_m]]-PP_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$P,function(mp){mp$var}))/
            M_momo
          CC_stat$var <-  CC_stat$var +
            (CC[[i_m]]-CC_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$C,function(mp){mp$var}))/
            M_momo
          UUSSTTAARR_stat$var <-  UUSSTTAARR_stat$var +
            (UUSSTTAARR[[i_m]]-UUSSTTAARR_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$U_star,function(mp){mp$var}))/
            M_momo
          UU_stat$var <-  UU_stat$var +
            (UU[[i_m]]-UU_stat$mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            do.call(cbind,lapply(momo[[i_m]]$results$U,function(mp){mp$var}))/
            M_momo
        }
        outr <- list(U=UU_stat,V=VV_stat,P=PP_stat,C=CC_stat,
                     U_star=UUSSTTAARR_stat,T=TT_stat)
        B_mean = B_var <- matrix(0,p,q)
        for(i_m in 1:M_momo){
          B_mean <-  B_mean + BB[[i_m]]/M_momo
        }
        for(i_m in 1:M_momo){
          B_var <-  B_var +
            (BB[[i_m]]-B_mean)^2*1/(M_momo-1)*(1+1/M_momo) +
            momo[[i_m]]$results$B[[r]]$var/M_momo
        }
        outr$B <- list(mean=B_mean,var=B_var)
      }else{
        outr <- NULL
      }
      out_total[[r]] <- outr
    }
  }
  out_total
}
