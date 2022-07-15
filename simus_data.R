
library(ddsPLS2)

n <- 20
sigma2 <- 0.1
p1 <- 10
p2 <- c(1,10,100)
p3 <- c(1,100,500)
n_test <- 1000
propsNA <- 0.45#c(0,0.1,0.3,0.45,0.6)
N_simus <- 100

doData <- function(n=100,p1,p2,p3,sigma2=0.5,pNA=0){
  p <- 2*p1+p2+p3
  q <- 3
  pPhi <- 3
  p_total <- p+pPhi+q
  Sigma <- diag(rep(1,p_total))
  mus <- rep(0,p_total)
  variables <- mvtnorm::rmvnorm(n,mus,Sigma)
  phis <- variables[,1:pPhi]
  noises <- variables[,-c(1:pPhi)]
  noiseY <- noises[,1:q]
  noiseX <- noises[,-c(1:q)]
  sdSignal <- sqrt(1-sigma2)
  sdNoise <- sqrt(sigma2)

  ## Y
  y <-cbind(sdSignal*phis[,1] + sdNoise*noiseY[,1],
            sdSignal*(phis[,1]+2*phis[,2])/sqrt(5) + sdNoise*noiseY[,2],
            noiseY[,3])

  ## X
  x <- matrix(NA,n,p)
  for(j in 1:p1){
    x[,j] <- sdSignal*phis[,1] + sdNoise*noiseX[,j]
    x[,p1+j] <- sdSignal*phis[,2] + sdNoise*noiseX[,j+p1]
  }
  for(j in 1:p2){
    x[,2*p1+j] <- sdSignal*phis[,3] + sdNoise*noiseX[,j+2*p1]
  }
  for(j in 1:p3){
    x[,2*p1+p2+j] <- noiseX[,j+2*p1+p2]
  }

  ## Remove data
  if(pNA>0){
    x[sample(1:prod(dim(x)),size = prod(dim(x))*pNA)]<-NA
  }

  ## Output
  list(x=x,y=y)
}

HDdoData <- function(n=100,p1,p2,p3,sigma2=0.5,pNA=0){
  p <- 2*p1+p2+p3
  q <- 3
  pPhi <- 3
  p_total <- p+pPhi+q
  # Sigma <- diag(rep(1,p_total))
  # mus <- rep(0,p_total)
  variables <- matrix(rnorm(n*p_total),nrow=n)
  phis <- variables[,1:pPhi]
  noises <- variables[,-c(1:pPhi)]
  noiseY <- noises[,1:q]
  noiseX <- noises[,-c(1:q)]
  sdSignal <- sqrt(1-sigma2)
  sdNoise <- sqrt(sigma2)

  ## Y
  y <-cbind(sdSignal*phis[,1] + sdNoise*noiseY[,1],
            sdSignal*(phis[,1]+2*phis[,2])/sqrt(5) + sdNoise*noiseY[,2],
            noiseY[,3])

  ## X
  x <- matrix(NA,n,p)
  for(j in 1:p1){
    x[,j] <- sdSignal*phis[,1] + sdNoise*noiseX[,j]
    x[,p1+j] <- sdSignal*phis[,2] + sdNoise*noiseX[,j+p1]
  }
  for(j in 1:p2){
    x[,2*p1+j] <- sdSignal*phis[,3] + sdNoise*noiseX[,j+2*p1]
  }
  for(j in 1:p3){
    x[,2*p1+p2+j] <- noiseX[,j+2*p1+p2]
  }

  ## Remove data
  if(pNA>0){
    x[sample(1:prod(dim(x)),size = prod(dim(x))*pNA)]<-NA
  }

  ## Output
  list(x=x,y=y)
}

doNameTest <- function(n,p1,p2,p3,sigma2,HD=F){
  out <- paste("../data/Simus/Koh_lanta/datas_test/",
               paste(c("n","p1","p2","p3","sigma2"),c(n,p1,p2,p3,sigma2),sep="_",collapse = "__"),
               ".RData",sep=""
  )
  if(HD){
    out <- paste("../data/Simus/Koh_lanta/datas_testHD/",
                 paste(c("n","p1","p2","p3","sigma2"),c(n,p1,p2,p3,sigma2),sep="_",collapse = "__"),
                 ".RData",sep=""
    )
  }
  out
}

doNameTrain <- function(n,p1,p2,p3,sigma2,pNA,id,HD=F){
  out <- paste("../data/Simus/Koh_lanta/datas/",
               paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
               ".RData",sep=""
  )
  if(HD){
    out <- paste("../data/Simus/Koh_lanta/datasHD/",
                 paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
                 ".RData",sep=""
    )
  }
  out
}

doNameRMSES <- function(n,p1,p2,p3,sigma2,pNA,id,HD=T){
  out <- paste("../data/Simus/Koh_lanta/results/RMSES__",
               paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
               ".csv",sep=""
  )
  if(HD){
    out <- paste("../data/Simus/Koh_lanta/resultsHD/RMSES__",
                 paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
                 ".csv",sep=""
    )
  }
  out
}

doNameMODELS <- function(n,p1,p2,p3,sigma2,pNA,id,HD=T){
  out <- paste("../data/Simus/Koh_lanta/results/MODELS__",
               paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
               ".RData",sep=""
  )
  if(HD){
    out <- paste("../data/Simus/Koh_lanta/resultsHD/MODELS__",
                 paste(c("n","p1","p2","p3","sigma2","pNA","id"),c(n,p1,p2,p3,sigma2,pNA,id),sep="_",collapse = "__"),
                 ".RData",sep=""
    )
  }
  out
}

generate_data <- function(){
  paras_test <- expand.grid(n_test,p1,p2,p3,sigma2)
  colnames(paras_test) <- c("n","p1","p2","p3","sigma1")
  for(i in 1:nrow(paras_test)){
    cat(i);cat(" ")
    pa <- as.numeric(paras_test[i,])
    set.seed(1)
    data_test <- doData(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=0)
    name <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5])
    save(data_test,file = name)
  }


  paras <- expand.grid(n,p1,p2,p3,sigma2,propsNA,1:N_simus)
  colnames(paras) <- c("n","p1","p2","p3","sigma1","pNA","id")
  for(i in 1:nrow(paras)){
    cat(i);cat("/");cat(nrow(paras));cat(" ")
    pa <- as.numeric(paras[i,])
    set.seed(pa[7])
    data_train <- doData(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6])
    name <- doNameTrain(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7])
    save(data_train,file = name)
  }
}

do_stuff <- function(NCORES=1){

  propsNA <- 0.6

  paras <- expand.grid(n,p1,p2,p3,sigma2,propsNA,1:N_simus)
  colnames(paras) <- c("n","p1","p2","p3","sigma1","pNA","id")

  B <- 20
  B_koh <- 100
  n_lambdas <- 10
  M <- 5
  R_max <- 7

  methods <- c("Koh-Lanta (in ddsPLS)",
               "MI-NIPALS-PLS",
               "NIPALS-PLS",
               "Mean-PLS",
               "missMDA-PLS",
               "Koh-Lanta (in ddsPLS LD)")


  n_B <- nrow(paras)
  NCORES_w <- min(NCORES,n_B)
  lB <- (n_B-n_B%%NCORES_w)/NCORES_w
  id_parallel <- sample(c(unlist(lapply(1:NCORES_w,function(ii){rep(ii,lB)})),1:(n_B%%NCORES_w)),n_B,replace = F)
  i_B=i <- 1
  cl <- makeCluster(NCORES_w)
  registerDoParallel(cl)
  RES <- foreach(i_B=1:NCORES_w,.packages = c("ddsPLS2","missMDA")) %do% {
    source("MI_PLS.R")
    source("kohLanta2_fun.R")
    source("simus_data.R")
    id_b <- which(id_parallel==i_B)
    paras_b <- paras[id_b,,drop=F]
    cat(i_B);cat(" ")
    for(i in 1:nrow(paras_b)){
      pa <- as.numeric(paras_b[i,])
      cat("\n  ");cat(id_b[i]);cat(" p3=");cat(pa[4]);cat(" pNA=");cat(pa[6])
      RMSES <- matrix(NA,6,3)
      MODELS <- list()
      nameRMSES <- doNameRMSES(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD = F)
      nameMODELS <- doNameMODELS(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD = F)

      # test <- F
      # if((file.exists(nameRMSES) & file.exists(nameMODELS))){
      #   test <- TRUE
      #   RMSES <- read.csv(nameRMSES)
      #   if(ncol(RMSES)==4){
      #     coco <- ncol(RMSES)
      #     RMSES <- RMSES[,(coco-2):coco]}
      #   load(nameMODELS)
      # }
      testP <- (pa[4]==1 &pa[3]==1)|(pa[4]==500 &pa[3]==100)
      if(testP)
      {
        test <- file.exists(nameRMSES) & file.exists(nameMODELS)
        if(test)
        {
          name <- doNameTrain(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7])
          load(file = name)
          nameTest <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5])
          load(file = nameTest)
          if(T)
          {
            ## Koh-lanta
            Mkh <- M
            if(pa[6]==0){Mkh <- 1}
            MODEL_KohLanta <- KohLanta(data_train$x,data_train$y,B = B_koh,
                                       init = "nipals",
                                       M = Mkh,verbose=F,n_lambdas = n_lambdas)
            nameTest <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5])
            x_test <- data_test$x
            est <- predict(MODEL_KohLanta,x_test)
            RMSES[1,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)
            MODELS[[1]] <- MODEL_KohLanta

            ## MI PLS
            models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=M,proper = T,errMin = errMin)
            MODELS[[2]] <- models
            est <- list(est=predict_MI_PLS(models,data_test$x))
            RMSES[2,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)

            ## I PLS
            models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=1,proper=F,errMin=errMin)
            MODELS[[3]] <- models
            est <- list(est=predict_MI_PLS(models,data_test$x))
            RMSES[3,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)

            ## MEAN PLS
            x_mean <- data_train$x
            mus_x <- colMeans(x_mean,na.rm = T)
            na_x <- which(is.na(mus_x))
            if(length(na_x)>0){
              mus_x[na_x] <- 0
            }
            idNA <- which(is.na(data_train$x),arr.ind = T)
            idNAii <- which(is.na(data_train$x),arr.ind = F)
            x_mean[idNAii] <- mus_x[idNA[,2]]
            models <- MI_PLS(x_mean,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
            MODELS[[4]] <- models
            est <- list(est=predict_MI_PLS(models,data_test$x))
            RMSES[4,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)

            ## missMDA
            if(pa[6]>0){
              idNAALL <- which(is.na(apply(data_train$x,2,sd,na.rm=T)))
              if(length(idNAALL)>0){
                ncomp <- estim_ncpPCA(data_train$x[,-idNAALL,drop=F], ncp.min = 0, ncp.max = R_max,verbose = T)$ncp
              }else{
                ncomp <- estim_ncpPCA(data_train$x, ncp.min = 0, ncp.max = R_max,verbose = T)$ncp
              }

              if (ncomp!=0){
                if(length(idNAALL)>0){
                  x_imp <- MIPCA(data_train$x[,-idNAALL,drop=F], ncp = ncomp, nboot=M,scale = TRUE,method = "Regularized")
                }else{
                  x_imp <- MIPCA(data_train$x, ncp = ncomp, nboot=M,scale = TRUE,method = "Regularized")
                }
                models <- MI_PLS(as.matrix(x_imp$res.MI[[1]]),data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
                if(length(idNAALL)>0){
                  y_preds <- predict_MI_PLS(models,data_test$x[,-idNAALL,drop=F])/M
                }else{
                  y_preds <- predict_MI_PLS(models,data_test$x)/M
                }
                for(m in 2:M){
                  models <- MI_PLS(as.matrix(x_imp$res.MI[[m]]),data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
                  if(length(idNAALL)>0){
                    y_preds <- y_preds + predict_MI_PLS(models,data_test$x[,-idNAALL,drop=F])/M
                  }else{
                    y_preds <- y_preds + predict_MI_PLS(models,data_test$x)/M
                  }
                }
              }else{
                mus <- colMeans(data_train$x,na.rm = T)
                idNAMM <- which(is.na(mus))
                if(length(idNAMM)>0){
                  mus[idNAMM] <- 0
                }
                x_imp <- data_train$x
                x_imp[which(is.na(x_imp))] <- mus[which(is.na(x_imp),arr.ind = T)[,2]]
                models <- MI_PLS(x_imp,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
                y_preds <- predict_MI_PLS(models,data_test$x)
              }
              RMSES[5,] <- sqrt(colMeans((y_preds-data_test$y)^2))/apply(data_test$y,2,sd)
              MODELS[[5]] <- list(imp=x_imp,models=models)
            }else{
              x_imp <- NULL
              models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
              y_preds <- predict_MI_PLS(models,data_test$x)
              RMSES[5,] <- sqrt(colMeans((y_preds-data_test$y)^2))/apply(data_test$y,2,sd)
              MODELS[[5]] <- list(imp=x_imp,models=models)
            }
          }
          ## Koh-lanta LD
          Mkh <- M
          if(pa[6]==0){Mkh <- 1}
          MODEL_KohLantaLD <- KohLanta(data_train$x,data_train$y,B = B_koh,LD = T,
                                       init = "nipals",
                                       M = Mkh,verbose=F,n_lambdas = n_lambdas)
          est <- predict(MODEL_KohLantaLD,data_test$x)
          RMSES[6,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)
          MODELS[[6]] <- MODEL_KohLantaLD

          write.csv(RMSES,file=nameRMSES)
          save(MODELS,file=nameMODELS)
        }
        else{
          RMSES <- read.csv(nameRMSES)[,-1]
          if(nrow(RMSES)==5)
          {
            cat("    Je fais LD   \n")
            ## Koh-lanta LD
            Mkh <- M
            if(pa[6]==0){Mkh <- 1}
            MODEL_KohLantaLD <- KohLanta(data_train$x,data_train$y,B = B_koh,LD = T,
                                         M = Mkh,verbose=F,n_lambdas = n_lambdas)
            est <- predict(MODEL_KohLantaLD,data_test$x)
            RMSES <- rbind(RMSES,sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd))
            MODELS[[6]] <- MODEL_KohLantaLD

            write.csv(RMSES,file=nameRMSES)
            save(MODELS,file=nameMODELS)
          }
        }
      }
    }
  }
  stopCluster(cl)
}

HDgenerate_data <- function(){
  paras_test <- expand.grid(n_test,p1,p2,p3,sigma2)
  colnames(paras_test) <- c("n","p1","p2","p3","sigma1")
  for(i in 1:nrow(paras_test)){
    cat(i);cat(" ")
    pa <- as.numeric(paras_test[i,])
    set.seed(1)
    data_test <- HDdoData(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=0)
    name <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],HD=T)
    save(data_test,file = name)
  }


  paras <- expand.grid(n,p1,p2,p3,sigma2,propsNA,1:N_simus)
  colnames(paras) <- c("n","p1","p2","p3","sigma1","pNA","id")
  for(i in 1:nrow(paras)){
    cat(i);cat("/");cat(nrow(paras));cat(" ")
    pa <- as.numeric(paras[i,])
    set.seed(pa[7])
    data_train <- HDdoData(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6])
    name <- doNameTrain(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD=T)
    save(data_train,file = name)
  }
}

HDdo_stuff <- function(NCORES=1){


  library(ddsPLS2)

  n <- 100
  sigma2 <- 0.1
  p1 <- 10
  p2 <- 100
  p3 <- 20000
  n_test <- 1000
  propsNA <- c(0,0.1,0.3,0.6)
  N_simus <- 100

  paras <- expand.grid(n,p1,p2,p3,sigma2,propsNA,1:N_simus)
  colnames(paras) <- c("n","p1","p2","p3","sigma1","pNA","id")

  B <- 20
  B_koh <- 20
  n_lambdas <- 10
  M <- 5
  R_max <- 7

  methods <- c("Koh-Lanta (in ddsPLS)",
               "MI-NIPALS-PLS",
               "NIPALS-PLS",
               "Mean-PLS",
               "missMDA-PLS",
               "Koh-Lanta (in ddsPLS LD)")


  n_B <- nrow(paras)
  NCORES_w <- min(NCORES,n_B)
  lB <- (n_B-n_B%%NCORES_w)/NCORES_w
  id_parallel <- sample(c(unlist(lapply(1:NCORES_w,function(ii){rep(ii,lB)})),1:(n_B%%NCORES_w)),n_B,replace = F)
  i_B=i <- 1
  cl <- makeCluster(NCORES_w)
  registerDoParallel(cl)
  RES <- foreach(i_B=1:NCORES_w,.packages = c("ddsPLS2","missMDA")) %do% {
    source("MI_PLS.R")
    source("kohLanta2_fun.R")
    source("simus_data.R")
    id_b <- which(id_parallel==i_B)
    paras_b <- paras[id_b,,drop=F]
    cat(i_B);cat(" ")
    for(i in 1:nrow(paras_b)){
      pa <- as.numeric(paras_b[i,])
      cat("\n  ");cat(id_b[i]);cat(" p3=");cat(pa[4]);cat(" pNA=");cat(pa[6])
      RMSES <- matrix(NA,6,3)
      MODELS <- list()
      nameRMSES <- doNameRMSES(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD=T)
      nameMODELS <- doNameMODELS(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD=T)

      # test <- F
      # if((file.exists(nameRMSES) & file.exists(nameMODELS))){
      #   test <- TRUE
      #   RMSES <- read.csv(nameRMSES)
      #   if(ncol(RMSES)==4){
      #     coco <- ncol(RMSES)
      #     RMSES <- RMSES[,(coco-2):coco]}
      #   load(nameMODELS)
      # }

      test <- (file.exists(nameRMSES) & file.exists(nameMODELS))

      if(!test){
        name <- doNameTrain(n=pa[1],p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],pNA=pa[6],id = pa[7],HD=T)
        load(file = name)
        nameTest <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5],HD=T)
        load(file = nameTest)

        ## Koh-lanta
        Mkh <- M
        if(pa[6]==0){Mkh <- 1}
        MODEL_KohLanta <- KohLanta(data_train$x,data_train$y,B = B_koh,
                                   M = Mkh,verbose=T,n_lambdas = n_lambdas)
        nameTest <- doNameTest(n=n_test,p1=pa[2],p2=pa[3],p3=pa[4],sigma2=pa[5])
        x_test <- data_test$x
        est <- predict(MODEL_KohLanta,x_test)
        RMSES[1,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)
        MODELS[[1]] <- MODEL_KohLanta

        ## MI PLS
        models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=M,proper = T,errMin = errMin)
        MODELS[[2]] <- models
        est <- list(est=predict_MI_PLS(models,data_test$x))
        RMSES[2,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)


        ## I PLS
        models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=1,proper=F,errMin=errMin)
        MODELS[[3]] <- models
        est <- list(est=predict_MI_PLS(models,data_test$x))
        RMSES[3,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)


        ## MEAN PLS
        x_mean <- data_train$x
        mus_x <- colMeans(x_mean,na.rm = T)
        na_x <- which(is.na(mus_x))
        if(length(na_x)>0){
          mus_x[na_x] <- 0
        }
        idNA <- which(is.na(data_train$x),arr.ind = T)
        idNAii <- which(is.na(data_train$x),arr.ind = F)
        x_mean[idNAii] <- mus_x[idNA[,2]]
        models <- MI_PLS(x_mean,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
        MODELS[[4]] <- models
        est <- list(est=predict_MI_PLS(models,data_test$x))
        RMSES[4,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)

        ## missMDA
        if(pa[6]>0){
          idNAALL <- which(is.na(apply(data_train$x,2,sd,na.rm=T)))
          if(length(idNAALL)>0){
            ncomp <- estim_ncpPCA(data_train$x[,-idNAALL,drop=F], ncp.min = 0, ncp.max = R_max,verbose = T)$ncp
          }else{
            ncomp <- estim_ncpPCA(data_train$x, ncp.min = 0, ncp.max = R_max,verbose = T)$ncp
          }

          if (ncomp!=0){
            if(length(idNAALL)>0){
              x_imp <- MIPCA(data_train$x[,-idNAALL,drop=F], ncp = ncomp, nboot=M,scale = TRUE,method = "Regularized")
            }else{
              x_imp <- MIPCA(data_train$x, ncp = ncomp, nboot=M,scale = TRUE,method = "Regularized")
            }
            models <- MI_PLS(as.matrix(x_imp$res.MI[[1]]),data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
            if(length(idNAALL)>0){
              y_preds <- predict_MI_PLS(models,data_test$x[,-idNAALL,drop=F])/M
            }else{
              y_preds <- predict_MI_PLS(models,data_test$x)/M
            }
            for(m in 2:M){
              models <- MI_PLS(as.matrix(x_imp$res.MI[[m]]),data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
              if(length(idNAALL)>0){
                y_preds <- y_preds + predict_MI_PLS(models,data_test$x[,-idNAALL,drop=F])/M
              }else{
                y_preds <- y_preds + predict_MI_PLS(models,data_test$x)/M
              }
            }
          }else{
            mus <- colMeans(data_train$x,na.rm = T)
            idNAMM <- which(is.na(mus))
            if(length(idNAMM)>0){
              mus[idNAMM] <- 0
            }
            x_imp <- data_train$x
            x_imp[which(is.na(x_imp))] <- mus[which(is.na(x_imp),arr.ind = T)[,2]]
            models <- MI_PLS(x_imp,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
            y_preds <- predict_MI_PLS(models,data_test$x)
          }
          RMSES[5,] <- sqrt(colMeans((y_preds-data_test$y)^2))/apply(data_test$y,2,sd)
          MODELS[[5]] <- list(imp=x_imp,models=models)
        }else{
          x_imp <- NULL
          models <- MI_PLS(data_train$x,data_train$y,R=R_max,M=1,proper = F,errMin = errMin)
          y_preds <- predict_MI_PLS(models,data_test$x)
          RMSES[5,] <- sqrt(colMeans((y_preds-data_test$y)^2))/apply(data_test$y,2,sd)
          MODELS[[5]] <- list(imp=x_imp,models=models)
        }

        ## Koh-lanta LD
        Mkh <- M
        if(pa[6]==0){Mkh <- 1}
        MODEL_KohLantaLD <- KohLanta(data_train$x,data_train$y,B = B_koh,LD = T,
                                     M = Mkh,verbose=F,n_lambdas = n_lambdas)
        est <- predict(MODEL_KohLantaLD,data_test$x)
        RMSES[6,] <- sqrt(colMeans((est$est-data_test$y)^2))/apply(data_test$y,2,sd)
        MODELS[[6]] <- MODEL_KohLantaLD

        print(RMSES)

        write.csv(RMSES,file=nameRMSES)
        save(MODELS,file=nameMODELS)
      }

    }
  }
  stopCluster(cl)
}

