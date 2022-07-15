getTPRFPR <- function(vec,obj,M,all){
  DIM <- length(all)
  predits <- vec
  PASpredits <- all[which(!(all %in% vec))]
  apredire <- obj
  aPASpredire <- all[which(!(all %in% obj))]
  TP <- length(intersect(predits,apredire))/M
  FN <- length(intersect(PASpredits,apredire))/M
  FP <- length(intersect(predits,aPASpredire))/M
  TN <- length(intersect(PASpredits,aPASpredire))/M
  c(TP/(TP+FN),TN/(TN+FP))
}

fifis <- list.files("../data/Simus/Koh_lanta/results/")
splispli <- strsplit(fifis,split = "__",fixed = T)
idRMSE <- which(unlist(lapply(splispli,function(fi){fi[1]=="RMSES"})))
idMo <- which(unlist(lapply(splispli,function(fi){fi[1]=="MODELS"})))

nPlo <- 100
isnPlo <- unlist(lapply(splispli,function(fi){fi[2]==paste("n",nPlo,sep="_")}))

idRMSE <- idRMSE[which(idRMSE %in% which(isnPlo))]
idMo <- idMo[which(idMo %in% which(isnPlo))]


paras <- matrix(NA,7)
methods <- c("Koh-Lanta (in ddsPLS)",
             "MI-NIPALS",
             "NIPALS",
             "Mean-PLS",
             "missMDA-PLS",
             "Koh-Lanta (in ddsPLS LD)")
methodsPLOT <- c("Koh-Lanta (in ddsPLS)",
                 "MI-NIPALS",
                 "NIPALS",
                 "Mean-PLS",
                 "missMDA-PLS",
                 "Koh-Lanta (in ddsPLS LD)")
OUT <- NULL
co <- NULL
for(i in 1:length(idRMSE)){
  para_i <- as.numeric(unlist(lapply(splispli[[idRMSE[i]]][-c(1,8)],
                                     function(ii){
                                       unlist(strsplit(
                                         ii,split = "_",fixed=T))[2]})))
  err <- read.csv(paste("../data/Simus/Koh_lanta/results/",
                        fifis[idRMSE[i]],sep=""))[,-1]
  idMeth <- 1:nrow(err)
  out <- rbind(cbind(methodsPLOT[idMeth],"y 1",err[,1]),
               cbind(methodsPLOT[idMeth],"y 2",err[,2]),
               cbind(methodsPLOT[idMeth],"y 3",err[,3]))
  outPara <- matrix(rep(para_i,nrow(out)),byrow = T,ncol = length(para_i))
  OUT <- rbind(OUT,cbind(outPara,out))
  # if(para_i[2]==10 & para_i[3]==1 & para_i[4]==1 & para_i[6]==0 & nrow(err)==6){
  #   if(err[6,1]>0.9) browser()
  # }
}

df <- data.frame(OUT)
names(df) <- c("n","p1","p2","p3","sigma1","pNA","method","y","rmse")
df$rmse <- as.numeric(df$rmse)

p2s <- unique(df$p2)
p3s <- unique(df$p3)
ys <- unique(df$y)

parass <- expand.grid(p2s,p3s)

cols <- RColorBrewer::brewer.pal(8,"Set1")
cols[6] <- "brown"
angle <- c(-45,45,90,-35,0,10)
density <- c(30,35,40,45,50,25)
nM <- length(methods)
nn <- length(unique(df$pNA))
pNAs <- sort(as.numeric(unique(df$pNA)))



df$method <- factor(df$method,levels=methodsPLOT[c(1,6,2,3,4,5)])

for(j in 1:3){
  yy <- ys[j]
  # pdf(file = paste("../data/Simus/Koh_lanta/",nPlo,"__Plot_",yy,".pdf",sep=""),   # The directory you want to save the file in
  #     width = 12, # The width of the plot in inches
  #     height = 12)
  postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__Plot_y",j,".eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height = 12, width = 12)
  par(mar=c(4,4,3,2))
  KK <- 8
  mama <- do.call(rbind,lapply(0:2*3,function(j){
    out <- matrix(rep(c(rep(j+1,KK),rep(j+2,KK),rep(j+3,KK)),KK),nrow = KK,byrow = T)
    if(j==0) out <- rbind(out,10)
    if(j==3) out <- rbind(out,11)
    out
  }))
  layout(mat = mama)
  for(i in 1:nrow(parass)){
    cat("\n ---------------")
    print(parass[i,]);cat("\n")
    p2i <- parass[i,1]
    p3i <- parass[i,2]
    id <- which(df$p3==p3i &df$p2==p2i &df$y==yy)
    print(table(df$method[id],df$pNA[id]))
    b <- boxplot(rmse~method*pNA,df[id,],las=2,
                 col="white",density=density[1:nM],angle=angle,xaxt="n",
                 xlab="",ylab="",ylim=c(0.3,1.3),
                 main=paste("n=",nPlo,", p1=",10,", p2=",p2i,", p3=",p3i,", ",yy,sep=""))
    axis(1,at = -3/2+(1:nn)*nM,labels = pNAs)
    title(xlab=expression(p["NA"]),line = 2.5)
    title(ylab="RMSE",line = 2.5)
    abline(v=(1:4)*nM+1/2,lty=3)
    abline(h=c(sqrt(0.1),1),lty=2,col=2)
    id_meth <- rep(1:nM,nn)
    coli_s <- cols[id_meth]
    for(ii in 1:ncol(b$stats)){
      rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
           ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           col=coli_s[ii])
    }
  }
  for(i in 1:2){
    par(mar=c(0,0,0,0))
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n",col="white")
    legend("center",legend = levels((df$method)),angle=angle,
           border=cols[1:nM],fill=cols[1:nM],density=density[1:nM],
           col = cols[1:nM],pch=15+(1:nM),
           ncol = 1,bty = "n",pt.cex=0.2,horiz = T)
  }
  dev.off()
}



OUTR <- NULL
SEL <- list()
Ms <- matrix(100,length(idMo),2)
PARAS <- matrix(NA,length(idMo),6)
TESTS <- rep(NA,length(idMo))
for(i in 1:length(idMo)){
  cat(i);cat(" ")
  para_i <- as.numeric(unlist(lapply(splispli[[idMo[i]]][-c(1,8)],function(ii){unlist(strsplit(ii,split = "_",fixed=T))[2]})))
  nameModel <- paste("../data/Simus/Koh_lanta/results/",fifis[idMo[i]],sep="")
  load(nameModel)
  TESTS[i] <- file.exists(nameModel)
  if(length(MODELS[[1]])>0){
    RMDA <- MODELS[[5]]$models[[1]]$R
    out <- rbind(cbind("Koh-Lanta (in ddsPLS)",unlist(lapply(MODELS[[1]]$models,function(m){m$R}))),
                 cbind("MI-NIPALS",unlist(lapply(MODELS[[2]],function(m){m$R})) ),
                 cbind("NIPALS",MODELS[[3]][[1]]$R),
                 cbind("Mean-PLS",0))
    if(!is.null(RMDA)){
      out <- rbind(out,cbind("missMDA-PLS",RMDA))
    }
    if(length(MODELS)>5 ){
      out <- rbind(out,cbind("Koh-Lanta (in ddsPLS LD)",MODELS[[6]]$models[[1]]$R))
    }
    outPara <- matrix(rep(para_i,nrow(out)),byrow = T,ncol = length(para_i))
    OUTR <- rbind(OUTR,cbind(outPara,out))

    ## Selection
    Mi <- length(MODELS[[1]]$models)
    selY_Koh_lanta <- unlist(lapply(1:Mi,function(m){MODELS[[1]]$models[[m]]$Selection$Y}))
    selX_Koh_lanta <- unlist(lapply(1:Mi,function(m){MODELS[[1]]$models[[m]]$Selection$X}))
    if(length(MODELS)>5){
      Mj <- length(MODELS[[6]]$models)
      selY_Koh_lanta_LD <- unlist(lapply(1:Mj,function(m){MODELS[[6]]$models[[m]]$Selection$Y}))
      selX_Koh_lanta_LD <- unlist(lapply(1:Mj,function(m){MODELS[[6]]$models[[m]]$Selection$X}))
    }else{
      Mj <- NA
      selY_Koh_lanta_LD=selX_Koh_lanta_LD <- NULL
    }
    Ms[i,] <- c(Mi,Mj)
    PARAS[i,] <- para_i
    SEL[[i]] <- list(HD=list(x=selX_Koh_lanta,y=selY_Koh_lanta),
                     LD=list(x=selX_Koh_lanta_LD,y=selY_Koh_lanta_LD))
  }else{
    Mi <- NA
    selY_Koh_lanta=selX_Koh_lanta <- NULL
    if(length(MODELS)>5){
      Mj <- length(MODELS[[6]]$models)
      selY_Koh_lanta_LD <- unlist(lapply(1:Mj,function(m){MODELS[[6]]$models[[m]]$Selection$Y}))
      selX_Koh_lanta_LD <- unlist(lapply(1:Mj,function(m){MODELS[[6]]$models[[m]]$Selection$X}))
    }else{
      Mj <- NA
      selY_Koh_lanta_LD=selX_Koh_lanta_LD <- NULL
    }
    Ms[i,] <- c(Mi,Mj)
    PARAS[i,] <- para_i
    SEL[[i]] <- list(HD=list(x=selX_Koh_lanta,y=selY_Koh_lanta),
                     LD=list(x=selX_Koh_lanta_LD,y=selY_Koh_lanta_LD))
  }
}


sel_obj_X <- 1:20
sel_obj_Y <- 1:2



TPR_FPR <- lapply(1:min(length(SEL),nrow(Ms)),function(i){
  cat(i);cat(" ")
  si <- SEL[[i]]
  Mi <- Ms[i,]
  out <- si
  All <- 1:(2*PARAS[i,2]+sum(PARAS[i,3:4]))
  out$HD$x <- getTPRFPR(si$HD$x,sel_obj_X,Mi[1],All)
  out$HD$y <- getTPRFPR(si$HD$y,sel_obj_Y,Mi[1],All)
  out$LD$x <- getTPRFPR(si$LD$x,sel_obj_X,Mi[2],All)
  out$LD$y <- getTPRFPR(si$LD$y,sel_obj_Y,Mi[2],All)
  out
})

parass <- expand.grid(as.numeric(p2s),as.numeric(p3s))
postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__SELX.eps",sep=""),
           horizontal = F, onefile = FALSE, paper = "special",
           height = 8, width = 8)
par(mfrow=c(3,3))
par(mar=c(4,4,2,0))
for(i in 1:nrow(parass)){
  pai <- as.numeric(parass[i,])
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",
       xlab="",ylab="",bty="n",
       main=paste("n=",nPlo,", p2=",pai[1],", p3=",pai[2],sep=""))
  X1 = X2=Y1=Y2 <- NULL
  title(xlab = "FPR",ylab="TPR",line = 2)
  for(j in 1:length(pNAs)){
    pnaj <- pNAs[j]
    id <- which(PARAS[,3]==pai[1] & PARAS[,4]==pai[2] & PARAS[,6]==pnaj)
    distris <- TPR_FPR[id]
    HDX <- colMeans(do.call(rbind,lapply(distris,function(dd){dd$HD$x})),na.rm = T)
    HDY <- colMeans(do.call(rbind,lapply(distris,function(dd){dd$HD$y})),na.rm = T)
    LDX <- colMeans(do.call(rbind,lapply(distris,function(dd){dd$LD$x})),na.rm = T)
    LDY <- colMeans(do.call(rbind,lapply(distris,function(dd){dd$LD$y})),na.rm = T)
    X1 <-c(X1,1-HDX[2])
    X2 <-c(X2,1-LDX[2])
    Y1 <-c(Y1,HDX[1])
    Y2 <-c(Y2,LDX[1])
    text(x = 1-HDX[2],y = HDX[1],col=cols[1],labels =as.character(pNAs[j]))
    text(1-LDX[2],LDX[1],col=cols[2],labels =as.character(pNAs[j]))
  }
  points(X1,Y1,type="l",lty=1,col=cols[1])
  points(X2,Y2,type="l",lty=2,col=cols[2])
  legend("right",bty="n",
         legend = c(methods[c(1,6)]),
         col=c(cols[c(1,2)]),lty=1:2)
}
dev.off()

dfR <- data.frame(OUTR)
names(dfR) <- c("n","p1","p2","p3","sigma1","pNA","method","R")
dfR$R <- as.numeric(dfR$R)

dfR$method <- factor(dfR$method,levels=methodsPLOT)


# pdf(file = paste("../data/Simus/Koh_lanta/",nPlo,"__R.pdf",sep=""),   # The directory you want to save the file in
#     width = 30, # The width of the plot in inches
#     height = 10)
postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__R.eps",sep=""),
           horizontal = F, onefile = FALSE, paper = "special",
           height = 12, width = 30)
par(mar=c(3,1,2,0))
# KK <- 2
# mama2 <- lapply(1:9,function(j){
#   out <- rbind(do.call(cbind,lapply(1:3,function(ii){matrix(1+ii+(j-1)*6,KK,KK)})),
#                do.call(cbind,lapply(1:3+3,function(ii){matrix(1+ii+(j-1)*6,KK,KK)})))
#   if(j %in% c(1,2,4,5,7,8)){
#     out <- cbind(rbind(out,1),1)
#   }else{
#     out <- rbind(out,1)
#   }
# })
# mama_2 <- rbind(cbind(mama2[[1]],mama2[[2]],mama2[[3]]),
#                 cbind(mama2[[4]],mama2[[5]],mama2[[6]]),
#                 cbind(mama2[[7]],mama2[[8]],mama2[[9]]))
mama_2 <- matrix(1:nrow(parass),nrow = 3,byrow = T)
layout(mat = mama_2)
# plot(1,1,col="white",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
for(i in 1:nrow(parass)){
  cat("\n ---------------")
  print(parass[i,]);cat("\n")
  p2i <- parass[i,1]
  p3i <- parass[i,2]
  id <- which(dfR$p3==p3i & dfR$p2==p2i)
  print(table(dfR$method[id],dfR$pNA[id]))
  coucou <- table(dfR$R[id],dfR$pNA[id],dfR$method[id])
  mi <- NULL
  nR <- nrow(coucou[,,1])
  methodoLev <- levels(dfR$method[id])
  for(ipna in 1:nn){
    for(iR in 1:nR){
      for(im in 1:nM){
        coucoui <- coucou[iR,ipna,im]
        denom <- sum(coucou[,ipna,im])
        mi <- rbind(mi,data.frame(list(freq=coucoui/denom,
                                       R=iR-1,
                                       pNA=pNAs[ipna],
                                       method=methodoLev[im],
                                       col=cols[im])))
      }
    }
  }
  mou <- matrix(mi$freq,ncol=1)
  par(cex=1)
  plot(0,0,xlim=c(1,length(mou)),ylim=c(0,1.2),col="white",
       xaxt="n",yaxt="n",bty="n",xlab="",ylab="",yaxt="n",
       main=paste("n=",nPlo,", p2=",p2i,", p3=",p3i,sep=""))
  for(ipna in 1:nn){
    x0 <- (ipna-1)*(nM)*(nR) + 2*nM +1
    x1 <- x0+nM
    rect(xleft = x0,xright = x1,ybottom = 0,
         ytop = 1,col = "gray",border = "white")
  }
  barplot(mou,col=mi$col,border = mi$col,angle=angle,density=density[1:nM],beside = T,add=T,yaxt="n")
  abline(v=(1:(nn-1))*(nM)*(nR)+1,lwd=2)
  abline(v=(1:(nn*nR))*nM+1,lty=2)
  par(cex=3/4)
  axis(2,at = (0:10)/10,(0:10)/10,las=2,tick = F,line = -1.5)
  par(cex=1/2)
  axis(1,line = -1,at = (1:(nn*nR))*nM-nM/4,labels = paste("R=",rep(0:(nR-1),nn),sep=""),tick = F)
  par(cex=3/4)
  axis(1,line = 0.5,at = (0:(nn-1)+1/2)*(nM)*(nR)+5,labels = paste("pNA=",pNAs,sep=""),tick = F)
  rect(xleft = 0,xright = 1e4,ybottom = 1,
       ytop = 2,col = "white",border = "white")
  legend("top",legend = levels(df$method),angle=angle,
         border=cols[1:nM],fill=cols[1:nM],density=density[1:nM],
         col = cols[1:nM],pch=15+(1:nM),
         ncol = 3,bty = "n",pt.cex=0.2,horiz = F)
}

dev.off()
