
NNN <- 6

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


cols <- RColorBrewer::brewer.pal(8,"Set1")
cols[6] <- "brown"
angle <- c(-45,45,90,-35,0,10)
density <- c(30,35,40,45,50,25)
nM <- length(methods)
nn <- length(unique(df$pNA))
pNAs <- sort(as.numeric(unique(df$pNA)))

p2s <- c(1,100)
p3s <- c(1,500)
ys <- unique(df$y)

parass <- expand.grid(p2s,p3s)[c(1,4),]

df$method <- factor(df$method,levels=methodsPLOT[c(1,6,2,3,4,5)])




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
  if(!is.null(MODELS[[1]]))
  {
    RMDA <- MODELS[[5]]$models[[1]]$R
    out <- rbind(cbind("Koh-Lanta (in ddsPLS)",unlist(lapply(MODELS[[1]]$models,function(m){m$R}))),
                 cbind("MI-NIPALS",unlist(lapply(MODELS[[2]],function(m){m$R})) ),
                 cbind("NIPALS",MODELS[[3]][[1]]$R),
                 cbind("Mean-PLS",0))
    if(!is.null(RMDA)){
      out <- rbind(out,cbind("missMDA-PLS",RMDA))
    }
    if(length(MODELS)>5 ){
      out <- rbind(out,cbind("Koh-Lanta (in ddsPLS LD)",unlist(lapply(MODELS[[6]]$models,function(m){m$R}))))
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
  }
  else
  {
    # browser()
    if(length(MODELS)>5 ){
      outou <- cbind("Koh-Lanta (in ddsPLS LD)",unlist(lapply(MODELS[[6]]$models,function(m){m$R})))
      out <- outou
      outPara <- matrix(rep(para_i,nrow(out)),byrow = T,ncol = length(para_i))
      OUTR <- rbind(OUTR,cbind(outPara,out))
      # if(para_i[3]==100 & para_i[4]==500 & para_i[6]==0.6 ) browser()
    }
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

dfR <- data.frame(OUTR)
names(dfR) <- c("n","p1","p2","p3","sigma1","pNA","method","R")
dfR$R <- as.numeric(dfR$R)

dfR$method <- factor(dfR$method,levels=methodsPLOT)



for(j in 1:3){
  yy <- ys[j]
  postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__Plot_y",j,"_CFEM.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height = 9, width = 12)
  par(mar=c(4,4,3,2))
  KK <- 8
  mama <- rbind(rbind(matrix(1,NNN,NNN),matrix(2,NNN,NNN)),3)
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
                 main=paste("p2=",p2i,", p3=",p3i,sep=""))
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
  for(i in 1){
    par(mar=c(0,0,0,0))
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n",col="white")
    df_id <- df[id,]
    methoLegend <- levels(df$method)
    legend("center",legend = methoLegend,angle=angle,
           border=legend_ok$col,fill=legend_ok$col,density=density[1:nM],
           col = legend_ok$col,pch=15+(1:nM),
           ncol = 6,bty = "n",pt.cex=0.2,horiz = F,cex=1.2)
  }
  dev.off()
}



postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__R_CFEM.eps",sep=""),
           horizontal = F, onefile = FALSE, paper = "special",
           height = 9, width = 15)
mama <- rbind(rbind(matrix(1,NNN,NNN),matrix(2,NNN,NNN)),3)

dfR_mem <- dfR
dfR$method <- factor(dfR$method,levels=levels(df$method))

layout(mama)
matricesRCORREC <- list()
for(i in 1:nrow(parass)){
  par(mar=c(3,1,2,0))
  cat("\n ---------------")
  print(parass[i,]);cat("\n")
  p2i <- parass[i,1]
  p3i <- parass[i,2]
  id <- which(dfR$p3==p3i & dfR$p2==p2i)
  print(table(dfR$method[id],dfR$pNA[id]))
  coucou <- table(dfR$R[id],dfR$pNA[id],dfR$method[id])
  matricesRCORREC[[i]] <- list()
  matricesRCORREC[[i]][[1]] =matricesRCORREC[[i]][[2]] = matrix(NA,nM,nn)
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
        if(iR==3|iR==2)
        {
          matricesRCORREC[[i]][[iR-1]][im,ipna] <- coucoui/denom
        }
      }
    }
  }
  # matricesRCORREC[[i]] <- maR2
  mou <- matrix(mi$freq,ncol=1)
  par(cex=1)
  plot(0,0,xlim=c(1,length(mou)),ylim=c(0,1.05),col="white",
       xaxt="n",yaxt="n",bty="n",xlab="",ylab="",yaxt="n",
       main=paste("p2=",p2i,", p3=",p3i,sep=""))
  for(ipna in 1:nn){
    x0 <- (ipna-1)*(nM)*(nR) + 2*nM +1
    x1 <- x0+nM
    rect(xleft = x0,xright = x1,ybottom = 0,
         ytop = 1,col = "gray",border = "white")
  }
  barplot(mou,col=mi$col,border = mi$col,angle=angle,density=density[1:nM],beside = T,add=T,yaxt="n")
  abline(v=(1:(nn-1))*(nM)*(nR)+1,lwd=2)
  abline(v=(1:(nn*nR))*nM+1,lty=2)
  par(cex=1)
  axis(2,at = (0:10)/10,(0:10)/10,las=2,tick = F,line = -1.5)
  par(cex=1)
  axis(1,line = -1,at = (1:(nn*nR))*nM-nM/4,labels = paste(rep(0:(nR-1),nn),sep=""),tick = F)
  par(cex=1)
  axis(1,line = 0.5,at = (0:(nn-1)+1/2)*(nM)*(nR)+5,labels = paste("pNA=",pNAs,sep=""),tick = F)
  rect(xleft = 0,xright = 1e4,ybottom = 1,
       ytop = 2,col = "white",border = "white")
}
par(mar=c(0,0,0,0))
plot(0,0,xlim=c(1,length(mou)),ylim=c(0,1.05),col="white",
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",yaxt="n")
legend_ok <- unique(mi[,c("method","col")])
legend("center",legend = (legend_ok$method),angle=angle,
       border=legend_ok$col,fill=legend_ok$col,density=density[1:nM],
       col = legend_ok$col,pch=15+(1:nM),
       ncol = 6,bty = "n",pt.cex=0.2,horiz = F,cex=1.0)
dev.off()






postscript(paste("../data/Simus/Koh_lanta/",nPlo,"__R_CFEM_light.eps",sep=""),
           horizontal = F, onefile = FALSE, paper = "special",
           height = 6, width = 10)
m <- matrix(1,4,4)
layout(rbind(rbind(cbind(m,2*m),cbind(3*m,4*m)),rep(5,8)))
par(mar=c(3,3,3,0))
barplot(matricesRCORREC[[1]][[1]]*100,beside=T,ylim=c(0,100),
        col=rep(cols[1:nM],nn ),main="p1=1, p3=1\n R=1")
abline(v=(nM+1)*(1:(nn-1))+0.5,lty=2,lwd=2)
axis(1,line = 0,at = (0:(nn-1))*(nM+1)+4,labels = paste("pNA=",pNAs,sep=""),tick = F)
barplot(matricesRCORREC[[1]][[2]]*100,beside=T,ylim=c(0,100),
        col=rep(cols[1:nM],nn ),main="p1=1, p3=1\n R=2")
abline(v=(nM+1)*(1:(nn-1))+0.5,lty=2,lwd=2)
axis(1,line = 0,at = (0:(nn-1))*(nM+1)+4,labels = paste("pNA=",pNAs,sep=""),tick = F)

barplot(matricesRCORREC[[2]][[1]]*100,beside=T,ylim=c(0,100),
        col=rep(cols[1:nM],nn ),main="p1=100, p3=500\n R=1")
abline(v=(nM+1)*(1:(nn-1))+0.5,lty=2,lwd=2)
axis(1,line = 0,at = (0:(nn-1))*(nM+1)+4,labels = paste("pNA=",pNAs,sep=""),tick = F)
barplot(matricesRCORREC[[2]][[2]]*100,beside=T,ylim=c(0,100),
        col=rep(cols[1:nM],nn ),main="p1=100, p3=500\n R=2")
abline(v=(nM+1)*(1:(nn-1))+0.5,lty=2,lwd=2)
axis(1,line = 0,at = (0:(nn-1))*(nM+1)+4,labels = paste("pNA=",pNAs,sep=""),tick = F)

par(mar=c(0,0,0,0))
plot(0,0,xlim=c(1,length(mou)),ylim=c(0,1.05),col="white",
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",yaxt="n")
legend_ok <- unique(mi[,c("method","col")])
legend("center",legend = (legend_ok$method),
       border=legend_ok$col,fill=legend_ok$col,
       col = legend_ok$col,pch=15+(1:nM),
       ncol = 6,bty = "n",pt.cex=0.2,horiz = F,cex=1.0)
dev.off()
