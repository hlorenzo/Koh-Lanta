plot_1 <- function(){
  ddff <- all_1_df
  ddff$method <- factor(ddff$method,levels = c("ddsPLS","MI-NIPALS","NIPALS","MEAN-NIPALS"))
  MM <- unique(ddff$method)
  nM <- length(MM)
  nn <- length(unique(ddff$n))
  NN <- 6
  layout(rbind(matrix(rep(1:3,NN),nrow = NN,byrow = T),c(1,1,1)*4))
  par(mar=c(4,2,3,1))
  for(j in 1:3){
    b<-boxplot(rmseSTD~method*nNA,ylim=c(0.1,1.4),
               ddff[which(ddff$y_j==paste("y",j)),],
               col="white",xaxt="n",
               xlab="Proportion of Missing Values",
               ylab="",outcol=cols[rep(1:nM,nn)],pch=15+rep(1:nM,nn),
               main=paste("y",j))
    axis(1,at = -3/2+(1:nn)*nM,labels = paste(props*100,"%",sep=""))
    abline(h=(0:15)/10,lty=2,col="gray")
    abline(v=+1/2+(0:5)*nM,col="gray")
    abline(h=c(0,1),col=cols[nM+1],lwd=1,lty=2)
    abline(h=sqrt(1-sqrt_1_minus_sig2^2),lty=1,col=cols[nM+2],lwd=1.5)
    text(5,0.1,labels = expression(sigma == 0.436),col=cols[nM+2])
    id_meth <- rep(1:nM,nn)
    coli_s <- cols[id_meth]
    for(ii in 1:length(coli_s)){
      rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
           ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           col=coli_s[ii])
    }
  }
  par(mar=c(0,0,0,0))
  ws <- 1/10
  plot(1:nM+1.1*ws,rep(0,nM),xlim=c(0,nM+1),ylim=c(-1,1)*0.5,
       col=coli_s,pch=15+(1:nM),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  legs <- c("Koh-Lanta (on ddsPLS)",
            "MI-NIPALS",
            "NIPALS",
            "MEAN-NIPALS")
  for(ii in 1:nM){
    rect(xleft = ii-ws/2,xright = ii+ws/2,
         ytop = ws,ybottom = -ws,
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],
         col=coli_s[ii])
    text(ii+ws/2,y = -0.4,labels = legs[ii],col=coli_s[ii])
  }
  # legend("center",legend = ,angle=angle,
  #        border=cols[1:nM],fill=cols[1:nM],density=density[1:nM],
  #        col = cols[1:nM],pch=15+(1:nM),
  #        ncol = 2,bty = "n",pt.cex=1)
}

plot_1_MNAR <- function(){
  ddff <- all_1_df
  ddff$method <- factor(ddff$method,levels = c("ddsPLS","MI-NIPALS","NIPALS","MEAN-NIPALS"))
  ddff$th <- factor(ddff$th)
  MM <- unique(ddff$method)
  nM <- length(MM)
  nn <- length(unique(ddff$th))
  NN <- 6
  layout(rbind(matrix(rep(1:3,NN),nrow = NN,byrow = T),c(1,1,1)*4))
  par(mar=c(4,2,3,1))
  for(j in 1:3){
    b<-boxplot(rmseSTD~method*th,ylim=c(0.1,1.4),
               ddff[which(ddff$y_j==paste("y",j)),],
               col="white",xaxt="n",
               xlab="Proportion of Missing Values",
               ylab="",outcol=cols[rep(1:nM,nn)],pch=15+rep(1:nM,nn),
               main=paste("y",j))
    axis(1,at = -3/2+(1:nn)*nM,labels = round(ths,2))
    abline(h=(0:15)/10,lty=2,col="gray")
    abline(v=+1/2+(0:5)*nM,col="gray")
    abline(h=c(0,1),col=cols[nM+1],lwd=1,lty=2)
    abline(h=sqrt(1-sqrt_1_minus_sig2^2),lty=1,col=cols[nM+2],lwd=1.5)
    text(5,0.1,labels = expression(sigma == 0.141),col=cols[nM+2])
    id_meth <- rep(1:nM,nn)
    coli_s <- cols[id_meth]
    for(ii in 1:length(coli_s)){
      rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
           ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           col=coli_s[ii])
    }
  }
  par(mar=c(0,0,0,0))
  ws <- 1/10
  plot(1:nM+1.1*ws,rep(0,nM),xlim=c(0,nM+1),ylim=c(-1,1)*0.5,
       col=coli_s,pch=15+(1:nM),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  legs <- c("Koh-Lanta (on ddsPLS)",
            "MI-NIPALS",
            "NIPALS",
            "MEAN-NIPALS")
  for(ii in 1:nM){
    rect(xleft = ii-ws/2,xright = ii+ws/2,
         ytop = ws,ybottom = -ws,
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],
         col=coli_s[ii])
    text(ii+ws/2,y = -0.4,labels = legs[ii],col=coli_s[ii])
  }
  # legend("center",legend = ,angle=angle,
  #        border=cols[1:nM],fill=cols[1:nM],density=density[1:nM],
  #        col = cols[1:nM],pch=15+(1:nM),
  #        ncol = 2,bty = "n",pt.cex=1)
}
