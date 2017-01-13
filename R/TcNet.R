  # Legendre Polynominals
  
  LgdP <- expression( tt,
                      ( 3* tt^2 - 1 )/2 , 
                      ( 5 * tt^3 - 3* tt )/2, 
                      ( 35 * tt^4 - 30 * tt^2 + 3)/8,
                      ( 63 * tt^5 - 70 * tt^3 + 15 * tt )/8,
                      ( 231 * tt^6 - 315 * tt^4 + 105 * tt^2 - 5)/16,
                      ( 429 * tt^7 - 693 * tt^5 + 315 * tt^3 - 35 * tt)/16,
                      ( 6435 * tt^8 - 12012 * tt^6 + 6930 * tt^4 - 1260 * tt^2 + 35)/128,
                      ( 12155 * tt^9 - 25740 * tt^7 + 18018 * tt^5 - 4620 * tt^3 + 315 * tt)/128,
                      ( 46189 * tt^10 - 109395 * tt^8 + 90090 * tt^6 - 30030 * tt^4 + 3465 * tt^2 - 63)/256 )
  
  GetMR <- function(rho,times)
  {
    MR <- matrix(1,length(times),length(times))
    for ( i in 1:length(times)){
      for(j in 1:length(times)){
        MR[i,j]= rho^(abs(times[j] - times[i]))
      }
    }
    return (MR)
  }
  GetMX <- function(times,r)
  {
    tnum = length(times) 
    X <- matrix(1,tnum,r+1)
    
    for(t in 1:tnum ){
      tt <- -1 + 2*(times[t] - times[1])/(times[tnum] - times[1])
      for(i in 1:r){
        X[t,i+1] <- eval(LgdP[i])
      }
    }
    return (X)
  }
  GetInitPij <- function(N,J)
  {
    P <- matrix(1/J,N,J)
    for (i in 1:N){
      P[i,] <- rnorm(J, mean=1/J, sd= 0.2 * 1/J )
      P[i,] <- P[i,]/sum(P[i,])
    }
    
    return (P)
  }
  GetMeanMatrix <- function(J,times,P,X,Asdata,InvMSigema)
  {
    m <- matrix(NA,J,length(times))
    N <- length(Asdata[,1])
    r <- length(X[1,])
    
    xInvSigema <- t(X) %*% InvMSigema
    xInvSigemax <- xInvSigema%*% X
    
    mU <- matrix(NA, J, ncol(X))
    for( j in 1:J){
      ud <- matrix(0, r, r)
      for( i in 1: N){
        ud <- ud + P[i,j]*xInvSigemax
      }
      ubd <- matrix(0, r, 1)
      for( i in 1: N){
        ubd <- ubd + P[i,j]*( xInvSigema %*% (Asdata[i,]) )
      }
      uj <- ginv(ud) %*% ubd
      m[j,] <- X %*% uj
      mU[j,] <- uj
    }
    return(list(M = m, U = mU))
  }
  GetNewSsquare <- function(Asdata,m,MR,times,P,J)
  {
    N <- length(Asdata[,1])
    
    InvMR <- ginv(MR)
    newSsquare <- 0
    for(i in 1:N){
      SumJ <- 0
      for(j in 1:J){
        yi_mj <- Asdata[i,]-m[j,]
        SumJ <- SumJ + P[i,j] * ((yi_mj) %*% InvMR %*% (yi_mj) )
      }
      newSsquare <- newSsquare + SumJ
    }
    newSsquare <- as.numeric(newSsquare/(length(times)*N))
    
    return(newSsquare)
  }
  GetNewRho.b <- function(rho,rhoDir)
  {
    
    newrho <- as.numeric(rho + 0.004*rhoDir)
    if (newrho > 0.99) newrho <- 0.99
    if (newrho < 0) newrho <- 0
    
    return (newrho)
  }
  GetNewRho<- function(Asdata,m,MR,times,P,J,rho,Ssquare)
  {
    N <- length(Asdata[,1])
    newrho <- 0
    for(i in 1:N){
      SumJ <- 0
      for(j in 1:J){
        yi_mj <- Asdata[i,]-m[j,]
        Item1 <- (1/(1 - rho*rho))*((yi_mj) %*% MR %*% (yi_mj) )
        Item2 <- 0
        for(k in 2:(length(times)-1) )
          Item2 <- Item2 + (yi_mj[k]^2)
        Item2 <- Item2 * rho
        Item3 <- 0
        for(k in 1:(length(times)-1) )
          Item2 <- Item3 + yi_mj[k] * yi_mj[k+1]
        SumJ <- SumJ + P[i,j] * (Item1 + Item2 - Item3)
      }
      newrho <- newrho + SumJ
    }
    newrho <- as.numeric(newrho/( (length(times)-1)* N * Ssquare))
    
    if(abs(newrho) >= 1) return( sign(newrho)*.5)
    else return(newrho)
  }
  GetLikelihood <- function(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times) 
  {
    N <- length(Asdata[,1])
    
    LogDetMSigema <- log(DetMSigema)/2
    LogM2Pi <- length(times)*log(2*pi)/2
    
    oneterm <- function(i, j) {
      f <- function(i,j)
        P[i,j]*(log(omiga[j]) - LogM2Pi - LogDetMSigema
                - ( ((Asdata[i,]-m[j,])) %*% InvMSigema %*% (Asdata[i,]-m[j,])) /2)
      mapply(f, i, j)
    }
    tmp <- outer(1:N, 1:J, oneterm)
    tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
    return(sum(tmp))
  }
  StepE <- function(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times)
  {
    TwoPiExp <- (2*pi) ^ ( length(times)/2 )
    TwoPiExp <- TwoPiExp * DetMSigema
    
    N <- length(Asdata[,1])
    
    tmp <- rep(0,N)
    for( i in 1:N){
      Fi <- rep(0,J)
      for( j in 1:J){
        yi_mj <- Asdata[i,]-m[j,]
        Fi[j] = exp( ( (yi_mj) %*% InvMSigema %*% (yi_mj) ) / -2) / TwoPiExp 
      }
      OmigaF <- omiga %*% Fi
      P[i,] <- (omiga * Fi) / OmigaF
      tmp[i] <- log(OmigaF)
    }
    tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
    Likelihood <- sum(tmp)
    
    return(list(P = P, Likelihood = Likelihood))
  }
  
  StepM.b <- function(Asdata,m,MR,times,Ssquare,P,rho,rhoDir,J,rpt)
  {
    newSsquare <- GetNewSsquare(Asdata,m,MR,times,P,J)
    if (rpt > 0)
      newrho <- GetNewRho.b(rho,rhoDir)
    else
      newrho <- rho
    
    return( c(newSsquare, newrho))
  }
  
  
RunEM.C <- function(Asdata,times,rho,Ssquare,X,P,MR,MSigema,InvMSigema,DetMSigema,omiga,m,J,r)
{
  IncLimit <- 3
  REPEAT_LIMIT <- 200
  LIKELIHOOD_DIFF <- 0.01
  
    rpt <- 1
    Likelihood <- -Inf
    
    rhoDir <- 1
    rhoIncCount <- 0
    
    while(TRUE){
      OldLikelihood <- Likelihood
      EResult <- StepE(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times)
      P <- EResult$P
      Likelihood <- EResult$Likelihood
      
      if( (abs(OldLikelihood - Likelihood) < LIKELIHOOD_DIFF) ){
        #cat("quit due to likelihood\n")
        #cat("LIKELIHOOD_DIFF:",LIKELIHOOD_DIFF,"\n")
        break
      }
      
      if( rpt >= REPEAT_LIMIT ){
        #cat("quit due to rpt\n")
        break
      }
      
      if ( Likelihood >= OldLikelihood){
        #rhoIncCount <- 0
      }else{
        #cat("decrease, limit:", IncLimit, "\n")
        rhoIncCount <- rhoIncCount + 1
        if (rhoIncCount >= IncLimit){
          rhoIncCount <- 0
          rhoDir <- rhoDir * -1
        }
      }
      
      newpars <- StepM.b(Asdata,m,MR,times,Ssquare,P,rho,rhoDir,J,rpt)
      #     cat("newpars:\n")
      #     print(newpars)
      Ssquare <- newpars[1]
      rho <- newpars[2]
      
      MR <- GetMR(rho,times)
      MSigema <- Ssquare * MR
      # print(MSigema)
      InvMSigema <- ginv(MSigema)
      DetMSigema <- (det(MSigema))^0.5
      
      N <- length(Asdata[,1])
      omiga <- colSums(P)/ N 
      
      rMeans <- GetMeanMatrix(J,times,P,X,Asdata,InvMSigema)
      m <- rMeans$M
      mU <- rMeans$U
      rpt <- rpt + 1
    }
    
    return(list(rho=rho,Ssquare=Ssquare,Likelihood=Likelihood,m=m,mU=mU,P=P))
}
  
  
  InitAndRunEM <- function(mExpression, times, J=5, r=4)  # J: Genes; r :order of legendre polynominal
  {
    rho <- 0.8
    Ssquare <- 20
    
    X <- GetMX(times,r)
    N <- length(mExpression[,1])
    P <- GetInitPij(N,J)
    
    MR <- GetMR(rho,times)
    MSigema <- Ssquare * MR
    #print(MSigema)
    InvMSigema <- ginv(MSigema)
    DetMSigema <- (det(MSigema))^0.5
    
    omiga <- colSums(P)/ N
    m <- GetMeanMatrix(J,times,P,X,mExpression,InvMSigema)$M
    
    EMResults <- RunEM.C(mExpression,times,rho,Ssquare,X,P,MR,MSigema,InvMSigema,DetMSigema,omiga,m,J,r) 
  }

GeneCluster <-
function(mExpression, times, NumberOfCluster,orderLOP)
{ 
  r=orderLOP
  J=NumberOfCluster
  
  if(!requireNamespace("MASS",quietly = TRUE)){
    stop("MASS needed for this function to work. Please install it.",
         call.=FALSE)
  }
  

  if (r>10) {warning("range of r is 1 to 10")}
  if (r<0) {warning("r must be greater than 0")}
  if (length(times)!= ncol(mExpression)) {warning(" length of time vector is different to column length ")}
  if (J>nrow(mExpression)) {warning("number of cluster larger than number of genes")}
  if (!is.numeric(mExpression)) {warning("data is not numeric ")}
  
     
  lgd=r
  
  
  a1 <- times
  LineColors <- colors()
  
  opar <- par(no.readonly = TRUE)
  par(col.axis="blue", mar=c(4, 4, 2.5, 0.25),pty = "m",lwd=1)
  
  Asdata <- mExpression
  
  tmp<-InitAndRunEM(mExpression, times, J,r)
  
  Pij <-tmp$P
  m <- tmp$m
  vu <-tmp$mU
  
  indx <- max.col(Pij)
  grp <- unique(indx)
  grp <- order(grp) # added by Yaqun Wang, 1/10/2017
  
  TransColors <- rainbow(length(grp),alpha = 0.12)
  SolidColors <- rainbow(length(grp),alpha = 0.8)
  
  vY <- c(min(mExpression), max(mExpression), rep(0,length(a1) - 2))
  plot(a1, vY , type = "n",xlab="Time", ylab="Expression")
  TitleStr <- sprintf("Summary of Gene Expression( Clusters: %d )",length(grp))
  title(TitleStr, font.main=3)
  
  for( i in 1:length(grp) ){
    y18 <- as.numeric(m[grp[i],])
    par(lwd = 2)
    lines(a1, y18, col=SolidColors[i])
    par(lwd = 1)
  }
  
  for( i in 1:length(grp) ){
    indxTemp <- which(indx == grp[i])
    
    plot(a1, vY , type = "n",xlab="Time", ylab="Expression")
    TitleStr <- sprintf("Gene Expression of Cluster%d, Genes:%d",grp[i],length(indxTemp)) # Modified by Yaqun wang, i => grp[i]
    title(TitleStr, font.main=1,cex.main = 1)

    
    for(j in 1: length(indxTemp) ){
      lines(a1,Asdata[indxTemp[j], ],col=TransColors[i])
    }
    par(lwd = 3)
    y18 <- as.numeric(m[grp[i],])
    lines(a1, y18, col=SolidColors[i])
    par(lwd = 1)
  }
  par(opar)
  
  return(list(MeanExpression=m,LOPCoefficient =vu,Classifications=indx))
}
GeneClusterBIC <-
function( mExpression, times, G = c(1:15), orderLOP)
{ 
  if(!requireNamespace("MASS",quietly = TRUE)){
    stop("MASS needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  r=orderLOP
  lgd=r
  if (r>10) {warning("range of r is 1 to 10")}
  if (r<0) {warning("r must be greater than 0")}
  if (length(times)!=ncol(mExpression)) {warning(" length of time vector is different to column length ")}
  if (!is.numeric(mExpression)) {warning("data is not numeric ")}
  
  EMResults <- matrix(0, length(G) , 4)
  i <- 1
  for(j in G){
  	cat("G: ", j, "\n")
    EMResults[i , 1] <- j
    tmp<-InitAndRunEM(mExpression, times, J=j,r=lgd)
    EMResults[i , 2:4] <- c(tmp$rho,tmp$Ssquare,tmp$Likelihood)
    i <- i + 1
  }
  m <- ncol(mExpression)
  N <- nrow(mExpression)
  
  BIC <- -2 * EMResults[ , 4] + log(N) * (EMResults[ , 1] * (m + 1) + 1)
  BIC <- rbind(EMResults[ , 1], BIC)
  
  markSmallest=FALSE
  startN=1
  par(mar=(c(4,4,2,2)+.1),cex=1.3)
  Nk <- ncol(BIC)
  plot(BIC[1, startN : Nk], BIC[2, startN : Nk] / 1e+3 , type="b",pch = 5 , col= "blue",
       lty = "solid",main="",xlab = "# of cluster", ylab = "BIC", cex=0.9)
    ind <- which( BIC[2, ] == min(BIC[2, ]) )
    points(BIC[1, ind], BIC[2, ind] / 1e+3, 
           col = "deeppink", pch = 8 )
  mtext(expression("x"*10^3),line = 0.3,adj = 0,cex=1.2)
  
  optimal=which(BIC[2,]==min(BIC[2,])) 
  
  return(list(BIC=BIC,optimal=optimal))

}
GetPt <- function(t, times, r)
{
  tnum <- length(times)
  Pt<- rep(1,r+1)
  tt <- -1 + 2 * (t - times[1])/(times[tnum] - times[1])
  for(i in 1:r){
    Pt[i+1] <- eval(LgdP[i])
  }
  return(Pt)
}

GetYHat <- function(xx, LOPCoefficient, times, r)
{  
  yhat <- matrix(0,nrow(LOPCoefficient),length(xx))
  for(i in 1:length(xx)){
    Pt <- GetPt(xx[i],times, r)
    yhat[,i] <- Pt %*% t(LOPCoefficient)
  }
  return(yhat)
}

GeneClusterInterp <-
function(LOPCoefficient, OriginalTime, outLen = 20)
{ 
  times=OriginalTime

  r <- length(LOPCoefficient[1,]) - 1
  newtime  <- seq(min(times),max(times), len=outLen)
  yHat <- GetYHat(newtime, LOPCoefficient, times, r)
  return(rbind(newtime,yHat))
}
ModifyMatrix <-function(LOPCoefficient,times,lowCut,upCut)
{
      
      cut=((max(times)-min(times))+1)*10
      output2=GeneClusterInterp(LOPCoefficient,times,cut)
      newoutput<-matrix(NA,nrow(LOPCoefficient),cut)
      for(i in 2:nrow(output2)){
        if((output2[i,1]> upCut )|(output2[i,1]< lowCut)){newoutput[i-1,]=output2[i,]}
        else{
          for(j in 1:(ncol(output2)-1)){
            if (((output2[i,j] < upCut )&(output2[i,j] > lowCut)) & ((output2[i,j+1] > upCut)|(output2[i,j+1] < lowCut))){
              newoutput[i-1,1:(cut-j)]=output2[i,(j+1):cut]
              break}
          }
        }
      }
      
      rownames(newoutput)=seq(1,nrow(LOPCoefficient))
      
      na_count<- apply(newoutput, 1, function(z) sum(is.na(z)))
      
      Modify=newoutput[!(na_count==cut),]
      
      Modify=na.omit(t(Modify))
      
      Modify=t(Modify)
      
      na_count=data.frame(na_count)
      
      return(list(na_count=na_count,Modify=Modify))
}
  
GeneClusterNet <-
function(mExpression,times, orderLOP,alpha1=0.5,alpha2=0.05, realign=F, cutoff=c(lowCut=-0.35,upCut=0.2),NumberOfCluster=0,sLabels=NULL)
{ 
  if(!requireNamespace("G1DBN",quietly = TRUE)){
    stop("G1DBN needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  if(!requireNamespace("igraph",quietly = TRUE)){
    stop("igraph needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  
  
  data=mExpression
  
  J=NumberOfCluster
  
  if(NumberOfCluster==0){J=GeneClusterBIC(data, times, G = c(1:15), orderLOP)$optimal}
  
  vu=GeneCluster(data, times, J, orderLOP)$LOPCoefficient
  
  if(realign){Xn2=ModifyMatrix(vu,times,cutoff[1],cutoff[2])$Modify}
  
  else {Xn2=GeneClusterInterp(vu, times, outLen = (max(times)-min(times)+1))[-1,]}
  
  kk=(Xn2[,1]>0)
  
  color=rep("yellow",length(kk))
  for(i in 1:length(kk)){
    if(kk[i]=="FALSE"){color[i]="red"} 
  }
  
  Xn=t(Xn2)
  S1 <- DBNScoreStep1(Xn, method='ls')
  S2 <- DBNScoreStep2(S1$S1ls, data=Xn, method='ls', alpha1=alpha1)
  G2 <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
  Step2InferredNet<- BuildNetwork(G2,1:nrow(Xn2))
  
  
#   if(NumberOfCluster==0) {sLabel=seq(1:nrow(Xn2))
#   }else {if(nrow(Xn2)==J){sLabel=sLabels}
#    else {sLabel=seq(1:nrow(Xn2))
#    }
#   }
  
  if(NumberOfCluster==0) {sLabel=seq(1:nrow(Xn2))
  }else {if( (nrow(Xn2)==J) & (length(sLabels) !=0 )){sLabel=sLabels} # Nodified by Yaqun 1/9/2017
  	else {sLabel=seq(1:nrow(Xn2))
  	}
  }
  
  
  if(realign){sLabel=row.names(Xn2)
  
  rownames(Step2InferredNet$AdjMatrix)=row.names(Xn2)
  colnames(Step2InferredNet$AdjMatrix)=row.names(Xn2)
  
  rownames(Step2InferredNet$Score)=row.names(Xn2)
  colnames(Step2InferredNet$Score)=row.names(Xn2)

  }
  
  g <- graph.adjacency(t(Step2InferredNet$AdjMatrix), mode ='directed')
  g <- simplify(g, remove.multiple = T, remove.loops = T)
  #V(g)$color <- color
  V(g)$color <- 'deepskyblue3'
  E(g)$color <- 'deeppink'
  E(g)$arrow.width <- 0.5
  plot(g, layout=layout.circle, vertex.label=sLabel, vertex.label.dist=0,vertex.size=15)
  title(main="DBN Inferred network") 
  
  
  return(list(Score=Step2InferredNet$Score ,AdjMatrix=Step2InferredNet$AdjMatrix))
}

