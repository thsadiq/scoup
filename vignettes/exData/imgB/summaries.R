# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Clear memory
rm(list=ls())

# ><>< # Specify data directory and filenames
pryDIR <- path.expand( file.path("~"))
dndsEst <- read.csv( file.path(pryDIR, "dnds.csv"), header=F, row.names=1)
omegaEst <- read.csv( file.path(pryDIR, "omega.csv") )
imgfile <- file.path(pryDIR, "FIG3.pdf")

# ><>< # Generate plot points
smmry <- function(abc){
    avg. <- mean(abc)
    s.e. <- sd(abc) / sqrt(length(abc))
    ans <- c(low=avg.-s.e., avg=avg., upp=avg.+s.e., se=s.e.)
    return(ans)
}
estNames <- colnames(omegaEst)
imposed <- apply(dndsEst, 1, smmry); dndsMat <- imposed[c(1,2,3),estNames]
inferred <- apply(omegaEst, 2, smmry); omegaMat <- inferred[1:3,estNames]

# ><>< # Estimate correlation coefficient of average estimates
svCor <- cor(dndsMat["avg",],omegaMat["avg",])
message( paste0("\nCorrelation Coefficient: ", sprintf("%.4f",svCor), "\n"))

# ><>< # Summarise standard errors
serror1 <- dndsMat["upp",] - dndsMat["avg",]
serror2 <- omegaMat["upp",] - omegaMat["avg",]
seTxt <- function(x.) paste(sprintf("%.4f",range(x.)), collapse=", ")
message( paste0("\nS.E. Range (dN/dS): ", seTxt(serror1), "\n"))
message( paste0("\nS.E. Range (omega): ", seTxt(serror2), "\n"))

# ><>< # Process Position
pdf(imgfile, height=15, width=15) 
  par(mar=c(0, 0, 0, 0))
  plot(0, 0, col="white", bty="n", xlab="", ylab="",
    xlim=c(.1,20.2), ylim=c(.32,20.8), xaxt="n", yaxt="n")

text(10.5, -0.15, expression(paste("OU Mean Reversion Rate, ",theta)), cex=2.5)
text(-.4, 10.5, expression(paste("dN/dS and ",omega)), cex=2.5, srt=90)
zWords <- expression(paste("OU Asymptotic Variance, ", Sigma^{2}))
text(10.5, 21.2, zWords, cex=2.5)

abX <- function(ygap) for(i in seq(1,20)) lines(rep(i,2), ygap, col="azure2")
abX(c(.7, 5.2)); abX(c(5.8, 10.2)); abX(c(10.8, 15.2)); abX(c(15.8, 20.15))

abY <- function(xgap) for(i in seq(1,20)) lines(xgap, rep(i,2), col="azure2")
abY(c(.8, 5.2)); abY(c(5.8, 10.2)); abY(c(10.8, 15.2)); abY(c(15.8, 20.2))

hLyne <- function(x1, x2, y0) lines(c(x1,x2), rep(y0,2), lwd=2)
phold <- sapply(seq(.8,16,by=5), function(a) hLyne(a, a+4.4, .8))
phold <- sapply(seq(.8,16,by=5), function(a) hLyne(a, a+4.4, 20.2))

vLyne <- function(y1, y2, x0) lines(rep(x0,2), c(y1,y2), lwd=2)
phold <- sapply(seq(.8,16,by=5), function(a) vLyne(a, a+4.4, .8))
phold <- sapply(seq(.8,16,by=5), function(a) vLyne(a, a+4.4, 20.19))

tickz <- function(syde, x1, x2, y1, y2, gp){
    if(syde==1) { sapply(seq(x1,x2,by=gp),
        function(a) lines(rep(a,2),c(y1,y2), lwd=2)) }else{
        sapply(seq(x1,x2,by=gp), function(a) lines(c(y1,y2),rep(a,2),lwd=2))}
}
tickz(2, 1, 20, .71, .78, 1)
sapply(seq(1,20,5) , function(a) tickz(1, a+.5, a+4, .65, .78, 1) )

htext <- function(strt, stop){ for(i in seq(strt,stop)){
    text(i+.5, .45, sprintf("%.2f",seq(.01,1,by=.33))[(i-strt)+1], cex=1.4) }
}
htext(1,5); htext(6,10); htext(11,15); htext(16,20)

vtext <- function(namz, bgn){ j <- 0; for(i in seq(bgn,bgn+4)) {
    j <- j+1; text(.8, i-.03, namz[j], cex=1.4, pos=2)}
}
n1s1 <- sprintf("%.2f",seq(0.0,1.0,.25)); vtext(n1s1, 1)
n1s2 <- sprintf("%.2f",seq(0.2,1.2,.25)); vtext(n1s2, 6)
n2s1 <- sprintf("%.2f",seq(0.0,1.0,.25)); vtext(n2s1, 11)
n2s2 <- sprintf("%.2f",seq(2.0,4.0,.50)); vtext(n2s2, 16)

# ><>< ===================== ><>< #
# ><>< Variables
# ><>< ===================== ><>< #
simSize <- 20
vNvSvec <- c(0, 5.00)
nsynVary <- c(0, 0.10)
eVary <- c(0.00, 0.34, 0.67, 1.00)
eThta <- c(0.01, 0.34, 0.67, 1.00)

# ><>< # Add Details to Plot
phold <- sapply(seq(1,length(eVary)), function(a){ xpoint <- seq(3,20,by=5)
    text(xpoint[a], 20.5, sprintf("%.2f",eVary[a]), cex=2, font=2)} )

initial <- 0
for(a in seq(1,length(nsynVary))){
    for(b in seq(1,length(vNvSvec))){
        if(nsynVary[a]==0){
            svar <- ifelse(vNvSvec[b]==0, 0, 1/vNvSvec[b])
        }else{
            svar <- ifelse(vNvSvec[b]==0, 0, round(nsynVary[a]/vNvSvec[b],2))}
        temp <- paste0("darwinTxt <- expression(paste(sigma[n]^{2}==",
            "\"", sprintf("%.2f",nsynVary[a]),"\", \"  &  \",",
            "sigma[s]^{2}==\"", sprintf("%.2f",svar), "\"))")
        eval(parse(text=temp))
        text(20.6,  3+(5*initial), darwinTxt, cex=1.7, srt=270)
        initial <- initial + 1
    }
}

midLyn <- function(a, b, h) lines(c(a,b),rep(h,2),lwd=3,lty=2,col="tomato")
useMidLyn <- function(yHeight){
    sapply(seq(0,15,by=5), function(a) midLyn(a+.8, a+5.2, yHeight)) }
useMidLyn(5); useMidLyn(9.25); useMidLyn(15)

# ><>< # Perform Plot Operation
sketch <- function(vect, spot, id="omg"){ extend <- 0.2
    if(id == "omg"){
        poynt <- 22; clr <- "gray60"; xID <- spot+extend
    }else{
        poynt <- 19; clr <- "gray5"; xID <- spot-extend
    }
    if((vect["upp"]-vect["low"]) > 1e-02){
        arrows(xID, vect["low"], xID, vect["upp"], .12, 90, 3, clr, lwd=3) }
    points(xID, vect["avg"], col=clr, pch=poynt, cex=2, lwd=2, bg="white")
}

shrink <- function(value, oldInt, newInt){
    newValue <- (value - oldInt[1]) * (diff(newInt) / diff(oldInt))
    newSpot <- newValue + newInt[1]
    return(newSpot)
}

# ><>< # Horizontal Axix Plot Points
hAxis <- function(pinD, pinE) (seq(1.5, 5, by=1)[pinE]) + ((pinD-1) * 5)

# ><>< # Plot Arrows
subplot <- function(vN, vS, oldWidth, newWidth){
    for(d in seq(1,length(eVary))){
        for(e in seq(1,length(eThta))){
            dLoc <- paste0("vN", vN, "vS", vS, "ouS", d, "ouT", e)
            dns <- dndsMat[,dLoc]; dns. <- shrink(dns, oldWidth, newWidth)
            omg <- omegaMat[,dLoc]; omg. <- shrink(omg, oldWidth, newWidth)
            sketch(omg., hAxis(d,e));   sketch(dns., hAxis(d,e), "dns")
        }
    }
}
subplot(1, 1, c(0,1.03), c(1,5.1))   # ><>< # vN1vS1: (vN=0, vS=0)
subplot(1, 2, c(.2,1.2), c(6,10.1))  # ><>< # vN1vS2: (vN=0, vS=5)
subplot(2, 1, c(0,1), c(11,15))      # ><>< # vN2vS1: (vN=0.1, vS=0)
subplot(2, 2, c(2,4), c(16,20))      # ><>< # vN2vS2: (vN=0.1, vS=0.02)

dev.off()

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #