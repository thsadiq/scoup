# ><>< ================================================================= ><>< #
# ><><                    Protein Sequence Simulation                    ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><   R code used to generate the outputs that are presented in the   ><>< #
# ><><                          JOSS manuscript                          ><>< #
# ><>< ================================================================= ><>< #

# ><>< # Path to Data Folder
databox <- "data/"

# ><>< # Bespoke Summary Function
smmry <- function(vect){
    avg. <- mean(vect)
    s.e. <- sd(vect) / sqrt(length(vect))
    ans <- c(low=avg.-s.e., avg=avg., upp=avg.+s.e., se=s.e.)
    return(ans)
}

# ><>< # Tidy Plot Function
splitarrow <- function(estmatrix, pnt, clr, adjst){
    size <- ncol(estmatrix)
    data <- apply(estmatrix, 2, smmry)
    putline <- vapply(seq(1,size), function(x){
        arrows(x-adjst, data["low",x], x-adjst,
            data["upp",x], 0.1, 90, 3, clr, lwd=2)
        return(0)}, FUN.VALUE=0)
    points(seq(1,8)-adjst, data["avg",],
        pch=pnt, bg="white", col=clr, cex=2, lwd=2)
}

# ><>< # Retrieve Data
s00dnds <- read.csv(paste0(databox,"s000dnds.csv"))
s01dnds <- read.csv(paste0(databox,"s010dnds.csv"))
s00paml <- read.csv(paste0(databox,"s000omega.csv"))
s01paml <- read.csv(paste0(databox,"s010omega.csv"))

# ><>< # Default `mar`
tags <- c("", "0.001", "0.005", "0.009",
    "0.030", "0.070", "0.100", "0.500", "0.900", "")

pdf("FIG2.pdf", width=15, height=7) 
    par(mar=c(4.2, 4.2, 0.15, 0.15), fig=c(0, 0.48, 0, 1))
    plot(0, 0, col="white", bty="n", xlab="", ylab="",
        xlim=c(0.8,8.2), ylim=c(0,0.9), xaxt="n", yaxt="n")
    splitarrow(s00dnds, 19, 'gray30',  0.2)
    splitarrow(s00paml, 22, 'gray70', -0.2)
    axis(1,c(4,5), tags[c(5,6)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    axis(1,seq(0,3),tags[seq(1,4)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    axis(1,seq(6,9),tags[seq(7,10)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    axis(2, seq(-.2,1,0.2), las=1, tck=-0.01, lwd=2, cex.axis=1.2, hadj=0.7)
    text(c(1.5,1.5), c(.78,.71), c(expression(omega),"dN/dS"), cex=2, pos=4)
    text(1.00, 0.85, expression(sigma[s]^2=="0.00"), cex=2, pos=4)
    mtext(expression(paste(omega,", dN/dS")), 2, 2.2, cex=2)
    points(c(1.3,1.3), c(.78,.71), pch=c(22,19), cex=2)
    mtext(expression(sigma[n]^2), 1, 3.3, cex=2)
    box(lwd=2, bty="7")

    par(fig=c(0.52, 1, 0, 1), new=TRUE)
    plot(0, 0, col="white", bty="n", xlab="", ylab="",
        xlim=c(0.8,8.2), ylim=c(0.9,10), xaxt="n", yaxt="n")
    splitarrow(s01dnds, 19, 'gray30',  0.2)
    splitarrow(s01paml, 22, 'gray70', -0.2)
    axis(1,c(4,5), tags[c(5,6)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    axis(1,seq(0,3),tags[seq(1,4)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    axis(1,seq(6,9),tags[seq(7,10)], tck=-0.01, lwd=2, padj=-.6, cex.axis=1.2)
    text(c(1.5,1.5), c(9.0,8.3), c(expression(omega),"dN/dS"), cex=2, pos=4)
    axis(2,seq(0,8,2),las=1, tck=-0.01, lwd=2, cex.axis=1.4, hadj=0.2)
    axis(2, 10, las=1, tck=-0.01, lwd=2, cex.axis=1.4, hadj=0.6)
    text(1.00, 9.7, expression(sigma[s]^2==0.01), cex=2, pos=4)
    mtext(expression(paste(omega,", dN/dS")), 2, 1.5, cex=2)
    points(c(1.3,1.3), c(9.0,8.3), pch=c(22,19), cex=2)
    axis(2, seq(8,12,2), rep("",3), tck=0, lwd=2)
    mtext(expression(sigma[n]^2), 1, 3.3, cex=2)
    box(lwd=2, bty="7")
dev.off()


# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #
