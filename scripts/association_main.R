
### Jinliang Yang
### updated, July 28th, 2014

### modified by Niranjani, Jan, 2020 : To utilize for Rice traits

library(ggplot2)
library(dplyr)
library(qqman)
############################
chistat <- function(snpdf){
  snpdf <- as.matrix(snpdf)
  stat <- apply(snpdf, 1, function(y){    
    x=1:3
    tab=matrix(y,byrow=T,ncol=2)
    o=glm(cbind(tab[,1],tab[,2])~x,
          family=binomial(link=logit))
    chisqStat=anova(o)[2,2]
    chisqStat
  })  
}

#### loop through all SNPs
get_chistat <- function(snps=snps, filter=50){
    ### snps: a data.frame format input file
    ### filter: minimum depth for a given site, default=50
    
    snps$low_total <- snps$low_ref + snps$low_alt
    snps$random_total <- snps$random_ref + snps$random_alt
    snps$high_total <- snps$high_ref + snps$high_alt
    snp50 <- subset(snps, high_total > filter & low_total > filter & random_total > filter)
    snp50 <- snps
    
    if (dim(snp50)[1] ==0) {
      message(sprintf(">>> There are ZERO input variants after depth filtering [ > %s ].", 
                      filter))
      stop("Terminating the programm")
    }

    message(sprintf("###>>> input [ %s ] variants, remaining [ %s ] after depth filtering [ > %s ]", 
                    nrow(snps), nrow(snp50), filter))
   
    out <- snp50[, 1:3]
    snpcol <- snp50[, c("high_ref", "high_alt", "low_ref", "low_alt", "random_ref", "random_alt")]
    
    ### add 1 for each column
    snpcol$high_ref <- snpcol$high_ref + 1
    snpcol$high_alt <- snpcol$high_alt + 1
    snpcol$low_ref <- snpcol$low_ref + 1
    snpcol$low_alt <- snpcol$low_alt + 1
    snpcol$random_ref <- snpcol$random_ref + 1
    snpcol$random_alt <- snpcol$random_alt + 1
    
    out$stat <- chistat(snpdf = snpcol)
    message(sprintf("###>>> DONE!"))
    return(out)
}


### function copied from GenABEL and modified by Jingliang Yang
estlambda <- function (data, plot = FALSE, proportion = 1, method = "regression", 
          filter = TRUE, df = 1, ...) 
{
  data <- data[which(!is.na(data))]
  if (proportion > 1 || proportion <= 0) 
    stop("proportion argument should be greater then zero and less than or equal to one")
  ntp <- round(proportion * length(data))
  if (ntp < 1) 
    stop("no valid measurements")
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10) 
    warning(paste("number of points is too small:", ntp))
  if (min(data) < 0) 
    stop("data argument has values <0")
  if (max(data) <= 1) {
    data <- qchisq(data, 1, lower.tail = FALSE)
  }
  if (filter) {
    data[which(abs(data) < 1e-08)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  out <- list()
  if (method == "regression") {
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
  }
  else if (method == "median") {
    out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
                                                      df)
    out$se <- NA
  }
  else if (method == "KS") {
    limits <- c(0.5, 100)
    out$estimate <- estLambdaKS(data, limits = limits, df = df)
    if (abs(out$estimate - limits[1]) < 1e-04 || abs(out$estimate - 
                                                       limits[2]) < 1e-04) 
      warning("using method='KS' lambda too close to limits, use other method")
    out$se <- NA
  }
  else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  if (plot) {
    lim <- c(0, max(data, ppoi, na.rm = TRUE))
    oldmargins <- par()$mar
    par(mar = oldmargins + 0.2)
    plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
         ylab = expression("Observed " ~ chi^2), ...)
    abline(a = 0, b = 1)
    #abline(a = 0, b = out$estimate, col = "red")
    par(mar = oldmargins)
  }
  out
}

### Jinliang Yang

get_qval <- function(pval=pval1, pcol="P1df", method="fdr"){
  qval <- p.adjust(pval[,pcol], method=method)
  return(qval)
}

xpgwas_qval <- function(stat=out, lambda=1){
    
    stat$pval <- pchisq(stat$stat/lambda, 1, lower.tail=FALSE)
    stat$qval <- get_qval(pval=stat, pcol="pval", method="fdr")
    #[1] 134 ?
    message(sprintf("###>>> [ %s ] significant sites using FDR < 0.05", nrow(subset(stat, qval < 0.05))))
    return(stat)
}


#### simple XP-GWAS plot

xpplot <- function(qval){
    qval$log10p <- -log10(qval$pval)
    plot(x=qval$pos, y=qval$log10p,  type="p", xaxt="n", pch=20, col="slateblue", bty="n",
         xlim=c(0, max(qval$pos)), ylim=c(0.5, 10 ), ylab="-log10(P-value)", xlab="Physical Position (bp)")
    cutoff <- min(subset(qval, qval < 0.05)$log10p)
    abline(h=cutoff, col="red", lty=2, lwd=2)
}  

xpplot_mine <- function(qval){
    qval$log10p <- -log10(qval$pval)
    b <- qval[,c(2,3,1,5)]
    names(b) <- c("CHR","BP","SNP","P")
    png(filename="ManhattanPlot.png",width = 480,height = 400)
    manhattan(b)
    dev.off()
    #ggsave(filename = "manhattan1.png",p1, width = 12, height = 8)
    
}

