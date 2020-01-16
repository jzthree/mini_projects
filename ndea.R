#Network differential enrichment analysis
#Zhou J, Park CY, Theesfeld CL, Wong AK, Yuan Y, Scheckel C, Fak JJ, Funk J, Yao K, Tajima Y, Packer A, Darnell RB, Troyanskaya OG. (2019). Whole-genome deep-learning analysis identifies contribution of noncoding mutations to autism risk. Nature Genetics.
#Note: the original publication contains an error in the NDEA method description: the weight W_ij(m) in the publication is the network edge score, not the network edge score divided by number of mutations per gene. We will correct this with the publisher.

ndea <- function(x, y, gene_indices_x, gene_indices_y, network, threshold=0, 
                 alternative='greater', mc.cores=10){
  network = network * (network>threshold)
  
  results <- do.call(rbind,
    mclapply(1:dim(network)[2], 
    function(i)
      weighted.t(x, y, 
                 as.numeric(network[gene_indices_x,i]),
                 as.numeric(network[gene_indices_y,i]), alternative=alternative), 
    mc.cores=mc.cores))
  return(results)
}



weighted.t <- function(x,y,weightx,weighty,alternative='greater'){
  wtd.var<-function(x,weights,na.rm=T){
    if (na.rm){
      s <- !is.na(x + weights)
      x <- x[s]
      weights <- weights[s]
    }
    as.numeric(stats::cov.wt(cbind(x), weights, method = "unbiased")$cov)
  }
  wtd.mean<-function(x,weights,na.rm=T){
    if (na.rm){
      s <- !is.na(x + weights)
      x <- x[s]
      weights <- weights[s]
    }
    sum(weights * x)/sum(weights)
  }
  tryCatch({
    validindx = is.finite(x) & is.finite(weightx)
    validindy = is.finite(y) & is.finite(weighty)
    mx <- wtd.mean(x[validindx], weightx[validindx], na.rm = TRUE)
    vx <- wtd.var(x[validindx], weightx[validindx], na.rm = TRUE)
    my <- wtd.mean(y[validindy], weighty[validindy], na.rm = TRUE)
    vy <- wtd.var(y[validindy], weighty[validindy], na.rm = TRUE)
    effective_N1 = (sum(weightx[validindx],na.rm=T)^2) /sum(weightx[validindx]^2 ,na.rm=T)
    effective_N2 = (sum(weighty[validindy],na.rm=T)^2) /sum(weighty[validindy]^2 ,na.rm=T)
    sxy <- sqrt((vx/effective_N1) + (vy/effective_N2))
    df = (((vx/effective_N1) + (vy/effective_N2))^2)/((((vx/effective_N1)^2)/(effective_N1 - 1)) +
                                                        ((vy/effective_N2)^2/(effective_N2 - 1)))
    return(c(t=(mx - my)/sxy,p.value=(1 - pt(abs((mx - my)/sxy), df)) * 2))},                            
    error=function(e)return(c(t=NA,p.value=NA)))
}
