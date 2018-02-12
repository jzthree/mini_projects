#This script provide function to do hypothesis testing on weighted correlation coefficients.
#Note that this is specifically for the case where weights are fixed and known, and are NOT 
#random variables.

require(wCorr)


#t-test similar to regular t-test except replacing the sample size N with
#the effective sample size (sum w_i)^2 / (sum w_i^2)
weighted.t.cor<-function(x,y,weight,method='Pearson'){
  validind = !is.na(x) & !is.na(y) & !is.na(weight)
  score = weightedCorr(x[validind],y[validind], weights = weight[validind],method =method)

  effective_N = (sum(weight[validind],na.rm=T)^2) /sum(weight[validind]^2 ,na.rm=T)
  t = score*sqrt((effective_N-2)/(1-score^2))
 
  return(c(cor=score,p.value= 2*pt(abs(t), effective_N-2, lower.tail = FALSE),dof = effective_N-2))
}

#z-test based on Fisher transformation, this is approximate but should be
#very accurate when the effective sample size is high.
#An advantage of z-test is it can be easily modified to test null hypothesis 
#of H0: r=r_0 where r_0 != 0
#z-test for Spearman's rho is based on https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient

weighted.z.scor<-function(x,y,weight,method='Spearman'){
  validind = !is.na(x) & !is.na(y) & !is.na(weight)
  score = weightedCorr(x[validind],y[validind], weights = weight[validind],method =method)
  if (method == 'Spearman')
    c = 1.06
  else
    c = 1
  effective_N = (sum(weight[validind],na.rm=T)^2) /sum(weight[validind]^2 ,na.rm=T)
  z = fisherz(score) 
  return(c(cor=score,p.value=2*pnorm(-abs((z)*sqrt((effective_N-3)/c))),dof = effective_N-3))
}
