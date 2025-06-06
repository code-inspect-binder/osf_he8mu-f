

###########################################
# Create a function to correct for rounding error in p-values
###########################################
###########################################
# Bound p-values to avoid errors that may come about if R rounds p to 0
# Credit to Lee Jussim for identifying this bug
# Credit to Uri Simonsohn for function
pbound=function(p) pmin(pmax(p,2.2e-16),1)

##Count non-missing
CountNM <-function(x) {
  sum( !is.na( x )) }

#compute estimated power using Schimmack's (2014) method
EstPower <-function (pvalues) {
    pvalues<-pmin(pmax(pvalues,2.2e-16),1-2.2e-16)
    EstPower <-pnorm(qnorm(1-pvalues/2) - qnorm(.975)) 
    return(EstPower)}


#PCurve related functions
#Adapted from http://www.p-curve.com/app4/ on 10/10/16. needed for pcurve based on our already coverted tests to pvalues
stouffer.P=function(pp) pnorm(sum(qnorm(pp),na.rm=TRUE)/sqrt(sum(!is.na(pp))))

#Taken directly from from http://www.p-curve.com/app4/ on 10/10/16. needed for pcurve
stouffer=function(pp) sum(qnorm(pp),na.rm=TRUE)/sqrt(sum(!is.na(pp)))

##############################
###Effect size calcuations
#############################

#Chi
wchiest <-function (chisq,n) {
  omegachi <- (chisq / n)^.5 }

#Estimate d from t and df
dest <-function (t,df) {
  r2 <- (t^2) / (t^2 + df)
 d.est <- (2*(r2^.5)) / (1-r2) }

#Estimate n2 from F and Df
cohenFfromN2 <- function(Fvalue,dfn,dfd) {
  n2 <- (dfn*Fvalue) / (dfn*Fvalue + dfd) 
  cohen.f2 <-n2/(1-n2) }

#Estimate r2 from T and df
R2fromT <-function (t,df) {
  r2 <- (t^2) / (t^2 + df)}

#Estimate r2 from F and df
R2fromF <- function(Fvalue,dfn,dfd) {
  n2 <- (dfn*Fvalue) / (dfn*Fvalue + dfd) }


r2toD <-function (r2) {
  d <-(2*(r2^.5)) / ((1-r2)^.5)
  return(d)
}


###########################################
##################APA Pvalue functions
###########################################
Prounder<-function(pvalue){
  if (is.na(pvalue)==TRUE){
    Prounder <- NA
  } else if (pvalue>.05) {
    Prounder <- NA
  } else if (pvalue<=.05 & pvalue>.01) {
    Prounder <- .05
  } else if (pvalue<=.01 & pvalue>.001) {
    Prounder <- .01
  } else if (pvalue<=.001) {
    Prounder <- .001
  } else 
    Prounder <- NA
  return(Prounder)
}

PvalueAPA6<-function(pvalue){
  if (is.na(pvalue)==TRUE){
    PvalueAPA6 <- NA
  } else if (pvalue<=.001) {
    PvalueAPA6 <- .001
  } else 
    PvalueAPA6 <- round(pvalue,3)
  return(PvalueAPA6)
}




###################################################
####Functions just for bootstrapping
###################################################

#####Theta functions for bootstrap. These get passed to BCa.Boot.CI 
Mean<- function(x,i){mean(x[i],na.rm = TRUE)} 
Median<- function(x,i){median(x[i],na.rm = TRUE)}

Peak = function (x,i)
{ 
  DensityDistro<-density(x[i],adjust=.5, bw = 'Sj',kernel='gaussian',na.rm = TRUE)
  PeakKernel<-DensityDistro$x[which(DensityDistro$y == max(DensityDistro$y))]
  return(PeakKernel)
}

Count <-function(x,i) {
  sum( !is.na( x[i] )) }


### Rindex for Exact values
r.index.calc.boot.Exact = function (p.exact,i)
{ 
  sig.values  <- ifelse(p.exact[i]<.05, 1, 0)
  p.sig <- mean(sig.values, na.rm = TRUE)
  
  zvalues<-qnorm(1-(p.exact[i]/2))
  EstPower <-pnorm((zvalues)-qnorm(0.975))
  med.est.pow <- median(EstPower,na.rm = TRUE)
  
  inf.rate = p.sig - med.est.pow
  r.ind = med.est.pow - inf.rate
  return(r.ind)
}


##########Bootstrap on Index
r.index.calc.boot = function (zvalues,i)
{
  EstPower <-pnorm((zvalues[i])-qnorm(0.975))
  med.est.pow <- median(EstPower,na.rm = TRUE)
  sig.values  <- ifelse((1-pnorm(zvalues[i]))<.05, 1, 0)
  p.sig <- mean(sig.values, na.rm = TRUE)
  
  
  inf.rate = p.sig - med.est.pow
  r.ind = med.est.pow - inf.rate
  return(r.ind)
}


##########Bootstrap on TIVA
tiva.calc = function (x,i)
{
  z.var = var(x[i], na.rm=T)
  z.df = length(x[i])-1
  tiva = z.var*z.df
  tiva.p = (1-pchisq(tiva,z.df))
  return(tiva)
}



###Function for 2000 BCa bootstrap CI. 
BCa.Boot.CI<- function(data,x,theta,LogT=FALSE,splitter, StatType=NULL){ 
  
  require(boot)
  set.seed(42) #increases replibility when other try to reproduce the same values as us 
  
  x <- eval(substitute(x), data)
  
  ###Set number of rows 
  if (missing(splitter)) {Nrows=1}
  else {
    splitter <- eval(substitute(splitter), data) 
    Nrows <- length(levels(splitter))}
  
  out<-matrix(0, Nrows, 4)
  
  
  ##Calcaute bootstrap for each split
  for (i in  1:Nrows)
  { 
    #set data for analysis
    if (missing(splitter)) {y<-x}
    else {y<-subset(x, splitter == levels(splitter)[i]) }
    
    ##tranform y if called
    if (LogT==FALSE) {y<-y}
    else {y<-log10(y)}
    
    Art<-theta(y)
    bootObj<-boot(data=y, statistic=theta, R=2000)
    Bootmean<-mean(bootObj$t)
    #BootSE<-sd(bootObj$t)
    BootCI<-boot.ci(bootObj, type="bca")
    LB<-BootCI$bca[4]
    UB<-BootCI$bca[5]
    out[i,]<-c(Art,Bootmean,LB,UB)
  }#end loop
  
  #name columns
  if (LogT==FALSE) {out<-out}
  else {out<-10^(out)-1}
  
  out<-data.frame(out)
  out$Analysis <- deparse(substitute(theta))
  
  if (missing(splitter)) {  
    names(out) <- c("Arithmetic","BootStat","LB","UP","Analysis")}
  else {  
    out$Split <- levels(splitter)
    names(out) <- c("Arithmetic","BootStat","LB","UP","Analysis","splitter")} 
  
  ##add Analysis Name
  if (missing(splitter)) {  
    orderAll<-c(5,1:4)
    orderBS<-c(5,2:4)
    orderAR<-c(5,1,3:4)
    orderCI<-c(5,3:4)}
  else {
    orderAll<-c(5,1:4,6)
    orderBS<-c(5,2:4,6)
    orderAR<-c(5,1,3:4,6)
    orderCI<-c(5,3:4,6)}
  
  ###Select output
  if (missing(StatType)) { out <- out[orderAll] }
  else if (StatType=='Both') { out <- out[orderAll] }
  else if (StatType=='BS') { out <- out[orderBS]}  
  else if (StatType=='AR') { out <- out[orderAR]}
  else if (StatType=='CI') { out <- out[orderCI]}
  else { out <- out }
  
  return(out)
}

