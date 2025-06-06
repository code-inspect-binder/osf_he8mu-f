
##########################################################################################################
###Script to analyze Status of our Science data Study 2
# 11/14/2016
# Matt Motyl & Alexander Demos, With help of Zach Melton, Alex Sokolovsky, and Carlos Salas 
##########################################################################################################
# Install necessary packages [ RUN ME ONCE< NOT EVERY time]
### Loaded into packrat so you can use the same version
install.packages("plyr")        # data manipulator package
install.packages("pwr")         # Load Power analysis formulas
install.packages("boot")        # Bootstrapping
install.packages("ggplot2")     # GGplot package for visualizing data
install.packages("scales")      # used in ggplot to get log scales to appear
install.packages("gridExtra")   # ggplot grid help
install.packages("reshape2")    # Helper functions
install.packages("knitr")       # Make tables
install.packages("np")          # Nonparametric Kernel Methods for Mixed Datatypes
install.packages("lsr")         # Chi square effect size
install.packages("matrixStats") # Weighted Median
install.packages("psych")       # Used for Rpb correlation


# Load necessary packages 
library(plyr)     
library(pwr)     
library(boot)   
library(ggplot2)
library(scales)   
library(gridExtra)
library(reshape2) 
library(knitr)
library(np)   
library(lsr)
# Note: Other packages load later based on use as to not interfer with each other. 

######################################################################################
######################################################################################
#### Set working directory for data: NOTE: User needs to change directory below  #####
######################################################################################
######################################################################################
setwd("C:/Users/AlexUIC/Dropbox/StatusofScience")

# Make sure to have the custom functions file in the same directory as the working directory
# this function contains our bootstapping function as well as the other useful functions
source("Custom_Functions.R")  # Load custom built function

######################################################################################
######################################################################################
############################ DATA MANAGEMENT #########################################
######################################################################################
######################################################################################

# Load data
ScienceStatus = read.csv("ScienceStatus_C.csv", row.names = 1)

# Rename 
names(ScienceStatus) = c("id","article.id","Keep","timestamp","coder","journal","year","jpsp.section",
                         "first.author","affiliation","specialization","citation",
                         "type","citations.number","study.number","sample.size","excluded.participants",
                         "expcor","design","conditions.number","predictors.number","dv.number",
                         "which.statistical.test","copied.statistic","df.numerator","df.denominator","test.statistic",
                         "p.exact","p.exact.value","drop.column", "null.predicted","sigtests.number.reported",
                         "sigtests.number.significant","total.num.studies","effect.size.reported","effect.size","footnotes.number",
                         "additional.analyses","covariates.number","coding.difficulty","comments")


# Some of the article were double coded. This removes them and leaves 1 coder per article. 
ScienceStatus<-subset(ScienceStatus,Keep==1)

#Added ID to above and covert to a factor to allow binding of subsets for later analysis
ScienceStatus$ID <-as.factor(ScienceStatus$id)

#Create new variables to store recoded data from old variables
#Recoding variables to be numeric... must convert to string before converting to numeric
ScienceStatus$sample.size_r = as.numeric(as.character(ScienceStatus$sample.size))
ScienceStatus$excluded.participants_r = as.numeric(as.character(ScienceStatus$excluded.participants))
ScienceStatus$predictors.number_r = as.numeric(ScienceStatus$predictors.number)
ScienceStatus$copied.stat_r = as.character(ScienceStatus$copied.statistic)
ScienceStatus$p.exact.value_r = as.numeric(as.character(ScienceStatus$p.exact.value))
ScienceStatus$sigtests.number.reported_r = as.numeric(as.character(ScienceStatus$sigtests.number.reported))
ScienceStatus$sigtests.number.significant_r = as.numeric(as.character(ScienceStatus$sigtests.number.significant))
ScienceStatus$footnotes.number_r = as.numeric(as.character(ScienceStatus$footnotes.number))
ScienceStatus$covariates.number_r = as.numeric(as.character(ScienceStatus$covariates.number))
ScienceStatus$year_r = as.factor(ScienceStatus$year)

##convert test statistic variable into numeric variable and store in new variable
ScienceStatus$test.statistic.values <- as.numeric(as.character(ScienceStatus$test.statistic))
ScienceStatus$test.statistic.values.abs <- abs(ScienceStatus$test.statistic.values)

# convert degrees of freedom into numeric variables and store in new variable
ScienceStatus$df.numerator.value = as.numeric(as.character(ScienceStatus$df.numerator))
ScienceStatus$df.denominator.value = as.numeric(as.character(ScienceStatus$df.denominator))

#Assigning numeric values to variable labels which were stored in CSV
ScienceStatus$coding.difficulty_r = revalue(ScienceStatus$coding.difficulty, c("Very Easy"="1",
                                           "Moderately Easy"="2", "Slightly Easy"="3", "Neither Difficult nor Easy"="4",
                                           "Slightly Difficult"="5", "Moderately Difficult"="6", "Very Difficult"="7"))
ScienceStatus$coding.difficulty_r = as.numeric(as.character(ScienceStatus$coding.difficulty_r))

# create two-level year variable
# when running analyses on these data, we'll have more stable estimates from ~10 years ago and the past 2 years
ScienceStatus$yearcat<-ScienceStatus$year_r
ScienceStatus$yearcat <- revalue(ScienceStatus$yearcat, c("2003"="2003-2004","2004"="2003-2004","2013"="2013-2014","2014"="2013-2014"))

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
######## Calculate Power, effect size, pvalue, covert over p to zvalues ############################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


#Create a new data frame in which to subset and rejoin the data later for each test type (t, F, r, chi)
ScienceStatus_J<-NULL

######################
###t tests
#####################
#Clean for T-test
ScienceStatus_Clean_T <- subset(ScienceStatus, which.statistical.test=="t-test (t)" & df.denominator.value !="NA" )
#calc effect size & Pvalue
ScienceStatus_Clean_T$t.test.dest<-dest(ScienceStatus_Clean_T$test.statistic.values,ScienceStatus_Clean_T$df.denominator.value)
ScienceStatus_Clean_T$pval <- pt(-abs(ScienceStatus_Clean_T$test.statistic.values),df=ScienceStatus_Clean_T$df.denominator.value)
ScienceStatus_Clean_T$t.test.r2<-R2fromT(ScienceStatus_Clean_T$test.statistic.values,ScienceStatus_Clean_T$df.denominator.value)



#T-power
ScienceStatus_Clean_T$t.test.Observed.power <-
  pwr.t.test(n = (ScienceStatus_Clean_T$df.denominator.value/2) +1, d = ScienceStatus_Clean_T$t.test.dest, sig.level = 0.05, power = NULL,
             type = c("two.sample"),
             alternative = "two.sided")$power

#Calculate Estimated Power
ScienceStatus_Clean_T$T.pvalue_calc<-pvalsT <- pt(-abs(ScienceStatus_Clean_T$test.statistic.values),df=ScienceStatus_Clean_T$df.denominator.value)
ScienceStatus_Clean_T$EstimatedPower_T<-EstPower(ScienceStatus_Clean_T$T.pvalue_calc)

CountNM(ScienceStatus_Clean_T$t.test.Observed.power)
CountNM(ScienceStatus_Clean_T$EstimatedPower_T)

##Assume small effect calculate N
for (i in 1:length(ScienceStatus_Clean_T)) {
  Route<-(ScienceStatus_Clean_T$design[i]=="Within")
  
  if (Route==TRUE) {
  ScienceStatus_Clean_T$t.test.N[i] <-
    pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
               type = c("paired"),
               alternative = "two.sided")$n
  } else {
    ScienceStatus_Clean_T$t.test.N[i] <-
    pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
               type = c("two.sample"),
               alternative = "two.sided")$n
  }
}

ScienceStatus_Clean_T$t.test.N<-ceiling(ScienceStatus_Clean_T$t.test.N)*2

#Remerge back to larger Dataset - not elagent, but problems from orginal files need to be tracked down.....
ScienceStatus_Clean_T_simple <- ScienceStatus_Clean_T[, c("ID","t.test.Observed.power", "t.test.r2","EstimatedPower_T","T.pvalue_calc","t.test.N")]

ScienceStatus_J<-join(ScienceStatus, ScienceStatus_Clean_T_simple, by = "ID")
rm("ScienceStatus_Clean_T_simple")
rm("ScienceStatus_Clean_T")

######################
###ANOVA
#####################
ScienceStatus_Clean_F <- subset(ScienceStatus, which.statistical.test=="ANOVA (F)" & df.numerator.value!="NA" & df.denominator.value !="NA" & test.statistic.values!="NA" | 
                                  which.statistical.test=="ANCOVA (F)" & df.numerator.value!="NA" & df.denominator.value !="NA" & test.statistic.values!="NA" |
                                  which.statistical.test=="ANCOVA"& df.numerator.value!="NA" & df.denominator.value !="NA" & test.statistic.values!="NA"|
                                  which.statistical.test=="General Linear Model"& df.numerator.value!="NA" & df.denominator.value !="NA" & test.statistic.values!="NA"|
                                  which.statistical.test=="GLM (Ftest)"& df.numerator.value!="NA" & df.denominator.value !="NA" & test.statistic.values!="NA") 



#F-effect size
ScienceStatus_Clean_F$cohen.F2<-cohenFfromN2(ScienceStatus_Clean_F$test.statistic.values,ScienceStatus_Clean_F$df.numerator.value,ScienceStatus_Clean_F$df.denominator.value)
ScienceStatus_Clean_F$F.test.r2<-R2fromF(ScienceStatus_Clean_F$test.statistic.values,ScienceStatus_Clean_F$df.numerator.value,ScienceStatus_Clean_F$df.denominator.value)


#f-power
ScienceStatus_Clean_F$F.test.Observed.power <-
  pwr.f2.test(u = ScienceStatus_Clean_F$df.numerator.value, v = ScienceStatus_Clean_F$df.denominator.value, f2 = ScienceStatus_Clean_F$cohen.F2, sig.level = 0.05, power = NULL)$power

hist(ScienceStatus_Clean_F$F.test.Observed.power)

#Calculate Estimated Power
ScienceStatus_Clean_F$F.pvalue_calc <- pf(ScienceStatus_Clean_F$test.statistic.values,ScienceStatus_Clean_F$df.numerator.value,ScienceStatus_Clean_F$df.denominator.value, lower.tail=F)
ScienceStatus_Clean_F$EstimatedPower_F<-EstPower(ScienceStatus_Clean_F$F.pvalue_calc)

#check Count
CountNM(ScienceStatus_Clean_F$F.test.Observed.power) - CountNM(ScienceStatus_Clean_F$EstimatedPower_F)

##Assume small effect calculate N
ScienceStatus_Clean_F$F.test.N<-NULL  

for (i in 1:length(ScienceStatus_Clean_F)) {
  Route<-(ScienceStatus_Clean_F$design[i]=="Within")
  
  if (Route==TRUE) {
    ScienceStatus_Clean_F$F.test.N[i] <-
      pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
                 type = c("paired"),
                 alternative = "two.sided")$n
  } else {
    ScienceStatus_Clean_F$F.test.N[i] <-
      pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
                 type = c("two.sample"),
                 alternative = "two.sided")$n
  }
}


ScienceStatus_Clean_F$F.test.N<-ceiling(ScienceStatus_Clean_F$F.test.N)*2

#Remerge back to larger Dataset - not elagent, but problems from orginal files need to be tracked down.....
ScienceStatus_Clean_F_simple <- ScienceStatus_Clean_F[, c("ID","F.test.Observed.power", "F.test.r2","EstimatedPower_F","F.pvalue_calc","F.test.N")]
ScienceStatus_J<-join(ScienceStatus_J, ScienceStatus_Clean_F_simple, by = "ID")
rm("ScienceStatus_Clean_F_simple")
rm("ScienceStatus_Clean_F")

######################
###Regression
#####################
#THis is a rought approx.  I cannot correct F/T correctly using the data extracted from the papers.  

ScienceStatus_Clean_rg <- subset(ScienceStatus, which.statistical.test=="Multiple Regression (t -- for individual predictor)" & df.denominator.value !="NA" & df.denominator.value > 1 & test.statistic.values!="NA")

#convert t to d
ScienceStatus_Clean_rg$rg.test.dest<-dest(ScienceStatus_Clean_rg$test.statistic.values,ScienceStatus_Clean_rg$df.denominator.value)
ScienceStatus_Clean_rg$rg.test.r2<-R2fromT(ScienceStatus_Clean_rg$test.statistic.values,ScienceStatus_Clean_rg$df.denominator.value)

#T-power
ScienceStatus_Clean_rg$rg.test.Observed.power <-
  pwr.t.test(n = ScienceStatus_Clean_rg$df.denominator.value, d = ScienceStatus_Clean_rg$rg.test.dest, sig.level = 0.05, power = NULL,
             type = c("two.sample"),
             alternative = "two.sided")$power


hist(ScienceStatus_Clean_rg$rg.test.Observed.power)

#Calculate Estimated Power
ScienceStatus_Clean_rg$rg.pvalue_calc<- pt(-abs(ScienceStatus_Clean_rg$test.statistic.values),df=ScienceStatus_Clean_rg$df.denominator.value)
ScienceStatus_Clean_rg$EstimatedPower_rg<-EstPower(ScienceStatus_Clean_rg$rg.pvalue_calc)

#check Count
CountNM(ScienceStatus_Clean_rg$rg.test.Observed.power) - CountNM(ScienceStatus_Clean_rg$EstimatedPower_rg)


##Assume small effect calculate N
ScienceStatus_Clean_rg$rg.test.N<-NULL  
for (i in 1:length(ScienceStatus_Clean_rg)) {
  Route<-(ScienceStatus_Clean_rg$design[i]=="Within")
  
  if (Route==TRUE) {
    ScienceStatus_Clean_rg$rg.test.N[i] <-
      pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
                 type = c("paired"),
                 alternative = "two.sided")$n
  } else {
    ScienceStatus_Clean_rg$rg.test.N[i] <-
      pwr.t.test(n = NULL, d = .430, sig.level = 0.05, power = .80,
                 type = c("two.sample"),
                 alternative = "two.sided")$n
  }
}

ScienceStatus_Clean_rg$rg.test.N<-ceiling(ScienceStatus_Clean_rg$rg.test.N)*2


#Remerge back to larger Dataset 
ScienceStatus_Clean_rg_simple <- ScienceStatus_Clean_rg[, c("ID","rg.test.Observed.power", "rg.test.r2","EstimatedPower_rg","rg.pvalue_calc","rg.test.N")]

ScienceStatus_J<-join(ScienceStatus_J, ScienceStatus_Clean_rg_simple, by = "ID")
rm("ScienceStatus_Clean_rg_simple")
rm("ScienceStatus_Clean_rg")
######################
###Correlation
#####################

ScienceStatus_Clean_r <- subset(ScienceStatus, which.statistical.test=="Basic Correlation (r)" & df.denominator.value !="NA" & test.statistic.values!="NA"& abs(test.statistic.values) <= 1  )
ScienceStatus_Clean_r$r.test.r2 <-(ScienceStatus_Clean_r$test.statistic.values)^2
summary(ScienceStatus_Clean_r$r.test.r2)

#r-power
ScienceStatus_Clean_r$r.test.Observed.power <- 
  pwr.r.test(n = ScienceStatus_Clean_r$df.denominator.value+2, r = abs(ScienceStatus_Clean_r$test.statistic.values), sig.level = 0.05, power = NULL,
             alternative = "two.sided")$power

hist(ScienceStatus_Clean_r$r.test.Observed.power)

#Calculate Estimated Power
##corrlation to t
ScienceStatus_Clean_r$rtot <-ScienceStatus_Clean_r$test.statistic.values / ((1-ScienceStatus_Clean_r$test.statistic.values^2) / (ScienceStatus_Clean_r$df.denominator.value))^.5
ScienceStatus_Clean_r$r.pvalue_calc <- pt(-abs(ScienceStatus_Clean_r$rtot),df=ScienceStatus_Clean_r$df.denominator.value)
ScienceStatus_Clean_r$EstimatedPower_r<-EstPower(ScienceStatus_Clean_r$r.pvalue_calc)


##Assume small effect calculate N
ScienceStatus_Clean_r$r.test.N<-pwr.r.test(n = NULL, r = .21, sig.level = 0.05, 
                                            power = .80, alternative = "two.sided")$n

ScienceStatus_Clean_r$r.test.N<-ceiling(ScienceStatus_Clean_r$r.test.N)

#Remerge back to larger Dataset 
ScienceStatus_Clean_r_simple <- ScienceStatus_Clean_r[, c("ID","r.test.Observed.power", "r.test.r2","EstimatedPower_r","r.pvalue_calc","r.test.N")]

ScienceStatus_J<-join(ScienceStatus_J, ScienceStatus_Clean_r_simple, by = "ID")
rm("ScienceStatus_Clean_r_simple")
rm("ScienceStatus_Clean_r")

######################
###Chi sqaure
#####################

ScienceStatus_Clean_chi <- subset(ScienceStatus, which.statistical.test=="Chi Square (X)" & df.denominator.value !="NA" & test.statistic.values!="NA" & sample.size_r!="NA")

#Calculate Estimated Power
ScienceStatus_Clean_chi$chi.pvalue_calc<-  pchisq(ScienceStatus_Clean_chi$test.statistic.values, ScienceStatus_Clean_chi$df.denominator.value, lower.tail=FALSE) 
ScienceStatus_Clean_chi$EstimatedPower_chi<-EstPower(ScienceStatus_Clean_chi$chi.pvalue_calc)


## w effect size
ScienceStatus_Clean_chi$Chi.effect<-wchiest(ScienceStatus_Clean_chi$test.statistic.values, ScienceStatus_Clean_chi$sample.size_r)

#r-power
ScienceStatus_Clean_chi$chi.test.Observed.power <- 
  pwr.chisq.test(w = ScienceStatus_Clean_chi$Chi.effect, N = ScienceStatus_Clean_chi$sample.size_r, df = ScienceStatus_Clean_chi$df.denominator.value, sig.level = 0.05, power = NULL)$power

hist(ScienceStatus_Clean_chi$chi.test.Observed.power )


##Assume small effect calculate N
ScienceStatus_Clean_chi$chi.test.N<-NULL

for (i in 1:length(ScienceStatus_Clean_chi)) {
    dfn<-ScienceStatus_Clean_chi$df.denominator.value[i]
    
    ScienceStatus_Clean_chi$chi.test.N[i]<-pwr.chisq.test(w = .2, N = NULL, 
                                                       df = dfn, sig.level = 0.05, power = .8)$N
}

ScienceStatus_Clean_chi$chi.test.N<-ceiling(ScienceStatus_Clean_chi$chi.test.N)


ScienceStatus_Clean_chi_simple <- ScienceStatus_Clean_chi[, c("ID","chi.test.Observed.power", "EstimatedPower_chi","chi.pvalue_calc","chi.test.N")]

ScienceStatus_J<-join(ScienceStatus_J, ScienceStatus_Clean_chi_simple, by = "ID")

rm("ScienceStatus_Clean_chi_simple")
rm("ScienceStatus_Clean_chi")



#######Merge back into one variable 

# Cohen's Oberved power (Post-hoc)
ScienceStatus_J$O.Power<-(rowSums(cbind(ScienceStatus_J$t.test.Observed.power, ScienceStatus_J$F.test.Observed.power, ScienceStatus_J$rg.test.Observed.power, 
                                        ScienceStatus_J$r.test.Observed.power,ScienceStatus_J$chi.test.Observed.power), na.rm = TRUE) + 
                            ifelse(is.na(ScienceStatus_J$t.test.Observed.power) & is.na(ScienceStatus_J$F.test.Observed.power) & 
                                     is.na(ScienceStatus_J$rg.test.Observed.power) & is.na(ScienceStatus_J$r.test.Observed.power) & 
                                     is.na(ScienceStatus_J$chi.test.Observed.power), NA, 0))

# Convert Pvalue to power (Post-hoc)
ScienceStatus_J$E.Power<-(rowSums(cbind(ScienceStatus_J$EstimatedPower_T, ScienceStatus_J$EstimatedPower_F, ScienceStatus_J$EstimatedPower_rg, 
                                        ScienceStatus_J$EstimatedPower_r,ScienceStatus_J$EstimatedPower_chi), na.rm = TRUE) + 
                            ifelse(is.na(ScienceStatus_J$EstimatedPower_T) & is.na(ScienceStatus_J$EstimatedPower_F) & 
                                     is.na(ScienceStatus_J$EstimatedPower_rg) & is.na(ScienceStatus_J$EstimatedPower_r) & 
                                     is.na(ScienceStatus_J$EstimatedPower_chi), NA, 0))

# Calcuate EXACT pvalue from the stats reported in the paper (Post-hoc)
ScienceStatus_J$Calc.Pvalue<-(rowSums(cbind(ScienceStatus_J$T.pvalue_calc, ScienceStatus_J$F.pvalue_calc, ScienceStatus_J$rg.pvalue_calc, 
                                            ScienceStatus_J$r.pvalue_calc,ScienceStatus_J$chi.pvalue_calc), na.rm = TRUE) + 
                                ifelse(is.na(ScienceStatus_J$T.pvalue_calc) & is.na(ScienceStatus_J$F.pvalue_calc) & 
                                         is.na(ScienceStatus_J$rg.pvalue_calc) & is.na(ScienceStatus_J$r.pvalue_calc) & 
                                         is.na(ScienceStatus_J$chi.pvalue_calc), NA, 0))

# Calcuate effect size from the stats reported in the paper (Post-hoc)
ScienceStatus_J$r2.effect<-(rowSums(cbind(ScienceStatus_J$t.test.r2, ScienceStatus_J$F.test.r2, ScienceStatus_J$rg.test.r2, 
                                          ScienceStatus_J$r.test.r2), na.rm = TRUE) + 
                              ifelse(is.na(ScienceStatus_J$t.test.r2) & is.na(ScienceStatus_J$F.test.r2) & 
                                       is.na(ScienceStatus_J$rg.test.r2) & is.na(ScienceStatus_J$r.test.r2), NA, 0))

# Calcuate N that would have been needed to detect a small effect effect (a priori power)
ScienceStatus_J$n.est<-(rowSums(cbind(ScienceStatus_J$t.test.N, ScienceStatus_J$F.test.N, ScienceStatus_J$rg.test.N, 
                                          ScienceStatus_J$r.test.N), na.rm = TRUE) + 
                              ifelse(is.na(ScienceStatus_J$t.test.N) & is.na(ScienceStatus_J$F.test.N) & 
                                       is.na(ScienceStatus_J$rg.test.N) & is.na(ScienceStatus_J$r.test.N), NA, 0))

################################################
###### Covert all pvalues into zscores #########
################################################

ScienceStatus_J$Calc.Pvalue<-pbound(ScienceStatus_J$Calc.Pvalue)

ScienceStatus_J$Sig <- ifelse(ScienceStatus_J$Calc.Pvalue<=.05, 1, 0)

ScienceStatus_J$Significant<-factor(ScienceStatus_J$Sig,
                                    levels=c(0,1),
                                    labels=c("Not Sig","Sig"))

ScienceStatus_J$Calc.Z<-qnorm(1-(ScienceStatus_J$Calc.Pvalue/2))
summary(ScienceStatus_J$Calc.Z)


##########################################
############# Convert pvalues ############
##########################################

#Convert all exact to rounded values
ScienceStatus_J$Round.Pvalue<-sapply(ScienceStatus_J$Calc.Pvalue, Prounder)

#Merge the two group (reported exactly - calculauted VS reported by rounding)
ScienceStatus_J$Report.Pvalue<-ifelse(ScienceStatus_J$p.exact=="Yes",ScienceStatus_J$Round.Pvalue,ScienceStatus_J$Calc.Pvalue)

####APA 6 format
ScienceStatus_J$APA.Pvalue<-sapply(ScienceStatus_J$Calc.Pvalue, PvalueAPA6)


##############################################################
###### Correct Article ID code for by Paper analyses #########
##############################################################

ScienceStatus_J <- ScienceStatus_J[order(ScienceStatus_J$citation, na.last=FALSE) , ]

article.id.corrected<-c(1,rep(0,nrow(ScienceStatus_J)-1))
for (i in 2:nrow(ScienceStatus_J)) {
  article.id.corrected[i]<-ifelse(ScienceStatus_J$citation[i]==ScienceStatus_J$citation[i-1],article.id.corrected[i-1],article.id.corrected[i-1]+1)
  
}

ScienceStatus_J$article.id.corrected<-article.id.corrected
###Correct.ID.errors

Study.Number<-ddply(ScienceStatus_J, .(citation), summarize,
                    StudyN=CountNM(article.id.corrected))

ScienceStatus_J<-merge(Study.Number,ScienceStatus_J, by = "citation")

########## Generate an ID per PAPER
length(unique(ScienceStatus_J$citation))
length(unique(ScienceStatus_J$article.id.corrected))

ScienceStatus_J$article.id<-ScienceStatus_J$article.id.corrected
CountNM(ScienceStatus_J$article.id)



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
##################################     ANALYSIS     ################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

####################################################################################################
################################# Descriptives on Sampling #########################################
####################################################################################################

### Number of Studies 
Counts.Studies.J<-ddply(ScienceStatus_J, .(journal,article.id), summarize,
                         N=max(StudyN,na.rm=TRUE))
Counts.Studies.J

Counts.Studies.J.2<-ddply(Counts.Studies.J, .(journal), summarize,
                           N=sum(N))
Counts.Studies.J.2

sum(Counts.Studies.J.2$N)

### Count number of journals before Removing PS non-social
Counts.Journal<-ddply(ScienceStatus_J, .(journal,article.id), summarize,
                      N=ifelse(CountNM(StudyN)>0,1,0))
sum(Counts.Journal$N)

Counts.Journal.2<-ddply(Counts.Journal, .(journal), summarize,
                      N=sum(N))
Counts.Journal.2
sum(Counts.Journal.2$N)






##############################################################################################################
################################################# SAMPLE SELECTOR ############################################
##############################################################################################################

# For the Article we removed all Psych Science papers where the first author did not specialize in Social / Personality

ScienceStatus_J$SP<-ifelse(ScienceStatus_J$journal!="Psych Science",1, ifelse(ScienceStatus_J$journal=="Psych Science" & ScienceStatus_J$specialization=="Social / Personality",1,0))
ScienceStatus_SP<-subset(ScienceStatus_J, SP==1)

### If you wish to overwrite our selection and choose your own adventure, plus select from among the code below.
### Please uncomment out the code you wish to select run it and continue your analysis below.
### Note if you skip this, your analysis will match our analysis. 
### While we did set a seed for out bootstrapping your final numbers may difer slightly in a small way from ours. 

########################################
################ Adventure 1: By Journal
########################################

## JESP: uncomment next line only
#ScienceStatus_SP<-subset(ScienceStatus_J, journal=="JESP")

## JPSP: uncomment next line only
#ScienceStatus_SP<-subset(ScienceStatus_J, journal=="JPSP")

## PSPB: uncomment next line only
#ScienceStatus_SP<-subset(ScienceStatus_J, journal=="PSPB")

## Psych Science: uncomment next line only
#ScienceStatus_SP<-subset(ScienceStatus_J, journal=="Psych Science")

########################################
################ Adventure 2: Only the the articles we removed from psych science
########################################

#ScienceStatus_SP<-subset(ScienceStatus_J$journal=="Psych Science" & ScienceStatus_J$specialization!="Social / Personality")

##############################################################################################################
########################################## Analysis resumes below ############################################
##############################################################################################################

### Count number of Papers relative to full sample
CountNM(ScienceStatus_SP$article.id)

Counts.Journal.SP<-ddply(ScienceStatus_SP, .(journal,article.id), summarize,
                      N=ifelse(CountNM(StudyN)>0,1,0))
sum(Counts.Journal.SP$N)


Counts.Journal.SP.N<-ddply(Counts.Journal.SP, .(journal), summarize,
                      N=sum(N))
Counts.Journal.SP.N

(Counts.Journal.SP.N$N/Counts.Journal.2$N)*100

### Count number of studies relative to full sample
Counts.Studies.SP<-ddply(ScienceStatus_SP, .(journal,article.id), summarize,
                         N=max(StudyN,na.rm=TRUE))
Counts.Studies.SP

Counts.Studies.SP.2<-ddply(Counts.Studies.SP, .(journal), summarize,
                         N=sum(N))
Counts.Studies.SP.2

sum(Counts.Studies.SP.2$N)

## Both years
summary(ScienceStatus_SP$expcor)
summary(ScienceStatus_SP$design)


##########################################################
########## Pvalue difference of reporting styles #######
##########################################################
P.Exact.Calcuated.All<-BCa.Boot.CI(ScienceStatus_SP,Calc.Pvalue,Mean,LogT=FALSE,splitter=p.exact, StatType="AR")
P.Rounded.Calcuated.All<-BCa.Boot.CI(ScienceStatus_SP,Round.Pvalue,Mean,LogT=FALSE,splitter=p.exact, StatType="AR")
P.Mixed.Calcuated.All<-BCa.Boot.CI(ScienceStatus_SP,Report.Pvalue,Mean,LogT=FALSE,splitter=p.exact, StatType="AR")
P.APA.Calcuated.All<-BCa.Boot.CI(ScienceStatus_SP,APA.Pvalue,Mean,LogT=FALSE,splitter=p.exact, StatType="AR")

Descriptives.Pvalue<-rbind(P.Exact.Calcuated.All,P.Rounded.Calcuated.All,P.Mixed.Calcuated.All,P.APA.Calcuated.All)

kable(Descriptives.Pvalue[Descriptives.Pvalue$splitter=="No",], digits=4)
kable(Descriptives.Pvalue[Descriptives.Pvalue$splitter=="Yes",], digits=4)

P.Exact.Calcuated.All.NS<-BCa.Boot.CI(ScienceStatus_SP,Calc.Pvalue,Mean,LogT=FALSE, StatType="AR")
P.Rounded.Calcuated.All.NS<-BCa.Boot.CI(ScienceStatus_SP,Round.Pvalue,Mean,LogT=FALSE, StatType="AR")
P.Mixed.Calcuated.All.NS<-BCa.Boot.CI(ScienceStatus_SP,Report.Pvalue,Mean,LogT=FALSE, StatType="AR")
P.APA.Calcuated.All.NS<-BCa.Boot.CI(ScienceStatus_SP,APA.Pvalue,Mean,LogT=FALSE, StatType="AR")

Descriptives.Pvalue.NS<-rbind(P.Exact.Calcuated.All.NS,P.Rounded.Calcuated.All.NS,P.Mixed.Calcuated.All.NS,P.APA.Calcuated.All.NS)

kable(Descriptives.Pvalue.NS, digits=4)


####################################################################################################
################################# Table 3 Metrics ##################################################
####################################################################################################

##Sig tests
Median.sigtests.number.reported<-BCa.Boot.CI(ScienceStatus_SP,sigtests.number.reported_r,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Mean.sigtests.number.reported<-BCa.Boot.CI(ScienceStatus_SP,sigtests.number.reported_r,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

ScienceStatus_SP$sigtests.number.reported_r2<-ScienceStatus_SP$sigtests.number.reported_r
ScienceStatus_SP$sigtests.number.reported_r2[ScienceStatus_SP$sigtests.number.reported_r2==0] <- NA
ScienceStatus_SP$within.study.sig<-(ScienceStatus_SP$sigtests.number.significant_r/ScienceStatus_SP$sigtests.number.reported_r2)*100
within.study.sig<-BCa.Boot.CI(ScienceStatus_SP,within.study.sig,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
Sig.Crit<-BCa.Boot.CI(ScienceStatus_SP,Sig*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

#Reporting practices
effect.size.Per<-BCa.Boot.CI(ScienceStatus_SP,effect.size.reported*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

ScienceStatus_SP$Pexact<-ifelse(ScienceStatus_SP$p.exact=="Yes",1,0)
ExactPvalue<-BCa.Boot.CI(ScienceStatus_SP,Pexact*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

#Rounding
ScienceStatus_SP$Rounding<-ScienceStatus_SP$p.exact.value_r - ScienceStatus_SP$Calc.Pvalue
ScienceStatus_SP$RoundingDown1<-ifelse(ScienceStatus_SP$Rounding < 0,1,0)
ScienceStatus_SP$RoundingDown2<-ifelse(ScienceStatus_SP$Rounding < -.004,1,0)
ExactPvalueRoundDown1<-BCa.Boot.CI(ScienceStatus_SP,RoundingDown1*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
ExactPvalueRoundDown2<-BCa.Boot.CI(ScienceStatus_SP,RoundingDown2*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

#Excluded
ScienceStatus_SP$excluded.participants_Percent<-(ScienceStatus_SP$excluded.participants_r/ScienceStatus_SP$sample.size_r)*100
excluded.part.Per<-BCa.Boot.CI(ScienceStatus_SP,excluded.participants_Percent,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

ScienceStatus_SP$excluded.participants_Studies<-(ScienceStatus_SP$excluded.participants_r>0)/(ScienceStatus_SP$excluded.participants_r>=0)*100
excluded.Studies.Per<-BCa.Boot.CI(ScienceStatus_SP,excluded.participants_Studies,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

#Footnotes + other reports
footnotes.Per<-BCa.Boot.CI(ScienceStatus_SP,footnotes.number_r,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")


ScienceStatus_SP$additional.analyses_r <-ifelse(ScienceStatus_SP$additional.analyses=="Yes",1,0)
additional.analyses.data<-ddply(ScienceStatus_SP, .(article.id,yearcat), summarize,.drop=TRUE,
                                additional.analyses_r = max(additional.analyses_r,na.rm = FALSE))

additional.analyses<-BCa.Boot.CI(additional.analyses.data,additional.analyses_r*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")


###Number of studies per paper
total.num.studies.data<-ddply(ScienceStatus_SP, .(article.id,yearcat), summarize,.drop=TRUE,
                              total.num.studies = max(total.num.studies,na.rm = FALSE))

total.num.studies.Mean<-BCa.Boot.CI(total.num.studies.data,total.num.studies,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

predictors.number<-BCa.Boot.CI(ScienceStatus_SP,as.numeric(predictors.number),Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
conditions.number<-BCa.Boot.CI(ScienceStatus_SP,as.numeric(conditions.number),Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
covariates.number<-BCa.Boot.CI(ScienceStatus_SP,as.numeric(covariates.number),Mean,LogT=FALSE,splitter=yearcat, StatType="AR")


###Tables
Descriptives<-rbind(Median.sigtests.number.reported,
                    Mean.sigtests.number.reported,
                    within.study.sig,
                    Sig.Crit,
                    ExactPvalue,
                    ExactPvalueRoundDown1,
                    ExactPvalueRoundDown2,
                    effect.size.Per,
                    excluded.Studies.Per,
                    excluded.part.Per,
                    footnotes.Per,
                    additional.analyses,
                    total.num.studies.Mean,
                    predictors.number,
                    conditions.number,
                    covariates.number)

kable(Descriptives[Descriptives$splitter=="2003-2004",], digits=2)
kable(Descriptives[Descriptives$splitter=="2013-2014",], digits=2)


####################################################################################################
####################################################################################################
####################################################################################################
################################### Replicability Metrics ##########################################
####################################################################################################
####################################################################################################
####################################################################################################

########################################################################### 
###########################################################################  
################################  Pcurve Analysis #########################
###########################################################################  
########################################################################### 

##########################################################
######### old Version of Pcurve Type Analysis
##########################################################

#Per study
Calcuated.atp04.study<-ddply(ScienceStatus_SP, .(yearcat), summarize,
                           less.than.04 = sum(Calc.Pvalue<.04, na.rm=TRUE),
                           less.than.05 = sum(Calc.Pvalue>=.04 &  Calc.Pvalue<.05, na.rm=TRUE))

Calcuated.atp04.study
Calcuated.atp04.study.Bitest.03<-binom.test(as.numeric(Calcuated.atp04.study[1,2:3]),alternative = c("two.sided"))
Calcuated.atp04.study.Bitest.03
Calcuated.atp04.study.Bitest.13<-binom.test(as.numeric(Calcuated.atp04.study[2,2:3]),alternative = c("two.sided"))
Calcuated.atp04.study.Bitest.13
chisq.test(Calcuated.atp04.study[1:2,2:3],simulate.p.value = FALSE)
cramersV(Calcuated.atp04.study[1:2,2:3] )

Calcuated.atp04.study.P<-Calcuated.atp04.study
Calcuated.atp04.study.P[1,2:3]<-(Calcuated.atp04.study.P[1,2:3]/as.numeric(sum(Calcuated.atp04.study.P[1,2:3])))*100
Calcuated.atp04.study.P[2,2:3]<-(Calcuated.atp04.study.P[2,2:3]/as.numeric(sum(Calcuated.atp04.study.P[2,2:3])))*100
Calcuated.atp04.study.P


#Per paper
Calcuated.med.pvalue<-ddply(ScienceStatus_SP, .(article.id,yearcat), summarize,
                             med.pvalue = median(Calc.Pvalue, na.rm=TRUE))

Calcuated.atp04.paper<-ddply(Calcuated.med.pvalue, .(yearcat), summarize,
                             less.than.04 = sum(med.pvalue<.04, na.rm=TRUE),
                             less.than.05 = sum(med.pvalue>=.04 &  med.pvalue<.05, na.rm=TRUE))
Calcuated.atp04.paper
Calcuated.atp04.paper.Bitest.03<-binom.test(as.numeric(Calcuated.atp04.paper[1,2:3]),alternative = c("two.sided"))
Calcuated.atp04.paper.Bitest.03
Calcuated.atp04.paper.Bitest.13<-binom.test(as.numeric(Calcuated.atp04.paper[2,2:3]),alternative = c("two.sided"))
Calcuated.atp04.paper.Bitest.13
chisq.test(Calcuated.atp04.paper[1:2,2:3],simulate.p.value = FALSE)
cramersV(Calcuated.atp04.paper[1:2,2:3] )

Calcuated.atp04.paper.P<-Calcuated.atp04.paper
Calcuated.atp04.paper.P[1,2:3]<-(Calcuated.atp04.paper.P[1,2:3]/as.numeric(sum(Calcuated.atp04.paper.P[1,2:3])))*100
Calcuated.atp04.paper.P[2,2:3]<-(Calcuated.atp04.paper.P[2,2:3]/as.numeric(sum(Calcuated.atp04.paper.P[2,2:3])))*100
Calcuated.atp04.paper.P


############ Bootstrap % evidentiary per paper Orginal P-curve
Calcuated.med.pvalue$EV<-ifelse(Calcuated.med.pvalue$med.pvalue>.05,NA,ifelse(Calcuated.med.pvalue$med.pvalue<.04,1,0))
PerEvidence.Year.Orginal<-BCa.Boot.CI(Calcuated.med.pvalue,EV*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
PerEvidence.Year.Orginal



##########################################################
######### New Version of Pcurve Analysis
##########################################################

################ Binomal testing

#Per study Calculations
Calcuated.ptest.chi<-ddply(ScienceStatus_SP, .(yearcat), summarize,
                           less.than.025 = sum(Calc.Pvalue<.025, na.rm=T),
                           less.than.05 = sum(Calc.Pvalue>=.025 &  Calc.Pvalue<.05, na.rm=T))

#One ways
###2003-2004
Calcuated.ptest.chi[1,2:3]/sum(Calcuated.ptest.chi[1,2:3])
PerStudy.03.Pcurve<-binom.test(as.numeric(Calcuated.ptest.chi[1,2:3]),alternative = c("two.sided"))


###2013-2014
Calcuated.ptest.chi[2,2:3]/sum(Calcuated.ptest.chi[2,2:3])
PerStudy.13.Pcurve<-binom.test(as.numeric(Calcuated.ptest.chi[2,2:3]),alternative = c("two.sided"))

###Two-way
chisq.test(Calcuated.ptest.chi[1:2,2:3],simulate.p.value = FALSE)
cramersV(Calcuated.ptest.chi[1:2,2:3] )
sum(Calcuated.ptest.chi[1:2,2:3] )
Calcuated.atp04


###########################################################################  
#### Calcuate Pcurve using Pcurve Scripts based on CALCULATED PER PAPER ###
###########################################################################  

#clean pvalues that were non-significant (as per Pcurve)
ScienceStatus_SP$Calc.Pvalue.Clean<-ifelse(ScienceStatus_SP$Calc.Pvalue>=0.05, NA, ScienceStatus_SP$Calc.Pvalue)


#Pcurve by paper - Copies from pcurve website on 10/10/2016
Zppr=NULL
Zppr.half=NULL
p.Zppr=NULL
p.Zppr.half=NULL

for (i in 1:max(ScienceStatus_SP$article.id, na.rm=TRUE)) {
    p<-as.numeric(subset(ScienceStatus_SP,ScienceStatus_SP$article.id==i, select=c(Calc.Pvalue.Clean))$Calc.Pvalue.Clean)
  
    if (sum(is.na(p)==TRUE)==0) {
    
    #2.1 Right Skew, Full p-curve
    ppr=as.numeric(ifelse(p<.05,20*p,NA))            #If p<.05, ppr is 1/alpha*p-value, so 20*pvalue, otherwise missing. 
    ppr=pbound(ppr)                                #apply pbound function to avoid 0
   #2.2 Right Skew, half p-curve
    ppr.half=as.numeric(ifelse(p<.025,40*p,NA))    #If p<.05, ppr is 40*pvalue, otherwise missing. 
    ppr.half=pbound(ppr.half)
    #3.1 Convert pp-values to Z scores, using Stouffer function above
    Zppr[i] = stouffer(ppr)                        #right skew  - this is a Z value from Stouffer's test
    Zppr.half[i]= stouffer(ppr.half)              #right skew, half p-curve - idem 
    #3.2 Overall p-values from Stouffer test
    p.Zppr[i] =pnorm(Zppr[i])	
    p.Zppr.half[i] =pnorm(Zppr.half[i])
  } else {
    Zppr[i] = NA        
    Zppr.half[i]= NA     
    p.Zppr[i] =NA	
    p.Zppr.half[i] =NA
  }
  
}

########## Is there evidentiary value per paper
Evidence<-ifelse(p.Zppr < .05 & p.Zppr.half < .05, 1,0)
article.id<-1:max(ScienceStatus_SP$article.id, na.rm=TRUE)
PcurveE<-as.data.frame(cbind(article.id,Evidence))

ScienceStatus_SP<-merge(PcurveE,ScienceStatus_SP, by = "article.id")
PerEvidence.ID<-ddply(ScienceStatus_SP, .(article.id,yearcat), summarize,
                      PerEvidence = mean(Evidence, na.rm=T),
                      StudyN = mean(StudyN, na.rm=T))
PerEvidence.ID

############ Bootstrap % evidentiary per paper
PerEvidence.Year<-BCa.Boot.CI(PerEvidence.ID,PerEvidence,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
PerEvidence.Year


#################################
## Check if pcurve is related to number of studies reported.
#################################

require(psych)
PerEvidence.ID.C<-subset(PerEvidence.ID, PerEvidence>=0)
biserial(PerEvidence.ID.C$StudyN,PerEvidence.ID.C$PerEvidence)

rbs <- function(d, i){
  d2 <- d[i,]
  return(biserial(d2$StudyN, d2$PerEvidence))
}

bootcorr <- boot(PerEvidence.ID.C[3:4], rbs, 2000)
bootcorr
boot.ci(bootcorr, type = "bca")

PerEvidence.ID.C.03<-PerEvidence.ID.C[PerEvidence.ID.C$yearcat=="2003-2004",]
PerEvidence.ID.C.13<-PerEvidence.ID.C[PerEvidence.ID.C$yearcat=="2013-2014",]
biserial(PerEvidence.ID.C.03$StudyN,PerEvidence.ID.C.03$PerEvidence)
biserial(PerEvidence.ID.C.13$StudyN,PerEvidence.ID.C.13$PerEvidence)
detach("package:psych", unload=TRUE)


########################################################
########################################################
#################### PCurve Graph ######################
########################################################
########################################################

#######################################################
## Calculated from statistic and DF
#######################################################

######### By Paper
P.By.Paper<-ddply(ScienceStatus_SP, .(article.id,yearcat), summarize,
                  Pmean.calc = stouffer.P(Calc.Pvalue.Clean),
                  PMedian.calc = median(Calc.Pvalue.Clean))

#By mean
Calc.pcurve.paper<-ddply(P.By.Paper, .(article.id,yearcat), summarize,
                         less.than.01 = sum(Pmean.calc<.01, na.rm=T),
                         less.than.02 = sum(Pmean.calc>=.01 &  Pmean.calc<.02, na.rm=T),
                         less.than.03 = sum(Pmean.calc>=.02 &  Pmean.calc<.03, na.rm=T),
                         less.than.04 = sum(Pmean.calc>=.03 &  Pmean.calc<.04, na.rm=T),
                         less.than.05 = sum(Pmean.calc>=.04 &  Pmean.calc<.05, na.rm=T))

Calc.pcurve.paper.mean<-cbind(Calc.pcurve.paper$article.id,Calc.pcurve.paper$yearcat, Calc.pcurve.paper[3:7]/rowSums(Calc.pcurve.paper[3:7]))
colnames(Calc.pcurve.paper.mean)[2] <- "yearcat"

Calc.pcurve.paper.percent<-ddply(Calc.pcurve.paper.mean, .(yearcat), summarize,
                                 less.than.01 = mean(less.than.01,na.rm=TRUE),
                                 less.than.02 = mean(less.than.02,na.rm=TRUE),
                                 less.than.03 = mean(less.than.03,na.rm=TRUE),
                                 less.than.04 = mean(less.than.04,na.rm=TRUE),
                                 less.than.05 = mean(less.than.05,na.rm=TRUE))
Calc.pcurve.paper.percent


######### By Median
Calc.pcurve.paper.med<-ddply(P.By.Paper, .(article.id,yearcat), summarize,
                         less.than.01 = sum(PMedian.calc<.01, na.rm=T),
                         less.than.02 = sum(PMedian.calc>=.01 &  PMedian.calc<.02, na.rm=T),
                         less.than.03 = sum(PMedian.calc>=.02 &  PMedian.calc<.03, na.rm=T),
                         less.than.04 = sum(PMedian.calc>=.03 &  PMedian.calc<.04, na.rm=T),
                         less.than.05 = sum(PMedian.calc>=.04 &  PMedian.calc<.05, na.rm=T))

Calc.pcurve.paper.median<-cbind(Calc.pcurve.paper.med$article.id,Calc.pcurve.paper.med$yearcat, Calc.pcurve.paper.med[3:7]/rowSums(Calc.pcurve.paper.med[3:7]))
colnames(Calc.pcurve.paper.median)[2] <- "yearcat"

Calc.pcurve.paper.percent.Median<-ddply(Calc.pcurve.paper.median, .(yearcat), summarize,
                                 less.than.01 = mean(less.than.01,na.rm=TRUE),
                                 less.than.02 = mean(less.than.02,na.rm=TRUE),
                                 less.than.03 = mean(less.than.03,na.rm=TRUE),
                                 less.than.04 = mean(less.than.04,na.rm=TRUE),
                                 less.than.05 = mean(less.than.05,na.rm=TRUE))

## Melt results prepare for graph

ByMean.Pcurve.Calc <- melt(Calc.pcurve.paper.percent)
ByMedian.Pcurve.Calc <- melt(Calc.pcurve.paper.percent.Median)

ByMean.Pcurve.Calc$P.Type<-"By Mean P-curve"
ByMedian.Pcurve.Calc$P.Type<-"By Median P-curve"


pcurve.Graph<-rbind(ByMean.Pcurve.Calc,ByMedian.Pcurve.Calc)

##By Type of pvalue
PcurveGraph<-ggplot(pcurve.Graph, aes(variable, value*100, group=P.Type)) + 
  facet_grid( ~ yearcat)+
  geom_point(aes(color=P.Type,shape=P.Type),size=1.75)+
  geom_line(aes(color=P.Type, linetype=P.Type), size=.75)+
  geom_hline(yintercept = 20,colour = "black",linetype="dashed", size=.5)+
  xlab("P-values")+ylab('Percentage of P-values')+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  scale_x_discrete(labels=c(".01",".02",".03",".04",".05"))+
  scale_color_manual(values=c("black","gray"))+
  annotate("text", x = 4, y = 23, label = "Null of zero effect",size=2)+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
PcurveGraph

ggsave(PcurveGraph,file="Figure5.tiff",scale = 1, dpi = 300,width=4, height=2.5,units = c("in"))


#####################################################################################
############################ Z Curve ##############################################
#####################################################################################

##Calculate n 
CountNM(ScienceStatus_SP$Calc.Z)/ length(ScienceStatus_SP$Calc.Z)

###Calcuate Mean, Median, Peak
Z.Mean.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
Z.Median.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Z.Mode.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,Peak,LogT=FALSE,splitter=yearcat, StatType="AR")

#########
#### Table 
########
Z.Forest.UnTrans<-rbind(Z.Mean.UnTrans,Z.Median.UnTrans,Z.Mode.UnTrans)

kable(Z.Forest.UnTrans, digits=2)

#########
#### Forest Plots 
#########
Z.Forest.Plot<-ggplot(Z.Forest.UnTrans, aes(Arithmetic, Analysis, colour = splitter, shape = splitter))+
  geom_point(size=3) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .2),size = 1,alpha=.65)+
  ylab("")+xlab('Z-Score')+
  #scale_x_continuous(breaks = seq(0, .25, .05),limits=c(0,.25))+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
Z.Forest.Plot

#########
#####Density Plot 
#########

Z.plot<-ggplot(data = ScienceStatus_SP, aes(Calc.Z,..scaled..,fill=yearcat,color=yearcat))+
  facet_wrap( ~ yearcat, ncol=2)+
  coord_cartesian(ylim = c(0,1))+
  stat_density(geom = "area", position = "stack",adjust=.5, bw = 'Sj',kernel='gaussian',trim = FALSE,alpha=.65,na.rm=TRUE)+
  ylab("Density")+
  #xlab('Z-Score')+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none") 
Z.plot

Z.grid<-grid.arrange(Z.plot,Z.Forest.Plot, nrow=2, ncol=1,heights = c(3, 4))
ggsave(Z.grid,file="Figure6.tiff",dpi=300,scale = 1, width=4, height=4,units = c("in"))


### Test denesity plots differences
x <- subset(ScienceStatus_SP, yearcat =="2003-2004" & Calc.Z!="NA", select =Calc.Z  )
y <- subset(ScienceStatus_SP, yearcat =="2013-2014" & Calc.Z!="NA", select =Calc.Z  )
z.D.Entropy<-npunitest(as.vector(x$Calc.Z),as.vector(y$Calc.Z),boot.num=200)
z.D.Entropy


########################################################################### 
###########################################################################  
################################  RIndex ##################################
###########################################################################  
########################################################################### 

###Calcuate Mean, Median, Peak
Rindex.Reults<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,r.index.calc.boot,LogT=FALSE,splitter=yearcat, StatType="AR")
Rindex.Reults

########################################################################### 
###########################################################################  
################################  TIVA ####################################
###########################################################################  
########################################################################### 


TIVA.results<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,tiva.calc,LogT=FALSE,splitter=yearcat, StatType="AR")
TIVA.results

### Normalize TIVA to make it power for Figure 9. 
TIVA.power<-function (x,i)
{
  SSNum = var(x[i], na.rm=T)
  SSDen = var(runif(length(x[i]), min=min(x[i],na.rm=T),max=max(x[i],na.rm=T)))
  
  R2<-SSNum/(SSDen+SSNum)
  cohend<-r2toD(R2)
  TivaP<-pwr.t.test(n = 84, d = cohend, sig.level = 0.05, power = NULL,
                    type = c("two.sample"),
                    alternative = "two.sided")$n
  return(TivaP)
}


TIVA.N.results<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,TIVA.power,LogT=FALSE,splitter=yearcat, StatType="AR")
TIVA.N.results


###################################################################################################
###################################################################################################  
#########################  Metrics for power: effect size, sample size ############################
###################################################################################################  
################################################################################################### 

############################################################################
####### A prior power - Study by Study assuming small effect size ##########
############################################################################
ScienceStatus_SP$Sample.N.test<-ifelse((ScienceStatus_SP$sample.size_r >= ScienceStatus_SP$n.est)==TRUE, 1,
       ifelse((ScienceStatus_SP$sample.size_r > ScienceStatus_SP$n.est)==FALSE,0,NA))


PerStudyAbove.D<-ddply(ScienceStatus_SP, .(yearcat,design), summarize,
                  PerStudyAbove= mean(Sample.N.test,na.rm=TRUE),
                  n = CountNM(Sample.N.test))

PerStudyAbove.J<-ddply(ScienceStatus_SP, .(yearcat,journal), summarize,
                     PerStudyAbove= mean(Sample.N.test,na.rm=TRUE),
                     n = CountNM(Sample.N.test))


P.Mean<-BCa.Boot.CI(ScienceStatus_SP,Sample.N.test,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
P.Mean     

### Based on N

ScienceStatus_SP$Sample.N.test.Scraped<-ifelse((ScienceStatus_SP$sample.size_r >= 172)==TRUE, 1,
                                      ifelse((ScienceStatus_SP$sample.size_r > 172)==FALSE,0,NA))

PerStudyAbove.Scraped<-ddply(ScienceStatus_SP, .(yearcat), summarize,
                     PerStudyAbove= mean(Sample.N.test.Scraped,na.rm=TRUE))

CountNM(ScienceStatus_SP$Sample.N.test.Scraped)

P.Mean.scraped.all<-BCa.Boot.CI(ScienceStatus_SP,Sample.N.test.Scraped,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

P.Mean.scraped.matched<-BCa.Boot.CI(subset(ScienceStatus_SP, Sample.N.test!="NA"),Sample.N.test.Scraped,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")

P.Mean.scraped.all
P.Mean.scraped.matched



############################################################################
#### Analysis on Sample Size
############################################################################

CountNM(ScienceStatus_SP$sample.size_r)
min(ScienceStatus_SP$sample.size_r, na.rm = TRUE)
max(ScienceStatus_SP$sample.size_r, na.rm = TRUE)

###Calcuate Mean, Median, Peak
SS.Mean.Trans<-BCa.Boot.CI(ScienceStatus_SP,sample.size_r,Mean,LogT=TRUE,splitter=yearcat, StatType="AR")
SS.Median.Trans<-BCa.Boot.CI(ScienceStatus_SP,sample.size_r,Median,LogT=TRUE,splitter=yearcat, StatType="AR")
SS.Mode.Trans<-BCa.Boot.CI(ScienceStatus_SP,sample.size_r,Peak,LogT=TRUE,splitter=yearcat, StatType="AR")


#########
#### Table for Sample Size
########
SS.Forest.Trans<-rbind(SS.Mean.Trans,SS.Median.Trans,SS.Mode.Trans)

kable(SS.Forest.Trans, digits=2)

#########
#### Forest Plots on Sample Size
########

SS.Forest.Plot<-ggplot(SS.Forest.Trans, aes(Arithmetic, Analysis, colour = splitter, shape = splitter))+
  geom_point(size=3) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .2),size = 1,alpha=.65)+
  ylab("")+xlab('Sample Size')+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
SS.Forest.Plot

#########
#####Density Plot for sample size
#########

SS.plot.Den<-ggplot(data = ScienceStatus_SP, aes(log10(sample.size_r),..scaled..,fill=yearcat,color=yearcat))+
  facet_wrap( ~ yearcat, ncol=2)+
  stat_density(geom = "area", position = "stack",adjust=.5, bw = 'Sj',kernel='gaussian',trim = FALSE,alpha=.65,na.rm=TRUE)+
  ylab("Density")+
  xlab(expression(Sample~Size~Log[10]))+
  scale_x_continuous(labels = math_format(10^.x))+ 
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none")  +
  annotation_logticks(sides = "b")
SS.plot.Den

SS.grid<-grid.arrange(SS.plot.Den,SS.Forest.Plot, nrow=2, ncol=1,heights = c(3, 4))
ggsave(SS.grid,file="Figure7.tiff",scale = 1, dpi = 300,width=4, height=4,units = c("in"))


### Test denesity plots differences
x <- subset(ScienceStatus_SP, yearcat =="2003-2004" & sample.size_r!="NA", select =sample.size_r  )
y <- subset(ScienceStatus_SP, yearcat =="2013-2014" & sample.size_r!="NA", select =sample.size_r  )
x <-log10(x)
y <-log10(y)

SS.D.Entropy<-npunitest(as.vector(x$sample.size_r),as.vector(y$sample.size_r),boot.num=200)
SS.D.Entropy



#####################################################################################
#####Observed Power
#####################################################################################

##Calculate n 
CountNM(ScienceStatus_SP$O.Power)/ length(ScienceStatus_SP$E.Power)

###Calcuate Mean, Median, Peak
Op.Mean.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,O.Power,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
Op.Median.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,O.Power,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Op.Mode.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,O.Power,Peak,LogT=FALSE,splitter=yearcat, StatType="AR")

#########
#### Table for Post-Hoc Power
########
OP.Forest.Trans<-rbind(Op.Mean.UnTrans,Op.Median.UnTrans,Op.Mode.UnTrans)

kable(OP.Forest.Trans, digits=2)

#########
#### Forest Plots on Post-Hoc Power
########

Op.Forest.Plot<-ggplot(OP.Forest.Trans, aes(Arithmetic, Analysis, colour = splitter, shape = splitter))+
  geom_point(size=3) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .2),size = 1,alpha=.65)+
  ylab("")+xlab('Observed Power')+
  scale_x_continuous(breaks = seq(.7, 1, .05),limits=c(.7,1))+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
Op.Forest.Plot

#########
#####Density Plot on observed Power
#########

Power.plot.Observed<-ggplot(data = ScienceStatus_SP, aes(O.Power,..scaled..,fill=yearcat,color=yearcat))+
  facet_wrap( ~ yearcat, ncol=2)+
  coord_cartesian(ylim = c(0,1))+
  stat_density(geom = "area", position = "stack",adjust=.5, bw = 'Sj',kernel='gaussian',trim = FALSE,alpha=.65,na.rm=TRUE)+
  ylab("Density")+xlab('Observed Power')+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none") 
Power.plot.Observed

OP.grid<-grid.arrange(Power.plot.Observed,Op.Forest.Plot, nrow=2, ncol=1,heights = c(3, 4))
ggsave(OP.grid,file="Figure8.tiff",dpi=300,scale = 1, width=4, height=4,units = c("in"))


### Test denesity plots differences
x <- subset(ScienceStatus_SP, yearcat =="2003-2004" & O.Power!="NA", select =O.Power  )
y <- subset(ScienceStatus_SP, yearcat =="2013-2014" & O.Power!="NA", select =O.Power  )

OP.D.Entropy<-npunitest(as.vector(x$O.Power),as.vector(y$O.Power),boot.num=200)
OP.D.Entropy


########################################################################
########################### Figure 9 ##################################
########################################################################

Figure9Plot = read.csv("Figure9Data.csv")

Figure9Plot$Metric<-factor(Figure9Plot$Metric.1, 
                  levels=c(1,2,3,4,5,6,7,8),
                  labels = c("Ambitious\nP-Curve","Orginal\nP-Curve",
                             "N-Pact","Estimated\n A priori\nPower\nfor d = .43","Observed\nPost hoc\nPower","Z-Curve","R-Index","TIVA"))

FP.Plot<-ggplot(Figure9Plot, aes(x=reorder(Metric, -EV), y=EV, fill=yearcat,group=yearcat)) + 
  coord_cartesian(ylim = c(.25,1))+
  geom_bar(stat = "identity", position="dodge")+
  geom_errorbar(aes(ymin=LB, ymax=UB),position = position_dodge(width = 0.90),width=.2)+
  geom_hline(yintercept = .8)+
  ylab("Evidentiary Value in Power Units")+
  xlab("")+
  scale_fill_manual(values=c("gray40", "gray80"))+
  theme_bw()+
  theme(text = element_text(size=9),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
FP.Plot

ggsave(FP.Plot,file="Figure9.tiff",scale = 1, dpi = 300,width=5.5, height=4,units = c("in"))





####################################################################################
################################# Table 4 ##########################################
####################################################################################


SS.Median.Trans.SP<-BCa.Boot.CI(ScienceStatus_SP,sample.size_r,Median,LogT=TRUE,splitter=yearcat, StatType="AR")
Op.Median.UnTrans.SP<-BCa.Boot.CI(ScienceStatus_SP,O.Power,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Ep.Median.UnTrans.SP<-BCa.Boot.CI(ScienceStatus_SP,E.Power,Median,LogT=FALSE,splitter=yearcat, StatType="AR")

TIVA.results.SP<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,tiva.calc,LogT=FALSE,splitter=yearcat, StatType="AR")
Z.Median.UnTrans.SP<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Rindex.Reults.SP<-BCa.Boot.CI(ScienceStatus_SP,Calc.Z,r.index.calc.boot,LogT=FALSE,splitter=yearcat, StatType="AR")

PCurve.SP.Calculated<-BCa.Boot.CI(PerEvidence.ID,PerEvidence*100,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
PCurve.SP.Calculated

##manually add Pcurve Per study
Pcurve.S03<-data.frame(Analysis="pstudy",
                       Arithmetic=round(as.numeric(PerStudy.03.Pcurve$estimate)*100,2),
                       LB=round(PerStudy.03.Pcurve$conf.int[1]*100,2),
                       UP=round(PerStudy.03.Pcurve$conf.int[2]*100,2),
                       splitter="2003-2004")
Pcurve.S13<-data.frame(Analysis="pstudy",
                       Arithmetic=round(as.numeric(PerStudy.13.Pcurve$estimate)*100,2),
                       LB=round(PerStudy.13.Pcurve$conf.int[1]*100,2),
                       UP=round(PerStudy.13.Pcurve$conf.int[2]*100,2),
                       splitter="2013-2014")

Pcurve.perStudy<-rbind(Pcurve.S03,Pcurve.S13)

###Table
Table4.SP<-rbind(Pcurve.perStudy,PCurve.SP.Calculated,SS.Median.Trans.SP, Op.Median.UnTrans.SP,Ep.Median.UnTrans.SP,
                 TIVA.results.SP,Z.Median.UnTrans.SP,Rindex.Reults.SP)

kable(Table4.SP[Table4.SP$splitter=="2003-2004",], digits=2)
kable(Table4.SP[Table4.SP$splitter=="2013-2014",], digits=2)





###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
############################ Supplemental Analysis ################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


############################################################################
########################### Study Designs ##################################
############################################################################


# Of All studies
Design03.A<-summary(subset(ScienceStatus_SP,yearcat=="2003-2004", select=design))
Design13.A<-summary(subset(ScienceStatus_SP,yearcat=="2013-2014", select=design))

# Of studies that are experimental in nature
Design03.E<-summary(subset(ScienceStatus_SP,yearcat=="2003-2004" & expcor=="Experimental (i.e., all IVs were manipulated)" |
                 expcor=="Quasi-experimental (i.e., at least 1 IV is manipulated and at least 1 IV is measured)", select=design))
Design13.E<-summary(subset(ScienceStatus_SP,yearcat=="2013-2014" & expcor=="Experimental (i.e., all IVs were manipulated)" |
                 expcor=="Quasi-experimental (i.e., at least 1 IV is manipulated and at least 1 IV is measured)", select=design))

kable(cbind(Design03.A,Design13.A))
kable(cbind(Design03.E,Design13.E))

############################################################################
########################### Number of Citations ############################
############################################################################

##### Effect size & sample size

##Calculate n 
CountNM(ScienceStatus_SP$citations.number)/ length(ScienceStatus_SP$citations.number)

###Calcuate Mean
Cite.Mean.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,citations.number,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
Cite.Mean.UnTrans

## Correlations
Cite.by.SSr2<-ddply(ScienceStatus_SP, .(yearcat), summarize,
                    Corr.k.ss = cor(citations.number,sample.size_r,"kendall",use="complete.obs"),
                    Corr.k.r  = cor(citations.number,r2.effect,"kendall",use="complete.obs"),
                    Corr.k.Ep  = cor(citations.number,E.Power,"kendall",use="complete.obs"),
                    Corr.k.Op  = cor(citations.number,O.Power,"kendall",use="complete.obs"),
                    N1=CountNM(sample.size_r*sample.size_r),
                    N2=CountNM(sample.size_r*r2.effect),
                    N3=CountNM(sample.size_r*E.Power),
                    N4=CountNM(sample.size_r*O.Power))

Cite.by.SSr2

library(psych)
RtoZtest.SS.K<-paired.r(Cite.by.SSr2$Corr.k.ss[1],Cite.by.SSr2$Corr.k.ss[2],NULL, Cite.by.SSr2$N1[1],Cite.by.SSr2$N1[2])
RtoZtest.R.K <- paired.r(Cite.by.SSr2$Corr.k.r[1],Cite.by.SSr2$Corr.k.r[2],NULL, Cite.by.SSr2$N2[1],Cite.by.SSr2$N2[2])
RtoZtest.Ep.K <-paired.r(Cite.by.SSr2$Corr.k.Ep[1],Cite.by.SSr2$Corr.k.Ep[2],NULL, Cite.by.SSr2$N3[1],Cite.by.SSr2$N2[2])
RtoZtest.Op.K <-paired.r(Cite.by.SSr2$Corr.k.Op[1],Cite.by.SSr2$Corr.k.Op[2],NULL, Cite.by.SSr2$N4[1],Cite.by.SSr2$N2[2])
detach(psych)

rbind(
  c(RtoZtest.SS.K$z, RtoZtest.R.K$z, RtoZtest.Ep.K$z,RtoZtest.Op.K$z),
  c(RtoZtest.SS.K$p, RtoZtest.R.K$p, RtoZtest.Ep.K$p,RtoZtest.Op.K$p)
)

##############################################################
################# Percentile analysis ########################
##############################################################

Cut.Points.Cite<-c(0.45,.55,.90)

Mean.Cite<-ddply(ScienceStatus_SP, .(yearcat), summarize, 
                 MC=quantile(citations.number,Cut.Points.Cite,na.rm = TRUE)) 

Mean.Cite

######## Nested ifelse based on cuts for each year

ScienceStatus_SP$Cite.Cgroup<-ifelse(ScienceStatus_SP$yearcat=="2003-2004" & ScienceStatus_SP$citations.number >= Mean.Cite$MC[1] & ScienceStatus_SP$citations.number <= Mean.Cite$MC[2],
                                     1,
                                     ifelse(ScienceStatus_SP$yearcat=="2013-2014" & ScienceStatus_SP$citations.number >= Mean.Cite$MC[4] & ScienceStatus_SP$citations.number <= Mean.Cite$MC[5],
                                            1,
                                            ifelse(ScienceStatus_SP$yearcat=="2003-2004" & ScienceStatus_SP$citations.number > Mean.Cite$MC[3],
                                                   2,
                                                   ifelse(ScienceStatus_SP$yearcat=="2013-2014" & ScienceStatus_SP$citations.number > Mean.Cite$MC[6],
                                                          2,NA))))

ScienceStatus_SP$Cite.Cgroup<-factor(ScienceStatus_SP$Cite.Cgroup)


SS.Median.03.Trans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),sample.size_r,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
SS.Median.13.Trans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),sample.size_r,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

Aprior.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),Sample.N.test*100,Mean,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
Aprior.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),Sample.N.test*100,Mean,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

R.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),r2.effect,Mean,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
R.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),r2.effect,Mean,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

Op.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),O.Power,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
Op.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),O.Power,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

Ep.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),E.Power,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
Ep.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),E.Power,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

Z.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),Calc.Z,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
Z.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),Calc.Z,Median,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

Rindex.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),Calc.Z,r.index.calc.boot,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
Rindex.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),Calc.Z,r.index.calc.boot,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

TIVA.Mean.03.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2003-2004"),Calc.Z,tiva.calc,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")
TIVA.Mean.13.UnTrans<-BCa.Boot.CI(subset(ScienceStatus_SP,Cite.Cgroup!="NA" & yearcat=="2013-2014"),Calc.Z,tiva.calc,LogT=FALSE,splitter=Cite.Cgroup, StatType="AR")

CiteConnect<-ddply(ScienceStatus_SP, .(article.id), summarize,
                   Cite = max(Cite.Cgroup))

PerEvidence.ID.Cite<-merge(CiteConnect,PerEvidence.ID, by = "article.id")
PerEvidence.ID.Cite$Cite<-factor(PerEvidence.ID.Cite$Cite)

A.PCurve.SP.03.Cite<-BCa.Boot.CI(subset(PerEvidence.ID.Cite,Cite!="NA" & yearcat=="2003-2004"),PerEvidence*100,Mean,LogT=FALSE,splitter=Cite, StatType="AR")
A.PCurve.SP.13.Cite<-BCa.Boot.CI(subset(PerEvidence.ID.Cite,Cite!="NA" & yearcat=="2013-2014"),PerEvidence*100,Mean,LogT=FALSE,splitter=Cite, StatType="AR")

Calcuated.med.pvalue.Cite<-merge(CiteConnect,Calcuated.med.pvalue, by = "article.id")
Calcuated.med.pvalue.Cite$Cite<-factor(Calcuated.med.pvalue.Cite$Cite)

O.PCurve.SP.03.Cite<-BCa.Boot.CI(subset(Calcuated.med.pvalue.Cite,Cite!="NA" & yearcat=="2003-2004"),EV*100,Mean,LogT=FALSE,splitter=Cite, StatType="AR")
O.PCurve.SP.13.Cite<-BCa.Boot.CI(subset(Calcuated.med.pvalue.Cite,Cite!="NA" & yearcat=="2013-2014"),EV*100,Mean,LogT=FALSE,splitter=Cite, StatType="AR")


###Tables


Cite.03<-rbind(O.PCurve.SP.03.Cite,A.PCurve.SP.03.Cite,SS.Median.03.Trans,Aprior.Mean.03.UnTrans,
                Op.Mean.03.UnTrans,TIVA.Mean.03.UnTrans,Z.Mean.03.UnTrans,Rindex.Mean.03.UnTrans)

kable(Cite.03[Cite.03$splitter=="1",], digits=2)
kable(Cite.03[Cite.03$splitter=="2",], digits=2)

Cite.13<-rbind(O.PCurve.SP.13.Cite,A.PCurve.SP.13.Cite,SS.Median.13.Trans,Aprior.Mean.13.UnTrans,
               Op.Mean.13.UnTrans,TIVA.Mean.13.UnTrans,Z.Mean.13.UnTrans,Rindex.Mean.13.UnTrans)

kable(Cite.13[Cite.13$splitter=="1",], digits=2)
kable(Cite.13[Cite.13$splitter=="2",], digits=2)


#####################################################################################
#####Effect Size
#####################################################################################

##Calculate n 
CountNM(ScienceStatus_SP$r2.effect)/ length(ScienceStatus_SP$article.id)


###Calcuate Mean, Median, Peak
R2.Mean.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,r2.effect,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
R2.Median.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,r2.effect,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
R2.Mode.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,r2.effect,Peak,LogT=FALSE,splitter=yearcat, StatType="AR")

r2toD(R2.Median.UnTrans[1,2])
r2toD(R2.Median.UnTrans[2,2])

##Change in sample need: 
Tpoint1 <-
  pwr.t.test(n = NULL, d = r2toD(R2.Median.UnTrans[1,2]), sig.level = 0.05, power = .8,
             type = c("two.sample"),
             alternative = "two.sided")

Tpoint2 <-
  pwr.t.test(n = NULL, d = r2toD(R2.Median.UnTrans[2,2]), sig.level = 0.05, power = .8,
             type = c("two.sample"),
             alternative = "two.sided")
Tpoint1$n
Tpoint2$n
(Tpoint2$n-Tpoint1$n)/(Tpoint1$n)*100



#########
#### Table 
########
R2.Forest.UnTrans<-rbind(R2.Mean.UnTrans,R2.Median.UnTrans,R2.Mode.UnTrans)

kable(R2.Forest.UnTrans, digits=2)

#########
#### Forest Plots 
#########
R2.Forest.Plot<-ggplot(R2.Forest.UnTrans, aes(Arithmetic, Analysis, colour = splitter, shape = splitter))+
  #facet_wrap( ~ PType, ncol=2)+
  geom_point(size=3) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .2),size = 1,alpha=.65)+
  ylab("")+
  xlab(expression(Estimated~Effect~Size~R^{2}))+
  scale_x_continuous(breaks = seq(0, .25, .05),limits=c(0,.25))+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
R2.Forest.Plot
#ggsave(R2.Forest.Plot,file="R2.Forest.Plot.png",scale = 1, width=3, height=3,units = c("in"))


#########
#####Density Plot 
#########

R2.plot<-ggplot(data = ScienceStatus_SP, aes(r2.effect,..scaled..,fill=yearcat,color=yearcat))+
  facet_wrap( ~ yearcat, ncol=2)+
  coord_cartesian(ylim = c(0,1))+
  stat_density(geom = "area", position = "stack",adjust=.5, bw = 'Sj',kernel='gaussian',trim = FALSE,alpha=.65,na.rm=TRUE)+
  ylab("Density")+
  xlab(expression(Estimated~Effect~Size~R^{2}))+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none") 
R2.plot

R2.grid<-grid.arrange(R2.plot,R2.Forest.Plot, nrow=2, ncol=1,heights = c(3, 4))
ggsave(R2.grid,file="FigureS2.tiff",scale = 1, dpi=300,width=4, height=4,units = c("in"))


### Test denesity plots differences
x <- subset(ScienceStatus_SP, yearcat =="2003-2004" & r2.effect!="NA", select =r2.effect  )
y <- subset(ScienceStatus_SP, yearcat =="2013-2014" & r2.effect!="NA", select =r2.effect  )
R2.D.Entropy<-npunitest(as.vector(x$r2.effect),as.vector(y$r2.effect),boot.num=200)
R2.D.Entropy



#####################################################################################
#####Post-Hoc Power
#####################################################################################
cor(ScienceStatus_SP$E.Power,ScienceStatus_SP$O.Power,use = "pairwise.complete.obs",method = c("spearman"))

CountNM(ScienceStatus_SP$E.Power*ScienceStatus_SP$O.Power)
##Calculate n 
CountNM(ScienceStatus_SP$O.Power)/ length(ScienceStatus_SP$article.id)

###Calcuate Mean, Median, Peak
Ep.Mean.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,E.Power,Mean,LogT=FALSE,splitter=yearcat, StatType="AR")
Ep.Median.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,E.Power,Median,LogT=FALSE,splitter=yearcat, StatType="AR")
Ep.Mode.UnTrans<-BCa.Boot.CI(ScienceStatus_SP,E.Power,Peak,LogT=FALSE,splitter=yearcat, StatType="AR")


#########
#### Table for Post-Hoc Power
########
EP.Forest.Trans<-rbind(Ep.Mean.UnTrans,Ep.Median.UnTrans,Ep.Mode.UnTrans)
kable(EP.Forest.Trans, digits=2)

#########
#### Forest Plots on Post-Hoc Power
########

Ep.Forest.Plot<-ggplot(EP.Forest.Trans, aes(Arithmetic, Analysis, colour = splitter, shape = splitter))+
  geom_point(size=3) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .2),size = 1,alpha=.65)+
  ylab("")+xlab('Post-Hoc Power')+
  scale_x_continuous(breaks = seq(.7, 1, .05),limits=c(.7,1))+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  theme_bw()+
  theme(text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
Ep.Forest.Plot
#ggsave(Ep.Forest.Plot,file="Ep.Forest.Plot.png",scale = 1, width=3, height=3,units = c("in"))


#########
#####Density Plot on Post-Hoc Power
#########

Power.plot.PostHoc<-ggplot(data = ScienceStatus_SP, aes(E.Power,..scaled..,fill=yearcat,color=yearcat))+
  facet_wrap( ~ yearcat, ncol=2)+
  coord_cartesian(ylim = c(0,1))+
  stat_density(geom = "area", position = "stack",adjust=.5, bw = 'Sj',kernel='gaussian',trim = FALSE,alpha=.65,na.rm=TRUE)+
  ylab("Density")+xlab('Post-Hoc Power')+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none") 
Power.plot.PostHoc

PH.grid<-grid.arrange(Power.plot.PostHoc,Ep.Forest.Plot, nrow=2, ncol=1,heights = c(3, 4))
ggsave(PH.grid,file="FigureS1.tiff",dpi=300,scale = 1, width=4, height=4,units = c("in"))


### Test denesity plots differences
x <- subset(ScienceStatus_SP, yearcat =="2003-2004" & E.Power!="NA", select =E.Power  )
y <- subset(ScienceStatus_SP, yearcat =="2013-2014" & E.Power!="NA", select =E.Power  )
Ph.D.Entropy<-npunitest(as.vector(x$E.Power),as.vector(y$E.Power),boot.num=200)
Ph.D.Entropy



############################################################################
#################### A prior power - Per Journal ###########################
############################################################################

SS.Median.Trans.0304<-BCa.Boot.CI(subset(ScienceStatus_SP, yearcat=="2003-2004"),sample.size_r,Median,LogT=TRUE,splitter=journal, StatType="AR")
SS.Median.Trans.1314<-BCa.Boot.CI(subset(ScienceStatus_SP, yearcat=="2013-2014"),sample.size_r,Median,LogT=TRUE,splitter=journal, StatType="AR")

Pj1<-matrix(0,4,6)
nsize1<-SS.Median.Trans.0304$Arithmetic
nsize2<-SS.Median.Trans.0304$LB
nsize3<-SS.Median.Trans.0304$UP
nsize4<-SS.Median.Trans.1314$Arithmetic
nsize5<-SS.Median.Trans.1314$LB
nsize6<-SS.Median.Trans.1314$UP

for (i in 1:4) {
  
  Pj1[i,1]<-pwr.t.test(n = nsize1[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  Pj1[i,2]<-pwr.t.test(n = nsize2[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  Pj1[i,3]<-pwr.t.test(n = nsize3[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  
  Pj1[i,4]<-pwr.t.test(n = nsize4[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  Pj1[i,5]<-pwr.t.test(n = nsize5[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  Pj1[i,6]<-pwr.t.test(n = nsize6[i]/2, d = .430, sig.level = 0.05, power = NULL,
                       type = c("two.sample"),
                       alternative = "two.sided")$power
  
}

FV.Plot.1<-rbind(cbind(SS.Median.Trans.0304,data.frame(Pj1[1:4,1:3])),
                 cbind(SS.Median.Trans.1314,data.frame(Pj1[1:4,4:6])))
FV.Plot.1$yearcat<-c(rep("2003-2004",4),rep("2013-2014",4))


### Calcuate WT median based on number of journals
library(matrixStats)
WTmedian03<-weightedMedian(FV.Plot.1$X1[1:4],PerStudyAbove.J$n[1:4])
WTmedian13<-weightedMedian(FV.Plot.1$X1[5:8],PerStudyAbove.J$n[5:8])
hline.data <- data.frame(z = c(WTmedian03,WTmedian13), yearcat = c("2003-2004","2013-2014")) 
detach("package:matrixStats", unload=TRUE)


pd <- position_dodge(0.1)

FigureS3<-ggplot(data = FV.Plot.1, aes(x=reorder(splitter,X1),y=(X1)))+
  facet_wrap( ~ yearcat,nrow=2)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=X2, ymax=X3), colour="black", width=.2, position=position_dodge(.9)) +
  xlab("")+
  scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1))+
  ylab("Statistical power to detect a d = .430")+
  geom_hline(aes(yintercept = z), hline.data)+
  coord_flip(ylim = c(0,1))+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","black"))+
  theme_bw()+
  theme(text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="none") 
FigureS3

ggsave(FigureS3,file="FigureS3.tiff",scale = 1, dpi = 300,width=4, height=4,units = c("in"))








