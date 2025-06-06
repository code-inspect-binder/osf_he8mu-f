library(foreign)
library(ggplot2)
library(boot)
library(knitr)
library(reshape2)
library(plyr)

setwd("C:/AlexFiles/SugerSync/UIC/ScienceProject/ForPaper1")


Study1<-read.spss("Perceptions_of_Replicability_deidentified.sav",use.value.labels = FALSE, to.data.frame = TRUE,
                  max.value.labels = Inf, trim.factor.names = FALSE,
                  trim_values = TRUE, reencode = NA, use.missings = to.data.frame)


####################################################################
################Figure 1. Perceived Replicability
####################################################################

#######now
JPSP1<-BCa.Boot.CI(Study1,estrepnow_1,Mean,LogT=FALSE, StatType="AR")
PS1<-BCa.Boot.CI(Study1,estrepnow_2,Mean,LogT=FALSE, StatType="AR")
JESP1<-BCa.Boot.CI(Study1,estrepnow_3,Mean,LogT=FALSE, StatType="AR")
PSPB1<-BCa.Boot.CI(Study1,estrepnow_4,Mean,LogT=FALSE, StatType="AR")

#bind
PowerNow<-rbind(JPSP1,PS1,JESP1,PSPB1)
PowerNow$time<-"Now"
PowerNow$Analysis<-(c("JPSP","PS","JESP","PSPB"))

######10 years ago
JPSP2<-BCa.Boot.CI(Study1,estreppast_1,Mean,LogT=FALSE, StatType="AR")
PS2<-BCa.Boot.CI(Study1,estreppast_2,Mean,LogT=FALSE, StatType="AR")
JESP2<-BCa.Boot.CI(Study1,estreppast_3,Mean,LogT=FALSE, StatType="AR")
PSPB2<-BCa.Boot.CI(Study1,estreppast_4,Mean,LogT=FALSE, StatType="AR")

#bind
PowerPast<-rbind(JPSP2,PS2,JESP2,PSPB2)
PowerPast$time<-"10 yrs ago"
PowerPast$Analysis<-(c("JPSP","PS","JESP","PSPB"))

#bind now and 10 years ago
Perceived.Power<-rbind(PowerNow,PowerPast)

###Tables of values 
kable(Perceived.Power, digits=2)

##Reorder
c("JPSP","PSPB","PS","JESP")
c("JESP","PS","PSPB","JPSP")
Perceived.Power$Analysis2 <- factor(Perceived.Power$Analysis, levels=c("PS","JESP","PSPB","JPSP"))

##extract attri for graph
Labels1<-attributes(Study1$estrepnow_1)$value.labels

Labels1<-c("0-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%")

####forest plot
ForestPlot.Replicability<-ggplot(Perceived.Power, aes(Arithmetic, Analysis2, colour = time, shape = time))+
geom_point(size=1.5) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .3),size = .75,alpha=.65)+
  ylab("")+xlab('Perceived Replicability')+
  scale_fill_manual(values=c("black","grey40"))+
  scale_color_manual(values=c("black","grey40"))+
  #scale_color_manual(values=c("black","grey40"))+
  scale_x_continuous(limits =c(1, 10), breaks = c(1,2,3,4,5,6,7,8,9,10),labels = Labels1)+
  #scale_linetype_manual(values = c("dashed", "solid"))+ 
  theme_bw()+
  #labs(colour = "Cylinders")+
  theme(text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(color = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.title= element_blank(),
        legend.margin = unit(.001, "cm"),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")
ForestPlot.Replicability

ggsave(ForestPlot.Replicability,file="Figure1.tiff",scale = 1, width=5, height=2,units = c("in"))

#######################################################################
#################### Figure 2
#######################################################################

pow.freq<-BCa.Boot.CI(Study1,pow.freq,Mean,LogT=FALSE, StatType="AR")
fal.freq<-BCa.Boot.CI(Study1,fal.freq,Mean,LogT=FALSE, StatType="AR")
irr.freq<-BCa.Boot.CI(Study1,irr.freq,Mean,LogT=FALSE, StatType="AR")
exp.freq<-BCa.Boot.CI(Study1,exp.freq,Mean,LogT=FALSE, StatType="AR")
sel.freq<-BCa.Boot.CI(Study1,sel.freq,Mean,LogT=FALSE, StatType="AR")
sto.freq<-BCa.Boot.CI(Study1,sto.freq,Mean,LogT=FALSE, StatType="AR")
rou.freq<-BCa.Boot.CI(Study1,rou.freq,Mean,LogT=FALSE, StatType="AR")
exc.freq<-BCa.Boot.CI(Study1,exc.freq,Mean,LogT=FALSE, StatType="AR")
dvs.freq<-BCa.Boot.CI(Study1,dvs.freq,Mean,LogT=FALSE, StatType="AR")
eff.freq<-BCa.Boot.CI(Study1,eff.freq,Mean,LogT=FALSE, StatType="AR")
dat.freq<-BCa.Boot.CI(Study1,dat.freq,Mean,LogT=FALSE, StatType="AR")
reg.freq<-BCa.Boot.CI(Study1,reg.freq,Mean,LogT=FALSE, StatType="AR")
add.freq<-BCa.Boot.CI(Study1,add.freq,Mean,LogT=FALSE, StatType="AR")
pee.freq<-BCa.Boot.CI(Study1,pee.freq,Mean,LogT=FALSE, StatType="AR")

Frequency.Data<-rbind(pow.freq,fal.freq,irr.freq,exp.freq,sel.freq,sto.freq,rou.freq,
                 exc.freq,dvs.freq,eff.freq,dat.freq,reg.freq,add.freq,pee.freq)
          
Frequency.Data$Analysis<-c("Conduct power analyses","Falsify data","Falsely claim results unaffected by demographics","Report an unexpected finding were predicted",
"Selectively report studies that worked","Stop data collection early","Round off a p-value that is just over .05",
"Not report all conditions","Collect data on multiple DVs","Report effect sizes",
"Make data publicly available","Pre-register hypotheses","Decide to collect additional data after looking",
"Exclude data after looking")
Frequency.Data$Analysis


###Custom write the order
Frequency.Data$ReOrderer<-c(13,1,2,8,10,3,4,6,9,14,12,11,7,5)
Frequency.Data$Colorer<-c(1,0,0,0,0,0,0,0,0,1,1,1,0,0)
Labels2<-c("Never","Rarely","Sometimes","Often","Always")

research.practices<-ggplot(Frequency.Data, aes(Arithmetic, reorder(Analysis,-ReOrderer)))+
  geom_point(size=1.25) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .3),size = .75,alpha=.65)+
  ylab("")+xlab('')+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  scale_x_continuous(limits =c(1, 5), breaks = c(1,2,3,4,5),labels = Labels2)+
  #scale_linetype_manual(values = c("dashed", "solid"))+ 
  theme_bw()+
  #labs(colour = "Cylinders")+
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
research.practices



ggsave(research.practices,file="Figure2.tiff",scale = 1, width=6, height=2.75,units = c("in"))





#######################################################################
#################### Figure 3
#######################################################################

pow.acc<-BCa.Boot.CI(Study1,pow.acc,Mean,LogT=FALSE, StatType="AR")
fal.acc<-BCa.Boot.CI(Study1,fal.acc,Mean,LogT=FALSE, StatType="AR")
irr.acc<-BCa.Boot.CI(Study1,irr.acc,Mean,LogT=FALSE, StatType="AR")
exp.acc<-BCa.Boot.CI(Study1,exp.acc,Mean,LogT=FALSE, StatType="AR")
sel.acc<-BCa.Boot.CI(Study1,sel.acc,Mean,LogT=FALSE, StatType="AR")
sto.acc<-BCa.Boot.CI(Study1,sto.acc,Mean,LogT=FALSE, StatType="AR")
rou.acc<-BCa.Boot.CI(Study1,rou.acc,Mean,LogT=FALSE, StatType="AR")
exc.acc<-BCa.Boot.CI(Study1,exc.acc,Mean,LogT=FALSE, StatType="AR")
dvs.acc<-BCa.Boot.CI(Study1,dvs.acc,Mean,LogT=FALSE, StatType="AR")
eff.acc<-BCa.Boot.CI(Study1,eff.acc,Mean,LogT=FALSE, StatType="AR")
dat.acc<-BCa.Boot.CI(Study1,dat.acc,Mean,LogT=FALSE, StatType="AR")
reg.acc<-BCa.Boot.CI(Study1,reg.acc,Mean,LogT=FALSE, StatType="AR")
add.acc<-BCa.Boot.CI(Study1,add.acc,Mean,LogT=FALSE, StatType="AR")
pee.acc<-BCa.Boot.CI(Study1,pee.acc,Mean,LogT=FALSE, StatType="AR")

Frequency.acc<-rbind(pow.acc,fal.acc,irr.acc,exp.acc,sel.acc,sto.acc,rou.acc,
                      exc.acc,dvs.acc,eff.acc,dat.acc,reg.acc,add.acc,pee.acc)


Frequency.acc$Analysis<-c("Conduct power analyses","Falsify data","Falsely claim results unaffected by demographics","Report an unexpected finding were predicted",
                           "Selectively report studies that worked","Stop data collection early","Round off a p-value that is just over .05",
                           "Not report all conditions","Collect data on multiple DVs","Report effect sizes",
                           "Make data publicly available","Pre-register hypotheses","Decide to collect additional data after looking",
                           "Exclude data after looking")


Labels3<-c("Very\nunacceptable","Moderately\nunacceptable", "Slightly\nunacceptable",           
          "Uncertain", "Slightly\nacceptable","Moderately\nacceptable","Very\nacceptable")


Frequency.acc$ReOrderer<-c(13,1,2,8,10,3,4,6,9,14,12,11,7,5)
Frequency.acc$Colorer<-c(1,0,0,0,0,0,0,0,0,1,1,1,0,0)


research.acc<-ggplot(Frequency.acc, aes(Arithmetic, reorder(Analysis,-ReOrderer)),color=as.factor(Colorer))+
  geom_point(size=1.25,aes(Color=Colorer)) +
  geom_errorbarh(aes(xmax = UP, xmin = LB,height = .3),size = .75,alpha=.65)+
  ylab("")+xlab('')+
  scale_fill_manual(values=c("#ED1C24","#00AEEF"))+
  scale_color_manual(values=c("black","grey40"))+
  scale_x_continuous(limits =c(-3, 3), breaks = c(-3,-2,-1,0,1,2,3),labels = Labels3)+
  #scale_linetype_manual(values = c("dashed", "solid"))+ 
  geom_vline(xintercept = 0)+
  theme_bw()+
  #labs(colour = "Cylinders")+
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
research.acc


ggsave(research.acc,file="Figure3.tiff",scale = 1, width=7.0, height=2.75,units = c("in"))





#######################################################################
#################### Figure 4
#######################################################################

###I should loop this, but being lazy.....

  
Percent1 <-function(x,i) {
  Num<-sum(x[i]==1,na.rm = TRUE)
  Dem<-sum( !is.na( x[i] )) 
  Per<-(Num/Dem)*100}

pow.change1<-BCa.Boot.CI(Study1,pow.change,Percent1,LogT=FALSE, StatType="AR")
fal.change1<-BCa.Boot.CI(Study1,fal.change,Percent1,LogT=FALSE, StatType="AR")
irr.change1<-BCa.Boot.CI(Study1,irr.change,Percent1,LogT=FALSE, StatType="AR")
exp.change1<-BCa.Boot.CI(Study1,exp.change,Percent1,LogT=FALSE, StatType="AR")
sel.change1<-BCa.Boot.CI(Study1,sel.change,Percent1,LogT=FALSE, StatType="AR")
sto.change1<-BCa.Boot.CI(Study1,sto.change,Percent1,LogT=FALSE, StatType="AR")
rou.change1<-BCa.Boot.CI(Study1,rou.change,Percent1,LogT=FALSE, StatType="AR")
exc.change1<-BCa.Boot.CI(Study1,exc.change,Percent1,LogT=FALSE, StatType="AR")
dvs.change1<-BCa.Boot.CI(Study1,dvs.change,Percent1,LogT=FALSE, StatType="AR")
eff.change1<-BCa.Boot.CI(Study1,eff.change,Percent1,LogT=FALSE, StatType="AR")
dat.change1<-BCa.Boot.CI(Study1,dat.change,Percent1,LogT=FALSE, StatType="AR")
reg.change1<-BCa.Boot.CI(Study1,reg.change,Percent1,LogT=FALSE, StatType="AR")
add.change1<-BCa.Boot.CI(Study1,add.change,Percent1,LogT=FALSE, StatType="AR")
pee.change1<-BCa.Boot.CI(Study1,pee.change,Percent1,LogT=FALSE, StatType="AR")

P.1<-rbind(pow.change1,fal.change1,irr.change1,exp.change1,sel.change1,sto.change1,rou.change1,
                      exc.change1,dvs.change1,eff.change1,dat.change1,reg.change1,add.change1,pee.change1)

P.1$Rating<-"Decreased"
P.1$Analysis<-c("Conduct power analyses","Falsify data","Falsely claim results unaffected by demographics","Report an unexpected finding were predicted",
                  "Selectively report studies that worked","Stop data collection early","Round off a p-value that is just over .05",
                  "Not report all conditions","Collect data on multiple DVs","Report effect sizes",
                  "Make data publicly available","Pre-register hypotheses","Decide to collect additional data after looking",
                  "Exclude data after looking")

###############

Percent2 <-function(x,i) {
  Num<-sum(x[i]==2,na.rm = TRUE)
  Dem<-sum( !is.na( x[i] )) 
  Per<-(Num/Dem)*100}

pow.change2<-BCa.Boot.CI(Study1,pow.change,Percent2,LogT=FALSE, StatType="AR")
fal.change2<-BCa.Boot.CI(Study1,fal.change,Percent2,LogT=FALSE, StatType="AR")
irr.change2<-BCa.Boot.CI(Study1,irr.change,Percent2,LogT=FALSE, StatType="AR")
exp.change2<-BCa.Boot.CI(Study1,exp.change,Percent2,LogT=FALSE, StatType="AR")
sel.change2<-BCa.Boot.CI(Study1,sel.change,Percent2,LogT=FALSE, StatType="AR")
sto.change2<-BCa.Boot.CI(Study1,sto.change,Percent2,LogT=FALSE, StatType="AR")
rou.change2<-BCa.Boot.CI(Study1,rou.change,Percent2,LogT=FALSE, StatType="AR")
exc.change2<-BCa.Boot.CI(Study1,exc.change,Percent2,LogT=FALSE, StatType="AR")
dvs.change2<-BCa.Boot.CI(Study1,dvs.change,Percent2,LogT=FALSE, StatType="AR")
eff.change2<-BCa.Boot.CI(Study1,eff.change,Percent2,LogT=FALSE, StatType="AR")
dat.change2<-BCa.Boot.CI(Study1,dat.change,Percent2,LogT=FALSE, StatType="AR")
reg.change2<-BCa.Boot.CI(Study1,reg.change,Percent2,LogT=FALSE, StatType="AR")
add.change2<-BCa.Boot.CI(Study1,add.change,Percent2,LogT=FALSE, StatType="AR")
pee.change2<-BCa.Boot.CI(Study1,pee.change,Percent2,LogT=FALSE, StatType="AR")

P.2<-rbind(pow.change2,fal.change2,irr.change2,exp.change2,sel.change2,sto.change2,rou.change2,
           exc.change2,dvs.change2,eff.change2,dat.change2,reg.change2,add.change2,pee.change2)

P.2$Rating<-"No Change"
P.2$Analysis<-c("Conduct power analyses","Falsify data","Falsely claim results unaffected by demographics","Report an unexpected finding were predicted",
                "Selectively report studies that worked","Stop data collection early","Round off a p-value that is just over .05",
                "Not report all conditions","Collect data on multiple DVs","Report effect sizes",
                "Make data publicly available","Pre-register hypotheses","Decide to collect additional data after looking",
                "Exclude data after looking")

#################

Percent3 <-function(x,i) {
  Num<-sum(x[i]==3,na.rm = TRUE)
  Dem<-sum( !is.na( x[i] )) 
  Per<-(Num/Dem)*100}

pow.change3<-BCa.Boot.CI(Study1,pow.change,Percent3,LogT=FALSE, StatType="AR")
fal.change3<-BCa.Boot.CI(Study1,fal.change,Percent3,LogT=FALSE, StatType="AR")
irr.change3<-BCa.Boot.CI(Study1,irr.change,Percent3,LogT=FALSE, StatType="AR")
exp.change3<-BCa.Boot.CI(Study1,exp.change,Percent3,LogT=FALSE, StatType="AR")
sel.change3<-BCa.Boot.CI(Study1,sel.change,Percent3,LogT=FALSE, StatType="AR")
sto.change3<-BCa.Boot.CI(Study1,sto.change,Percent3,LogT=FALSE, StatType="AR")
rou.change3<-BCa.Boot.CI(Study1,rou.change,Percent3,LogT=FALSE, StatType="AR")
exc.change3<-BCa.Boot.CI(Study1,exc.change,Percent3,LogT=FALSE, StatType="AR")
dvs.change3<-BCa.Boot.CI(Study1,dvs.change,Percent3,LogT=FALSE, StatType="AR")
eff.change3<-BCa.Boot.CI(Study1,eff.change,Percent3,LogT=FALSE, StatType="AR")
dat.change3<-BCa.Boot.CI(Study1,dat.change,Percent3,LogT=FALSE, StatType="AR")
reg.change3<-BCa.Boot.CI(Study1,reg.change,Percent3,LogT=FALSE, StatType="AR")
add.change3<-BCa.Boot.CI(Study1,add.change,Percent3,LogT=FALSE, StatType="AR")
pee.change3<-BCa.Boot.CI(Study1,pee.change,Percent3,LogT=FALSE, StatType="AR")

P.3<-rbind(pow.change3,fal.change3,irr.change3,exp.change3,sel.change3,sto.change3,rou.change3,
           exc.change3,dvs.change3,eff.change3,dat.change3,reg.change3,add.change3,pee.change3)

P.3$Rating<-"Increased"
P.3$Analysis<-c("Conduct power analyses","Falsify data","Falsely claim results unaffected by demographics","Report an unexpected finding were predicted",
                "Selectively report studies that worked","Stop data collection early","Round off a p-value that is just over .05",
                "Not report all conditions","Collect data on multiple DVs","Report effect sizes",
                "Make data publicly available","Pre-register hypotheses","Decide to collect additional data after looking",
                "Exclude data after looking")


Research.change<-rbind(P.1,P.2,P.3)

####Stacked plot

#Manually reorder
ReOrderer<-c(13,1,2,8,10,3,4,6,9,14,12,11,7,5)

### CHANGE d3 to be this and run from d3 down the plot!

ReOrderer<-c(13,1,2,8,10,3,4,6,9,14,12,11,7,5)

d3 <- P.1[order(-ReOrderer,P.1$Analysis),]
P.1a <- arrange(transform(P.1,Analysis=factor(Analysis,levels=d3$Analysis)),Analysis)
P.2a <- arrange(transform(P.2,Analysis=factor(Analysis,levels=d3$Analysis)),Analysis)
P.3a <- arrange(transform(P.3,Analysis=factor(Analysis,levels=d3$Analysis)),Analysis)


Research.change<-rbind(P.1a,P.2a,P.3a)
Research.change$Rating <- factor(Research.change$Rating, levels = c("Decreased", "No Change","Increased")) 

####Stacked plot
Research.change.Plot<-ggplot(Research.change, aes(x=Analysis, y=Arithmetic, fill = Rating))+
  geom_bar(stat='identity',position = "stack")+coord_flip()+
  xlab("")+ylab('Percentage')+
  scale_fill_manual(values=c("black","gray80","gray40"))+
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
Research.change.Plot

ggsave(Research.change.Plot,file="Figure4.tiff",scale = 1, width=7.0, height=2.75,units = c("in"))


