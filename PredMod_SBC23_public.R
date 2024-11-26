
#### Ethnoarchaeological inductive predictive Model ####
#### Sorgenti del Brembo di Carona - v.2023 ####
#### E. Croce, F, Carrer ####

#### DATA UPLOAD ####
library(rgrass)
library(terra)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("D:/Archeologia/Carona/Analisi/Modelli/Modello_2023")

# Spatial data upload
Carona<-read_RAST(c("DTM_SBC","Slope","TWI","River_cost","Lakes_cost",
                  "North","East","Profile","Lithology","Permeability",
                  "Avalanche","MorFeat"))
Carona$Lithology<-as.factor(Carona$Lithology)
Carona$Permeability<-as.factor(Carona$Permeability)
Carona$Avalanche<-as.factor(Carona$Avalanche)
Carona$MorFeat<-as.factor(Carona$MorFeat)
Sites<-read_VECT("Summer_Farms_Points")
No_Sites<-read_VECT("Random90")

# Tables creation
sitestab<-data.frame(extract(Carona, Sites),name=rep(1,length(Sites)))
nositestab<-data.frame(extract(Carona, No_Sites),name=rep(0,length(No_Sites)))
tab<-rbind(sitestab,nositestab)
tab$Lithology<-as.factor(tab$Lithology)
tab$Permeability<-as.factor(tab$Permeability)
tab$Avalanche<-as.factor(tab$Avalanche)
tab$MorFeat<-as.factor(tab$MorFeat)
tab$evidence<-as.factor(tab$name)
tab_cat<-data.frame(tab$Lithology,tab$Permeability,tab$Avalanche,tab$MorFeat,tab$name)
colnames(tab_cat)<-c("Lithology","Permeability","Avalanche","MorFeat","name")
tab_num<-data.frame(tab$DTM_SBC,tab$Slope,tab$TWI,tab$River_cost,tab$Lakes_cost,tab$North,tab$East,tab$Profile,tab$name)
colnames(tab_num)<-c("DTM_SBC","Slope","TWI","River_cost","Lakes_cost","North","East","Profile","name")
sites_num<-subset(tab_num, name == 1)
sites_num<-select(sites_num,-name)
nosites_num<-subset(tab_num, name == 0)
nosites_num<-select(nosites_num,-name)
sites_cat<-subset(tab_cat, name == 1)
sites_cat<-select(sites_cat, -name)
nosites_cat<-subset(tab_cat, name == 0)
nosites_cat<-select(nosites_cat, -name)

#### COLLINEARITY ####
colround<-round(cor(select(tab_num, -name),method="pearson"),4)
as.data.frame(colround) %>% filter_all(any_vars(.>0.7 & .<1))

#### KOLMOGOROV-SMIRNOV TEST ####

library(Matching)

# Function to analyze all the record at once
# x = sites
# y = no-sites
# nsim = number of repetitions

ksb.col<-function(x,y,nsim){
  ksb_tab<-list()
  for(i in 1:ncol(x)){
    ksb_tab[[i]]<-ks.boot(x[,i],y[,i],nboots=nsim)
  }
  names(ksb_tab)<-colnames(x)
  return(lapply(ksb_tab,summary))
}

ks_file<-capture.output(ksb.col(sites_num,nosites_num,10000))

#### CHI-SQUARE TEST - SITES####

# Lithology
Lithology_sites <- data.frame ('Lithology' = sitestab$Lithology) %>%
  mutate(Lithology = factor(Lithology, levels = c('1','2','3','4'))) %>%
  count (Lithology, .drop = FALSE) 

Lithology_area <- data.frame (Carona$Lithology) %>%
  group_by (Lithology) %>%
  summarize (pixels = n()) %>%
  mutate (perc = pixels/sum(pixels))%>%
  mutate (sites = Lithology_sites$n)

chisq.test(x=Lithology_area$sites,p=Lithology_area$perc,simulate.p.value =  TRUE)

# Permeability
Permeability_sites <- data.frame ('Permeability' = sitestab$Permeability) %>%
  mutate(Permeability = factor(Permeability, levels = c('1','2','3'))) %>%
  count (Permeability, .drop = FALSE) 

Permeability_area <- data.frame (Carona$Permeability) %>%
  group_by (Permeability) %>%
  summarize (pixels = n()) %>%
  mutate (perc = pixels/sum(pixels))%>%
  mutate (sites = Permeability_sites$n)

chisq.test(x=Permeability_area$sites,p=Permeability_area$perc,simulate.p.value =  TRUE)

# Avalanche
Avalanche_sites <- data.frame ('Avalanche' = sitestab$Avalanche) %>%
  mutate(Avalanche = factor(Avalanche, levels = c('1','2','3'))) %>%
  count (Avalanche, .drop = FALSE) 

Avalanche_area <- data.frame (Carona$Avalanche) %>%
  group_by (Avalanche) %>%
  summarize (pixels = n()) %>%
  mutate (perc = pixels/sum(pixels))%>%
  mutate (sites = Avalanche_sites$n)

chisq.test(x=Avalanche_area$sites,p=Avalanche_area$perc,simulate.p.value =  TRUE)

# MorFeat
MorFeat_sites <- data.frame ('MorFeat' = sitestab$MorFeat) %>%
  mutate(MorFeat = factor(MorFeat, levels = c('1','2','3','4','5','6'))) %>%
  count (MorFeat, .drop = FALSE) 

MorFeat_area <- data.frame (Carona$MorFeat) %>%
  group_by (MorFeat) %>%
  summarize (pixels = n()) %>%
  mutate (perc = pixels/sum(pixels))%>%
  mutate (sites = MorFeat_sites$n)

chisq.test(x=MorFeat_area$sites,p=MorFeat_area$perc,simulate.p.value =  TRUE)

#### UNIVARIATE LOGISTIC REGRESSION ####

univar<-function(x,y){
  lrunivar<-list()
  for(i in x){
    f<-formula(paste("name","~",i))
    lrunivar[[i]]<-glm(f,data=y,family=binomial(logit))
  }
  return(lrunivar)
}

tab$Lithology<-relevel(tab$Lithology, ref='4')
tab$Avalanche<-relevel(tab$Avalanche, ref='3')

indep<-c("DTM_SBC","Slope","TWI","River_cost","Lithology","Avalanche")

glm_uni<-univar(indep,tab)

# Likelihood ratio test

ratio<-capture.output(1-(pchisq((sapply(glm_uni,function(x){x$null.deviance}))-
                                  (sapply(glm_uni,function(x){x$deviance})),1)))

#### MULTIVARIATE LOGISTIC REGRESSION ####

library(MASS)

glm_mul<-glm(name~DTM_SBC+Slope+TWI+River_cost+Lithology+Avalanche,
             data=tab,family=binomial(logit))

glmfile<-capture.output(summary(glm_mul))

# Variance inflation factor

library(car)

vif_MUL<-vif(glm_mul)
vifile<-capture.output(vif_MUL)

# Akaike Information Criterion

AIC_glm<-stepAIC(glm_mul)
summary(AIC_glm)
aicfile<-capture.output(summary(AIC_glm))

# Bayesian Information Criterion

BIC_glm<-stepAIC(glm_mul,k=log(length(tab$name)))
bicfile<-capture.output(summary(BIC_glm))

# coefficients standardization (Vaughn & Crawford 2009:550)

library(fmsb)

bx<-summary(BIC_glm)$coef[-1,1]

R<-NagelkerkeR2(BIC_glm)
R<-R$R2

dtm<-sd(as.numeric(unlist(BIC_glm$model[2])))

slope<-sd(as.numeric(unlist(BIC_glm$model[3])))

sx<-cbind(dtm,slope)

sy<-sd(as.numeric(unlist(BIC_glm$model[1])))

(bx*sx*R)/sy

# AUC - ROC

library(pROC)

roc<-roc(BIC_glm$y,BIC_glm$fitted.values)

jpeg("AUC.jpeg",width = 2000, height = 2000, units = "px",
     res = 300)
plot.roc(roc,reuse.auc = T,legacy.axes = T,main="Area under the ROC curve")
dev.off()

write("\nArea under the ROC curve (AUC)",file="Multivariate_Regression.txt",append=T)
write(capture.output(roc),file="Multivariate_Regression.txt",append=T)
write("\nDiscriminatory Ability (HOSMER et. al. 2013 - Applied Logistic Regression - 3rd Ed.) \n",
      file="Multivariate_Regression.txt",append=T)
write("    0.5 = No better than chance", file="Multivariate_Regression.txt",append=T)
write("0.5-0.7 = Poor", file="Multivariate_Regression.txt",append=T)
write("0.7-0.8 = Acceptable", file="Multivariate_Regression.txt",append=T)
write("0.8-0.9 = Excellent", file="Multivariate_Regression.txt",append=T)
write("0.9-1.0 = Outstanding", file="Multivariate_Regression.txt",append=T)

#### RESIDUALS ANALYSIS (BIC) ####

BIC_resid<-residuals(BIC_glm)
bicresfile<-capture.output(summary(BIC_resid))

Pearson_resid<-residuals(BIC_glm, "pearson")
pearresfile<-capture.output(summary(Pearson_resid))

# residuals correlogram

site_coor<-crds(Sites, df=TRUE, list=FALSE)
nosite_coor<-crds(No_Sites, df=TRUE, list=FALSE)
coor<-as.matrix(rbind(site_coor,nosite_coor))

library(pgirmess)

resid_correl<-data.frame(correlog(coor,BIC_resid,method="Moran"))
nrow(resid_correl)
plot(correlog(coor,BIC_resid,method="Moran"))

#### MODEL VALIDATION ####

# Data Upload

predsurf<-read_RAST("Probability_Surface")

detach("package:MASS", unload = TRUE)

sitespred<-extract(predsurf, Sites)
sitespred<-select(sitespred,-ID)
colnames(sitespred)<-"pred_mod"

nositespred<-extract(predsurf, No_Sites)
nositespred<-select(nositespred,-ID)
colnames(nositespred)<-"pred_mod"

# Kolmogorov-Smirnov Test of predictive values

library(Matching)

kspredfile<-capture.output(summary(ks.boot(sitespred,nositespred,10000)))

# Gain Calculation

a<-length(sitespred[sitespred>0.5])
b<-length(sitespred[sitespred<0.5]) 
c<-length(nositespred[nositespred>0.5])
d<-length(nositespred[nositespred<0.5])
n<-a+b+c+d
n1<-a+b
n2<-c+d
t<-(n*((a*d)-(b*c))^2)/((n1*n2)*(a+c)*(b+d))
pval<-pchisq(t, df=1, lower.tail=F)


# Upload raster with predictive categories

Cat_Pred<-read_RAST(c("Prob01","Prob02","Prob03",
                      "Prob04","Prob05","Prob06",
                      "Prob07","Prob08","Prob09"))

# Kvamme's Gain

detach("package:Matching", unload = TRUE)
detach("package:MASS", unload = TRUE)

# sites in every predictive step
sitespredcat<-extract(Cat_Pred, Sites)
sitespredcat<-select(sitespredcat, -ID)
totsites<-colSums(sitespredcat)

# sites presence percentage for every predictive step
percsit<-round((totsites/length(Sites)),digits=2)

# prediction area by predictive steps
areapredcat<-freq(Cat_Pred,bylayer=T,usenames=TRUE)
areapredcat<-areapredcat %>% 
  group_by(layer) %>% 
  mutate(area = sum(count)) %>%
  filter (value > 0) %>%
  mutate (percarea = round(count/area,2))
percarea<-areapredcat$percarea

# Kvamme's gain calculation
Sites_KG<-(1-(percarea/percsit))
KG_sites<-rbind(percarea,percsit,Sites_KG)
colnames(KG_sites)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                      ">0.8",">0.9")
row.names(KG_sites)<-c("% Area","% Siti","Kvamme's Gain")
round(KG_sites, digits=2)

# Cumulative Prediction plot

# sites
percsit["Prob0"]<-1
percsit["Prob10"]<-0
percfull<-as.data.frame(percsit)
percprob<-rbind(percfull[10,],percfull[1,],percfull[2,],percfull[3,],
                percfull[4,],percfull[5,],percfull[6,],percfull[7,],
                percfull[8,],percfull[9,],percfull[11,])

# no sites
percnono<-round((right_ns/length(No_Sites)),digits=2) 
percnono
percnono["Prob0"]<-1
percnono["Prob10"]<-0
percfullns<-as.data.frame(percnono)
percprobns<-rbind(percfullns[11,],percfullns[1,],percfullns[2,],percfullns[3,],
                  percfullns[4,],percfullns[5,],percfullns[6,],percfullns[7,],
                  percfullns[8,],percfullns[9,],percfullns[10,])

# table
cutoff<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

tabperc<-cbind(cutoff,(percprob*100),(percprobns*100))
colnames(tabperc)<-c("cutoff","%sites","%random")

# plot

jpeg("Correct_Prediction.jpg",width=2000,height=2000,units="px",
     res=300)
par(oma= c(1, 0, 1, 0))
plot(tabperc[,2]~tabperc[,1],xlab="predicted site probability",ylab="percent
correct predictions",main="Correct Prediction %",sub="Sites in red & Non-Sites in black",
     col="red",pch = 19,axes=F,font.sub=3)
lines(tabperc[,2]~tabperc[,1],col="red")
par(new=TRUE)
plot(tabperc[,3]~tabperc[,1],xlab="",ylab="",pch = 19,axes=F)
lines(tabperc[,3]~tabperc[,1])
abline(v=0.3,lty=2,col="grey")
Axis(side = 1, at = seq(0,1,by = 0.1), labels=seq(0,1,by = 0.1), cex=0.2)
Axis(side = 1, at = seq(0,1,by = 0.05),labels=F)
Axis(side = 2, at = seq(0,100,by = 10), labels=seq(0,100,10),las=1)
Axis(side = 2, at = seq(0,100,by = 5), labels=F)
dev.off()

#### MODEL ASSESSEMENT WITH FIELD SURVEY SITES####

# Data Upload
setwd("D:/Archeologia/Carona/Analisi/Modelli/Modello_2023")

Survey<-(read_VECT("SBC_Survey_2023"))

# Sites by typology

mod.val<-function(x){
  classes<-list()
  y<-levels(factor(x$tipo))
  for(i in y){
    classes[[i]]<-subset(x,x$tipo==i)
  }
  count<-sapply(classes,function(x){(length(x$pred_mod))})
  pred<-sapply(classes,function(x){
    (length(x$pred_mod[x$pred_mod>0.30]))/(length(x$pred_mod))})
  unpred<-1-pred
  tab<-rbind(round(count,0),round(pred,2),round(unpred,2))
  row.names(tab)<-c("n. siti",">0.30","<0.30")
  count<-length(x$pred_mod)
  all<-rbind((length(x$pred_mod)),(length(x$pred_mod[x$pred_mod>0.30]))/
               (length(x$pred_mod)),
             1-((length(x$pred_mod[x$pred_mod>0.30]))/
                  (length(x$pred_mod))))
  colnames(all)<-"Totale"
  tab<-cbind(tab,round(all,2))
  return(tab)
}

surveycat<-mod.val(Survey)

# Subset Single Typology
Survey<-as.data.frame(Survey)
resto<-subset(Survey,tipo=="Altro")
bait<-subset(Survey,tipo=="Baita")
carbone<-subset(Survey,tipo=="Carbonaia")
enel<-subset(Survey,tipo=="Idroelettrico")
indef<-subset(Survey,tipo=="Indefinito")
muro<-subset(Survey,tipo=="Muro")
pont<-subset(Survey,tipo=="Ponte")
recinto<-subset(Survey,tipo=="Recinto")
ricovero<-subset(Survey,tipo=="Ricovero")
riparo<-subset(Survey,tipo=="Riparo")
rocco<-subset(Survey,tipo=="Roccolo")
spiazzo<-subset(Survey,tipo=="Spiazzo")
stalla<-subset(Survey,tipo=="Stalla")
milit<-subset(Survey,tipo=="Militare")
minaltro<-subset(Survey,tipo=="Minerario")
miniera<-subset(Survey,tipo=="Miniera")
reglana<-subset(Survey,tipo=="Reglana")
natrip<-subset(Survey,tipo=="Riparo Naturale")
viab<-subset(Survey,tipo=="Sentieri")
smarino<-subset(Survey,tipo=="Smarino")  

# Kolmogorov-Smirnov test (all sites)
predsurvey<-data.frame(Survey$pred_mod)

library(Matching)

# x=valore predittivo punti del survey
# y=superficie predittiva
# z=maschera
# nsim=ripetizione


ks.survey<-function(x,y,z,nsim){
  randpts<-spatSample((mask(y,z,inverse=T)),(100*nrow(x))/40,method='random',
                                   replace=F,na.rm=T,values=TRUE,as.df=TRUE)
  return(summary(ks.boot(x,randpts,nboots=nsim)))
}

ks_survey<-capture.output(ks.survey(predsurvey,predsurf,Survey_vect,10000))

# Kolmogorov-Smirnov test (by typology)

# x=punti del survey (df)
# y=superficie predittiva
# z=maschera (tutti i siti del survey)
# nsim=ripetizione

library(Matching)

class.survey<-function(x,y,z,nsim){
  predval<-list()
  b<-levels(factor(x$tipo))
  randpts<-spatSample((mask(y,z,inverse=T)),(100*nrow(x))/40,method='random',
                      replace=F,na.rm=T,values=TRUE,as.df=TRUE)
  for(i in b){
    sur<-subset(x$pred_mod,x$tipo==i)
    predval[[i]]<-ks.boot(sur,randpts,nboots=nsim)
  }
  return(lapply(predval,summary))
}

ks_surclass<-capture.output(class.survey(Survey,predsurf,Survey_vect,10000))

# Kvamme's Gain Calculation (all sites)

detach("package:Matching", unload = TRUE)
detach("package:MASS", unload = TRUE)

# survey sites in every predictive step
surveypredcat<-extract(Cat_Pred, Survey_vect)
surveypredcat<-select(surveypredcat, -ID)
surveypredcat[is.na(surveypredcat)]<-0
totsurvey<-colSums(surveypredcat)

# sites presence percentage for every predictive step
percsur<-round((totsurvey/nrow(Survey)),digits=2)


# prediction area by predictive steps
areasurpc<-freq(Cat_Pred,bylayer=T,usenames=TRUE)
areasurpc<-areasurpc %>% 
  group_by(layer) %>% 
  mutate(area = sum(count)) %>%
  filter (value > 0) %>%
  mutate (percarea = round(count/area,2))
percareas<-areasurpc$percarea

# calculation
Survey_KG<-(1-(percareas/percsur))

KG_survey<-rbind(percareas,percsur,Survey_KG)
colnames(KG_survey)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                      ">0.8",">0.9")
row.names(KG_survey)<-c("% Area","% Siti","Kvamme's Gain")
round(KG_survey, digits=2)

# Kvamme's Gain calculation (by typology)

# x= survey (vect)
# y= catgorized predmod maps
# z= survey (d.f.)

kg.survey<-function(x,y){
  KGval<-list()
  b<-levels(factor(x$tipo))
  for(i in b){
    sub<-subset(x,x$tipo==i)
    predcat<-extract(y, sub, ID=FALSE)
    predcat[is.na(predcat)]<-0
    percsur<-round(((colSums(predcat))/nrow(sub)),digits=2)
    areapredcat<-freq(y,bylayer=T,usenames=TRUE)
    areapredcat<-areapredcat %>% 
      group_by(layer) %>% 
      mutate(area = sum(count)) %>%
      filter (value > 0) %>%
      mutate (percareas = round(count/area,2))
    percareas<-areapredcat$percareas
    KG<-(1-(percareas/percsur))
    TabKG<-rbind(percareas,percsur,KG)
    colnames(TabKG)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                       ">0.8",">0.9")
    row.names(TabKG)<-c("% Area","% Siti","Kvamme's Gain")
    KGval[[i]]<-round(TabKG, digits=2)
  }
  return(KGval)
}

kg_catsur<-kg.survey(Survey_vect,Cat_Pred)

#end
