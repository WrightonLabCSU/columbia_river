library(ggplot2)
library(vegan)
library(grid)
library(MASS)
library(ggsn)
library(RColorBrewer)
library(extrafont)

#For both of these data sets, the rows need to be samples
#columns need to be either OTUS for virus, or variables for CHEM



#################
#Normalized, JGI only for viruses first


setwd('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/Columbia_River_Hyporheic_Zone/NMDS/v2')

virus<-read.table('viral_genomes.txt', header=TRUE, sep="\t", check.names=TRUE)

host<-read.table('host_genomes.txt', header=TRUE, sep="\t", check.names=TRUE)

CHEM<-read.table('28samp_geochem.txt', header=TRUE, sep="\t", check.names=TRUE)

CHEMvar<-read.table('28samp_geochem_variables.txt', header=TRUE, sep="\t", check.names=TRUE)


###############MRPP and ANOSIM
CHEMvar$core[which(CHEMvar$Site==0)] <- "N"
CHEMvar$core[which(CHEMvar$Site==1)] <- "S"
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==1)] <- '0_10'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==2)] <- '10_20'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==3)] <- '20_30'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==4)] <- '30_40'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==5)] <- '40_50'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==6)] <- '50_60'
CHEMvar$depthrange_cm[which(CHEMvar$depthrange_cm==1.5)] <- '0_30'
CHEMvar$core <- as.factor(CHEMvar$core)

#establish samples as rownames. These need to be in the exact same order, otherwise the scaling of the vectors will be wrong.
rownames(virus)<-virus[,1]
virus<-virus[,-1]
#virus<-t(virus)
rownames(host)<-host[,1]
host<-host[,-1]
#host<-t(host)
rownames(CHEM)<-CHEM[,1]
CHEM<-CHEM[,-1]
#CHEM<-t(CHEM)
rownames(CHEMvar)<-CHEMvar[,1]
CHEMvar<-CHEMvar[,-1]
#CHEMvar<-t(CHEMvar)

log.chem<-log(CHEM) #log of chem values
sqrt.chem<-sqrt(CHEM) #sqrt of chem values

#viral 
Ord_dist_vir <-metaMDSdist(virus, distance = "bray", autotransform = T, noshare = 0.1, trace = 1, trymax=500)
mrpp(virus, CHEMvar$depthrange_cm, permutations=999, distance="bray")
anosim(Ord_dist_vir, CHEMvar$depthrange_cm, permutations = 999)

mrpp(virus, CHEMvar$core, permutations=999, distance="bray")
anosim(Ord_dist_vir, CHEMvar$core, permutations = 999)

#host stats
Ord_dist_host <-metaMDSdist(host, distance = "bray", autotransform = FALSE, noshare = 0.1, trace = 1)
mrpp(host, CHEMvar$depthrange_cm, permutations=999, distance="bray")
anosim(Ord_dist_host, CHEMvar$depthrange_cm, permutations = 999)

mrpp(host, CHEMvar$core, permutations=999, distance="bray")
anosim(Ord_dist_host, CHEMvar$core, permutations = 999)


#to size the dots based on the PC1 of all the chemical data (like a bubble plot in primer), let's apply a PCA scale
chem.pca_log<-prcomp(na.omit(log.chem), center=TRUE, scale.=TRUE)
chem.pca_sqrt<-prcomp(na.omit(sqrt.chem), center=TRUE, scale.=TRUE)


#In order to use the PC1 scale, let's establish a variable with all of the PCA scale coordinates by sample
chem.variables_log<-predict(chem.pca_log)
chem.variables_sqrt<-predict(chem.pca_sqrt)

chem.variables.PC1_log<-chem.variables_log[,1]
chem.variables.PC1_sqrt<-chem.variables_sqrt[,1]

#This establishes the coordinates for the samples and all vOTUS onto an non-metric dimensional scaling
ordvir<-metaMDS(virus)
ordvir # the stress for this should be around <.1 
stressplot(ordvir) #scatter seems close to the line, demonstrating that it is a good fit. All good!

ordhost<-metaMDS(host)
ordhost # the stress for this should be around <.1 
stressplot(ordhost) #scatter seems close to the line, demonstrating that it is a good fit. All good!



#envfit here will fit the vectors onto the ordination established above
fit_logvir <- envfit(ordvir, log.chem, perm = 999, na.rm=TRUE)
scores(fit_logvir, "vectors")

fit_sqrtvir <- envfit(ordvir, sqrt.chem, perm = 999, na.rm=TRUE)
scores(fit_sqrtvir, "vectors")

fit_loghost <- envfit(ordhost, log.chem, perm = 999, na.rm=TRUE)
scores(fit_loghost, "vectors")

fit_sqrthost <- envfit(ordhost, sqrt.chem, perm = 999, na.rm=TRUE)
scores(fit_sqrthost, "vectors")


#Lets plot this first to make sure it works, just in vegan
plot(ordvir)
plot(fit_logvir)
plot(fit_logvir, p.max = 0.05, col = "green")

plot(ordvir)
plot(fit_sqrtvir)
plot(fit_sqrtvir, p.max = 0.05, col = "green")

plot(ordhost)
plot(fit_sqrthost)
plot(fit_sqrthost, p.max = 0.05, col = "green")

plot(ordhost)
plot(fit_loghost)
plot(fit_loghost, p.max = 0.05, col = "green")



##Can also plot it like this with the names so that it looks nicer. -JRR
ordiplot(ordvir,type="points")
orditorp(ordvir,display="sites",cex=1,air=.2)
plot(fit_logvir)
plot(fit_logvir, p.max = 0.05, col = "green")

ordiplot(ordhost,type="points")
orditorp(ordhost,display="sites",cex=1,air=.2)
plot(fit_loghost)
plot(fit_loghost, p.max = 0.05, col = "green")










####################################################
####################################################
####################################################
#Now to do a PROCRUSTES analysis on the two MDS that were created

vare.proc<-procrustes(ordvir, ordhost, symmetric = FALSE)
vare.proc
summary(vare.proc)
plot(vare.proc, kind=1, type='text')
plot(vare.proc, kind=2)
residuals(vare.proc)

protest(ordvir, ordhost, scores = 'sites', symmetric=FALSE, permutations=10000)

