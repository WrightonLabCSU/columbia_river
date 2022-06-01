#3D NMDS for continuous variables
#This does not work well with NA values. I think I can make it ignore them but for now I will remove them.

library(readxl)
library(ggplot2)
library(vegan)
library(vegan3d)
library(grid)
library(MASS)
library(RColorBrewer)
library(purrr)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(mds)

##read in feature table with species as columns and samples as rows

hostmag<-read.table('7.2_NMDS_29samp_MAG_metaP_pergeonme_presabs.txt', header=T, sep="\t", check.names=TRUE)
rownames(hostmag)<-hostmag[,1]
hostmag<-hostmag[,-1]


##read in chemistry (nona)
chem = read.csv('7.5_29_samples_geochemistry.csv', sep = ',', header = TRUE, check.names = T)
rownames(chem)=chem[,1]
chem = chem[,-1]

log_chem=log(chem[,3:19])

#Make the variables in the chem file characters so ggplot2 doesn't yell.
chem$nh4 <- as.character(log_chem$nh4)
chem$cn <- as.character(log_chem$cn)
chem$ca <- as.character(log_chem$ca)
chem$n_per <- as.character(log_chem$n_per)
chem$al <- as.character(log_chem$al)
chem$s_per <- as.character(log_chem$s_per)

##NMDS in 3D (k = 3 gives 3D)
NMDS_Bray_hostmag <-metaMDS(hostmag, distance = "bray", k =2,
                          autotransform = FALSE, noshare = 0.1, trace = 1, trymax = 500)
ord.hostmag = as.data.frame(scores((NMDS_Bray_hostmag), display="sites"))
ord.hostmag$sampleid=row.names(ord.hostmag)

stressplot(NMDS_Bray_hostmag)

chem$Site[which(chem$Site==0)] <- "N"
chem$Site[which(chem$Site==1)] <- "S"
chem$Depth[which(chem$Depth==1)] <- '.0-10'
chem$Depth[which(chem$Depth==2)] <- '.10-20'
chem$Depth[which(chem$Depth==3)] <- '.20-30'
chem$Depth[which(chem$Depth==4)] <- '.30-40'
chem$Depth[which(chem$Depth==5)] <- '.40-50'
chem$Depth[which(chem$Depth==6)] <- '.50-60'
chem$Depth[which(chem$Depth==7)] <- '.0-30'
chem$Site <- as.factor(chem$Site)

#This establishes the coordinates for the samples and all OTUS onto an non-metric dimensional scaling
Ord_dist <-metaMDSdist(hostmag, distance = "bray", autotransform = FALSE, noshare = 0.1, trace = 1)

p<-ggplot(ord.hostmag)+geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour=chem$Site))
p

q<-ggplot(ord.hostmag)+geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour=chem$Depth))
q

#Run mrpp and ANOSIM on Depth between samples
mrpp(hostmag, chem$Depth, permutations=999, distance="bray")
anosim(Ord_dist, chem$Depth, permutations = 999)

mrpp(hostmag, chem$Site, permutations=999, distance="bray")
anosim(Ord_dist, chem$Site, permutations = 999)


################################## By Site
##plot 3D NMDS in plotly
site = plot_ly(ord.hostmag, x = ~NMDS1, y= ~NMDS2, color = ~chem$Site, colors =c("#66c2a5","#fc8d62"), size=10 )%>%
  add_markers(text = ~paste(sampleid)) %>%
  layout(
    title='Site',
    scene = list(xaxis = list(title = 'NMDS1'),
                 yaxis = list(title = 'NMDS2')))
site

################################## By Depth
depth = plot_ly(ord.hostmag, x = ~NMDS1, y= ~NMDS2, color = ~chem$Depth, colors =c("#f1eef6","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#034e7b"), size=10 )%>%
  add_markers(text = ~paste(sampleid)) %>%
  layout(
    title='Depth',
    scene = list(xaxis = list(title = 'NMDS1'),
                 yaxis = list(title = 'NMDS2')))
depth

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################