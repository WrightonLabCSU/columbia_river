#3D NMDS for continuous variables
#This does not work well with NA values. I think I can make it ignore them but for now I will remove them.

library(readxl)
library(ggplot2)
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
library(vegan)

##read in feature table with species as columns and samples as rows

virusmag<-read.table('7.4_NMDS_29samp_vMAG_metaP_pergeonme_presabs.txt', header=T, sep="\t", check.names=TRUE)
rownames(virusmag)<-virusmag[,1]
virusmag<-virusmag[,-1]


##read in chemistry (nona)
chem = read.csv('/7.5_29_samples_geochemistry.csv', sep = ',', header = TRUE, check.names = T)
rownames(chem)=chem[,1]
chem = chem[,-1]

##NMDS in 3D (k = 3 gives 3D)
NMDS_Bray_virusmag <-metaMDS(virusmag, distance = "bray",
                          autotransform = F, noshare = 0.1, trace = 1, trymax = 500)
#sppscores(NMDS_Bray_virusmag)=virusmag

ord.virusmag = as.data.frame(scores((NMDS_Bray_virusmag), display="sites"))
ord.virusmag$sampleid=row.names(ord.virusmag)


plot(NMDS_Bray_virusmag)
stressplot(NMDS_Bray_virusmag)

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

#Run mrpp and ANOSIM on Depth between samples (BROKEN)
#This establishes the coordinates for the samples and all OTUS onto an non-metric dimensional scaling
Ord_dist <-metaMDSdist(virusmag, distance = "bray", autotransform = FALSE, noshare = 0.1, trace = 1)

mrpp(virusmag, chem$Depth, permutations=999, distance="bray")
anosim(Ord_dist, chem$Depth, permutations = 999)

mrpp(virusmag, chem$Site, permutations=999, distance="bray")
anosim(Ord_dist, chem$Site, permutations = 999)


################################## By Site
##plot 3D NMDS in plotly
site = plot_ly(ord.virusmag, x = ~NMDS1, y= ~NMDS2, color = ~chem$Site, colors =c("#66c2a5","#fc8d62"), size=10 )%>%
  add_markers(text = ~paste(sampleid)) %>%
  layout(
    title='Site',
    scene = list(xaxis = list(title = 'NMDS1'),
                 yaxis = list(title = 'NMDS2')))
site

################################## By Depth
depth = plot_ly(ord.virusmag, x = ~NMDS1, y= ~NMDS2, color = ~chem$Depth, colors =c("#f1eef6","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#034e7b"), size = 10 )%>%
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
############################################################################################
