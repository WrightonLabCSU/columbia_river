##### Load relevant libraries
library(ggplot2)
library(pls)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Google Drive/University/Wrighton Lab; PhD/Viral Predation, More than meets the eye/paper_figures/figure_xx_viral_hosts/4.SPLS/rarefied_to_deep_host_and_virus_111_virus") #whereve you are going to work. VIP has to be here.

source("/Volumes/Macintosh HD/Users/josue.rodriguez/Google Drive/University/Wrighton Lab; PhD/Viral Predation, More than meets the eye/paper_figures/figure_xx_viral_hosts/4.SPLS/VIP.R") ## This is a custom script  to compute VIPs, has to be in the same directory as this script (or this path has to be changed)


##### Load OTU data. In my case, this contains all viruses with non zero abundances.

all_data_vOTU<-read.csv("host_virus_or_both.csv",header=T) #Input is abundance tables for either the hosts, the viruses, or a combination of both.

# By default, R doesn't remove the first column if it has an ID, so we do it instead
row.names(all_data_vOTU)<-all_data_vOTU[,1]
all_data_vOTU<-all_data_vOTU[,-1] 
all_data_vOTU<-as.matrix(all_data_vOTU) # force R to see it as a matrix
dim(all_data_vOTU) # this gives the matrix dimension, so 28 vOTUs and 30 samples

## now we transpose the matrix, because sPLS wants the samples as rows
all_data_vOTU<-t(all_data_vOTU) #transpose

##### Load metabolites/geochemistry data. Mine has 31 samples and these are rows.
metabolite<-read.csv("/Volumes/Macintosh HD/Users/josue.rodriguez/Google Drive/University/Wrighton Lab; PhD/Viral Predation, More than meets the eye/paper_figures/figure_xx_viral_hosts/4.SPLS/rarefied_to_deep_host_and_virus_111_virus/10_samples_geochemistry.csv",header=T)

#same as before, write in the row names correctly.
row.names(metabolite)<-metabolite[,1]
metabolite<-metabolite[,-1]
metabolite<-as.matrix(metabolite) # force R to see it as a matrix
dim(metabolite) # this gives the matrix dimension, so 19 metabolites and 33 samples


#Step 2: Predict values

th_r2<-0.1 # We will only look at the PLS if the correlation is better than 0.1
for (i in 1:ncol(metabolite)){ # We treat each metabolite independently
  parameter<-colnames(metabolite)[i]
  obs_values<-metabolite[,i] # these are the observed values we'll try to predict
  print(paste("Trying to predict ",parameter," --- ",i,sep=""))
  # We perform the sPLS, trying to predict our specific metabolite vector (metabolite[,i]) using our whole OTU table. Validation is LOO, so "Leave one out", i.e. we train model on n-1 samples and try to predict the value for the remaining one. the "method" argument is to chose the correct type of sPLS
  pls_result<-plsr(obs_values ~ all_data_vOTU, validation="LOO",method="oscorespls") 
  # Now we check the vector of r2 (sPLS tries to use different numbers of OTUs and provides a correlation between predicted and observed for each of them, so we get a vectore of r2 and not just one r2 value)
  r2_vector<-R2(pls_result)
  max<-0
  max_comp<--1
  for (j in 1:length(r2_vector$val)){
    if(!(is.na(r2_vector$val[j]))){
      if(r2_vector$val[j]>th_r2){
        if(r2_vector$val[j]>max){
          max<-r2_vector$val[j]
          max_comp<-r2_vector$comp[j]
        }
      }
    }
  }
  print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep=""))
}

# So here we print the highest r2 across all predictions.


#Step 3: Plotting the results

## So now we can look at the metabolites, and see if any can be predicted by the OTUs abundance
## We  set i to the corresponding column value based on results from above
i<-3 #this number corresponds to the "iteration" of the prediction. E.g. prediction 1, prediction 2, prediction 3, etc. In my case, the value with a significant SPLS value was the 6th, "CN". Change it to whatever yours is.

# And we regenerate the corresponding results
parameter<-colnames(metabolite)[i]
obs_values<-metabolite[,i] # these are the observed values we'll try to predict
pls_result<-plsr(obs_values ~ all_data_vOTU, validation="LOO",method="oscorespls") 
r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(!(is.na(r2_vector$val[j]))){
    if(r2_vector$val[j]>th_r2){
      if(r2_vector$val[j]>max){
        max<-r2_vector$val[j]
        max_comp<-r2_vector$comp[j]
      }
    }
  }
}


# Plotting predicted vs observed
df<-data.frame(x=obs_values,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_","-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted -- r2=",max)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()

# Next we checking the VIP (variable importance in projection), and output a table of the 100 highest values, this will tell us which OTU you need to know the abundance of to correctly predict the feature of interest
output<-paste("VIP_values_with_",parameter,".csv",sep="")
cat("Rank,OTU,VIP\n",file = output,append=FALSE)
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
for (k in 1:100){
  cat(paste(k,names(vip_components[k]),vip_components[k],"\n",sep=","),file=output,append=TRUE)
}
## Check the correlation between predicted and modeled (should be consistent with what plsr gave us)
cor.test(df$x,df$y)
## Alternatively, we can also use the built-in function "predplot" but I find it less pretty
## Can be good to double check the ggplot2 plot though (should be the same)
predplot(pls_result,ncomp=max_comp)

## Now we can also plot individual OTUs vs n_per based on the VIP list
for (k in 1:5){
  OTU<-unlist(names(vip_components[k]))
  vec_OTU<-all_data_vOTU[,OTU]
  df<-data.frame(x=vec_OTU,y=obs_values)
  print(ggplot(data=df) + geom_point(aes(x=x,y=y)) + xlab(OTU) + ylab(paste("Measured ",parameter,sep="")) + ggtitle(paste("Comparison of ",OTU," vs measured ",parameter)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")))
}
