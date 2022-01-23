#from lat_long_trebouxia21_MODIFIED_pmcUNIQUECLADESE_upd.R

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

######################
#ADD Sampling Error
######################
#read in data file
missing<-read.csv(file="/xxxx/MASTER_accession_list_INFOADDED_MODIFIED.csv",stringsAsFactors=FALSE)

#retain seqs that passed criteria
keeps<-missing[!is.na(missing$OTU_97.5_Raw_28_Oct_2019_Modified) &!is.na(missing$LatC),]

#number passing quality filters
nrow(keeps)

#looks like all w coords have climatic data
keeps[!is.na(keeps$LatC) & is.na(keeps$bio1),]

#might not be same for biomes (some could be close to water (looks like 72 accessions lack biome info)
nrow(keeps[!is.na(keeps$LatC) & is.na(keeps$Biome),])

#make temp variables in degrees C
cols<-paste("bio",1:11,sep="")
for(x in 1:nrow(keeps)){
	for(b in 1:length(cols)){
		keeps[x,cols[b]]<-keeps[x,cols[b]]/10
	}
}

keeps.f<-keeps[,c(paste("bio",1:19,sep=""),"OTU_97.5_Raw_28_Oct_2019_Modified")]
rownames(keeps.f)<-keeps$Organism_Acc_No
keeps.f$OTU_97.5_Raw_28_Oct_2019_Modified<-as.factor(keeps.f$OTU_97.5_Raw_28_Oct_2019_Modified)


#pooled variance function from Liam Revell: http://blog.phytools.org/2016/08/simple-function-to-compute-pooled.html
vars<-function(x) mean(x^2)-mean(x)^2

pooled.var<-function(y,x){
    k<-length(levels(x))
    N<-summary(x)
    sum(sapply(levels(x),function(l,y,x) vars(y[x==l]),y=y,x=x)*(N-1))/(sum(N)-k)
}

#calculate pooled variance for each variable
se.bc<-rep(NA,19)
names(se.bc)<-paste("SE_BIO",1:19,sep="")

for(x in 1:19){
	se.bc[x]<-pooled.var(y=keeps.f[,paste("bio",x,sep="")],x=keeps.f$OTU_97.5_Raw_28_Oct_2019_Modified)
}

#then follow Slater & Frischia (2019) and Garamszegi (2014) and compute sampling error as being the pooled variance/sqrt(sample size) for each OTU
taxa<-read.csv(file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv",stringsAsFactors=FALSE,row.names=1)
taxa[,(paste("SE_","Bio",1:19,sep=""))]<-NA

for(x in 1:19){
	for(z in 1:nrow(taxa)){
		taxa[z,(paste("SE_","Bio",x,sep=""))]<-sqrt(se.bc[x])/sqrt(taxa[z,"No_Sequences"])
	}
}

write.csv(taxa,file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES_wSE_corrected.csv",row.names=TRUE)

######################
#Adjust MCC tree (ladderize) for some analyses
######################

#read in ladderized MCC tree
tr<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre")
tr<-ladderize(tr,FALSE)
write.tree(tr,file="/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED_ladderized.tre")
