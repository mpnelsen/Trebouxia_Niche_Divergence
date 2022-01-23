#from change_rates4.R

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

require(phytools)
require(phangorn)
require(geiger)


#function to get age of tip or clade from Nelsen et al. 2020a Geobiology
#when a single tree
get.age<-function(tree,clade,position){
	node.depth.edgelength(tree)->node.ages
	max(node.ages)->depth
	depth-node.ages->new.node.ages
	if(length(clade)==1){
		tip.number<-which(tree$tip.label==clade)
		parent.node<-Ancestors(tree,tip.number,type="parent")
		age<-new.node.ages[parent.node]
	} else if(length(clade)>1){
		crown.node<-findMRCA(tree,clade)
		if(position=="crown"){
			age<-new.node.ages[crown.node]
			} else if(position=="stem"){
				stem.node<-Ancestors(tree,crown.node,type="parent")
				if(stem.node==0){
					age<-NA
					} else if(stem.node!=0){
						age<-new.node.ages[stem.node]
			}
		}
	}
return(age)
}



#identify path
path<-"/xxxx/"

#identify variables of interest
#voi<-c(1,5,6,12,16,17)
voi<-c(1,5,6)

#read modern data (Means)
taxa<-read.csv(file=paste(path,"OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE,row.names=1)

#Add columns for Max/min
taxa[,paste("Max_Bio",voi,sep="")]<-NA
taxa[,paste("Min_Bio",voi,sep="")]<-NA

#Get Max and Min Observed for OTUs from Nelsen et al. 2021 FEMS Microbiology Ecology
dat<-read.csv(file=paste(path,"MASTER_accession_list_INFOADDED_MODIFIED.csv",sep=""))

#subset to retain OTUs of interest
nrow(dat)
dat<-dat[dat$OTU_97.5_Raw_28_Oct_2019_Modified %in% rownames(taxa),]
dat<-dat[!is.na(dat$bio1),]
nrow(dat)
sum(taxa$No_Sequences)

for(x in 1:nrow(taxa)){
	subdat<-NA
	subdat<-dat[dat$OTU_97.5_Raw_28_Oct_2019_Modified %in% rownames(taxa)[x],]
	for(v in 1:length(voi)){
		taxa[x,paste("Min_Bio",voi[v],sep="")]<-min(subdat[paste("bio",voi[v],sep="")])
		taxa[x,paste("Max_Bio",voi[v],sep="")]<-max(subdat[paste("bio",voi[v],sep="")])
	}
}

#put on degrees celsius
taxa$Min_Bio1<-taxa$Min_Bio1/10
taxa$Min_Bio5<-taxa$Min_Bio5/10
taxa$Min_Bio6<-taxa$Min_Bio6/10
taxa$Max_Bio1<-taxa$Max_Bio1/10
taxa$Max_Bio5<-taxa$Max_Bio5/10
taxa$Max_Bio6<-taxa$Max_Bio6/10

#read mcc tree
mcc<-read.nexus(file=paste(path,"burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre",sep=""))

#read 1k trees
tr1k<-read.nexus(file=paste(path,"burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees",sep=""))

#for each tree
for(t in 1:1001){
	print(t)
	#subset tree of interest
	tr<-NULL
	trtr<-NULL
	if(t %in% 1:1000){
		tr<-tr1k[[t]]
	}
	if(t==1001){
		tr<-mcc
	}
	#go through variables of interest
	#for(var in c(1,5,6,12,16,17)){
	for(var in c(1,5,6)){
		#print(var)
		toplot<-NULL
		temp<-NULL
				
		#read results of previous fitContinuous analyses to identify best fitting model and get parameter estimates [found in Model_Fitting_Fig3 folder]
		toplot<-read.csv(file=paste(path,"xxxx/","bio",var,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
		
		#make a data frame to store info
		cols<-c("ModMeanVal","ModMinVal","ModMaxVal","Age","NodeNo","AncVal","ValChangeMean","AbsValChangeMean","RateMean","AbsRateMean","ValChangeMin","AbsValChangeMin","RateMin","AbsRateMin","ValChangeMax","AbsValChangeMax","RateMax","AbsRateMax")
		temp<-data.frame(matrix(nrow=nrow(taxa),ncol=c(length(cols))))
		rownames(temp)<-rownames(taxa)
		colnames(temp)<-cols

		#for each tip, get: modern mean value, max value and min value
		temp$ModMeanVal<-taxa[,paste("Mean_Bio",var,sep="")]
		temp$ModMaxVal<-taxa[,paste("Max_Bio",var,sep="")]
		temp$ModMinVal<-taxa[,paste("Min_Bio",var,sep="")]

		#for each tip, get: stem age
		for(z in 1:nrow(temp)){
			temp$Age[z]<-get.age(tr,rownames(temp)[z],"stem")
		}

		#transform the tree using best fit model of evolution previously identified
		if(toplot$Best[t] %in% "BM"){
			trtr<-rescale(tr,model="BM",sigsq=toplot$BM_sigsq[t])
		}
		if(toplot$Best[t] %in% "OU"){
			trtr<-rescale(tr,model="OU",alpha=toplot$OU_a[t],sigsq=toplot$OU_sigsq[t])
		}
		if(toplot$Best[t] %in% "EB"){
			trtr<-rescale(tr,model="EB",a=toplot$EB_a[t],sigsq=toplot$EB_sigsq[t])
		}

		#for each tip, get: mrca node number
		for(z in 1:length(trtr$tip.label)){
			temp$NodeNo[rownames(temp) %in% trtr$tip.label[z]]<-getParent(trtr,z)
		}

		#use fastAnc phytools or reconstruct ape
		vals<-NULL
		vals<-temp$ModMeanVal
		names(vals)<-rownames(temp)
		rec<-fastAnc(trtr,vals)

		#get new values
		for(z in 1:nrow(temp)){
			temp$AncVal[z]<-rec[names(rec) %in% temp$NodeNo[z]][[1]]
			
			#mean diffs
			temp$ValChangeMean[z]<-temp$AncVal[z]-temp$ModMeanVal[z]
			temp$AbsValChangeMean[z]<-abs(temp$ValChangeMean[z])
			temp$RateMean[z]<-temp$ValChangeMean[z]/temp$Age[z]
			temp$AbsRateMean[z]<-temp$AbsValChangeMean[z]/temp$Age[z]
			
			#Max val diffs
			temp$ValChangeMax[z]<-temp$AncVal[z]-temp$ModMaxVal[z]
			temp$AbsValChangeMax[z]<-abs(temp$ValChangeMax[z])
			temp$RateMax[z]<-temp$ValChangeMax[z]/temp$Age[z]
			temp$AbsRateMax[z]<-temp$AbsValChangeMax[z]/temp$Age[z]
			
			#Min val diffs
			temp$ValChangeMin[z]<-temp$AncVal[z]-temp$ModMinVal[z]
			temp$AbsValChangeMin[z]<-abs(temp$ValChangeMin[z])
			temp$RateMin[z]<-temp$ValChangeMin[z]/temp$Age[z]
			temp$AbsRateMin[z]<-temp$AbsValChangeMin[z]/temp$Age[z]
		}
		write.csv(temp,paste(path,"/xxxx/","BIO",var,"_Tree",t,"_RateChange.csv",sep=""))
	}
}	


#OK...now to summarize
#for each variable
require(coda)
for(var in c(1,5,6)){
	print(var)
	resmean<-NULL
	trresmean<-NULL
	#read results from the MCC tree (1001)
	resmean<-read.csv(file=paste(path,"xxxx/","BIO",var,"_Tree",1001,"_RateChange.csv",sep=""),row.names=1)
	resmean<-resmean[,c("ModMeanVal","ModMinVal","ModMaxVal","Age","AncVal","AbsRateMean","AbsRateMin","AbsRateMax")]
	colnames(resmean)<-c("ModMeanVal","ModMinVal","ModMaxVal","Age","AncVal","MCCAbsRateMean","MCCAbsRateMin","MCCAbsRateMax")
	resmean[,c("MeanAbsRate_Mean","MeanLower95Rate","MeanUpper95Rate","MinAbsRate_Mean","MinLower95Rate","MinUpper95Rate","MaxAbsRate_Mean","MaxLower95Rate","MaxUpper95Rate")]<-NA
	#there is a faster and more efficient way to do this, but use this for now
	for(otu in 1:nrow(resmean)){
		otuvals.mean<-NULL
		otuvals.min<-NULL
		otuvals.max<-NULL
		#get individual values for otu from each tree for the variable of interest
		for(t in 1:1000){
			trresmean<-NULL
			trresmean<-read.csv(file=paste(path,"/xxxx/","BIO",var,"_Tree",t,"_RateChange.csv",sep=""),row.names=1)
			otuvals.mean<-c(otuvals.mean,trresmean$AbsRateMean[rownames(trresmean) %in% rownames(resmean)[otu]])
			otuvals.min<-c(otuvals.min,trresmean$AbsRateMin[rownames(trresmean) %in% rownames(resmean)[otu]])
			otuvals.max<-c(otuvals.max,trresmean$AbsRateMax[rownames(trresmean) %in% rownames(resmean)[otu]])
		}
		#mean and 95% HPD using mean value
		resmean$MeanAbsRate_Mean[otu]<-mean(otuvals.mean)
		resmean$MeanLower95Rate[otu]<-HPDinterval(as.mcmc(otuvals.mean),0.95)[[1]]
		resmean$MeanUpper95Rate[otu]<-HPDinterval(as.mcmc(otuvals.mean),0.95)[[2]]
	
		#mean and 95% HPD using min value
		resmean$MinAbsRate_Mean[otu]<-mean(otuvals.min)
		resmean$MinLower95Rate[otu]<-HPDinterval(as.mcmc(otuvals.min),0.95)[[1]]
		resmean$MinUpper95Rate[otu]<-HPDinterval(as.mcmc(otuvals.min),0.95)[[2]]
	
		#mean and 95% HPD using max value
		resmean$MaxAbsRate_Mean[otu]<-mean(otuvals.max)
		resmean$MaxLower95Rate[otu]<-HPDinterval(as.mcmc(otuvals.max),0.95)[[1]]
		resmean$MaxUpper95Rate[otu]<-HPDinterval(as.mcmc(otuvals.max),0.95)[[2]]

	}
	write.csv(resmean,file=paste(path,"/xxxx/BIO",var,"_Summary_RateChange.csv",sep=""))
	print(head(resmean))
}	
	
path<-"/home/FM/mnelsen/trebouxia_phoebe/"

res1<-read.csv(file=paste(path,"/xxxx/BIO","1","_Summary_RateChange.csv",sep=""),row.names=1)
head(res1)

res5<-read.csv(file=paste(path,"/xxxx/BIO","5","_Summary_RateChange.csv",sep=""),row.names=1)
head(res5)

res6<-read.csv(file=paste(path,"/xxxx/BIO","6","_Summary_RateChange.csv",sep=""),row.names=1)
head(res6)

require(phytools)
mcc<-read.nexus(file="~/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre")
mcc$tip.label<-gsub("G","C",mcc$tip.label)
mcc$tip.label<-gsub("_97.5_","",mcc$tip.label)

plotTree(mcc,fsize=0.8)

res1<-read.csv(file="~/xxxx/BIO1_Summary_RateChange.csv",row.names=1,stringsAsFactors=FALSE)
rownames(res1)<-gsub("G","C",rownames(res1))
rownames(res1)<-gsub("_97.5_","",rownames(res1))
summary(res1[,c(6:17)])

res5<-read.csv(file="~/xxxx/BIO5_Summary_RateChange.csv",row.names=1,stringsAsFactors=FALSE)
rownames(res5)<-gsub("G","C",rownames(res5))
rownames(res5)<-gsub("_97.5_","",rownames(res5))
summary(res5[,c(6:17)])

res6<-read.csv(file="~/xxxx/BIO6_Summary_RateChange.csv",row.names=1,stringsAsFactors=FALSE)
rownames(res6)<-gsub("G","C",rownames(res6))
rownames(res6)<-gsub("_97.5_","",rownames(res6))
summary(res6[,c(6:17)])

vals<-cbind(res6[,"MCCAbsRateMean"],res5[,"MCCAbsRateMean"],res1[,"MCCAbsRateMean"])
rownames(vals)<-rownames(res1)

colnames(vals)<-c("BIO6","BIO5","BIO1")

pdf("~/xxxx/MCC_Mean_upd.pdf")
plotTree.barplot(mcc,vals,args.plotTree=list(fsize=0.4),args.barplot=list(col=c("royalblue","navajowhite4","olivedrab3"),border=NA,beside=TRUE,xlim=c(0,1.1),legend.text=TRUE,space=c(0,0.5),args.legend=list(x=1,y=length(mcc$tip.label)-50,bty="n")))
mtext("Absolute Rate of Change (\u00B0C/Myr)",1,at=0.5,line=2.5)
dev.off()

