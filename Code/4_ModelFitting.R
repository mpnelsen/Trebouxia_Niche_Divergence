#from lat_long_trebouxia20_MODIFIED.r and lat_long_trebouxia21_MODIFIED_pmcUNIQUECLADESE_upd.R

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

##########################
##########################
#############FIT MODELS of trait evolution
##########################

#fit models
require(geiger)

#read sample of trees that have been renamed and pruned to only linclude those w bioclim data
#trees<-read.nexus(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_symbiotic_algae_review/24june2018/working_dir_28oct2019/individual_97.5_clusters/repseq_alignments_by_clade/beast/outs_w_trebouxia_merged_profiled_10pcntburn_removed_1krandomlyselected_renamed_andretainedonlywbioclim.trees")
trees<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees")

#read back in summary data
#taxa<-read.csv(file="~/Desktop/stuff_otu_summary.csv",stringsAsFactors=FALSE,row.names=1)
taxa<-read.csv(file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES_wSE_corrected.csv",stringsAsFactors=FALSE,row.names=1)

#reduce to those that only have bioclim data
tbc<-taxa[!is.na(taxa$Mean_Bio1),]

#read and clean mcc tree
#tr<-read.nexus(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_symbiotic_algae_review/24june2018/working_dir_28oct2019/individual_97.5_clusters/repseq_alignments_by_clade/beast/mcc_medhts_outs_w_trebouxia_merged_profiled_10pcntburn_removed.tre")
tr<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre")

#add mcc tree to end of trees (make tree 1001)
trees<-c(trees,tr)

fitmodels<-function(trees,biocno,dat){
	#make empty data frame
	cols<-c("Tree","BM_lnL","BM_AICc","BM_sigsq","BM_z0","BM_rows","BM_convergence","BM_Delta_AICc","BM_AICc_Weight","OU_lnL","OU_AICc","OU_alpha","OU_sigsq","OU_z0","OU_rows","OU_convergence","OU_Delta_AICc","OU_AICc_Weight","EB_lnL","EB_AICc","EB_a","EB_sigsq","EB_z0","EB_rows","EB_convergence","EB_Delta_AICc","EB_AICc_Weight","Best","Best_AICc_Weight")
	df<-as.data.frame(matrix(nrow=1001,ncol=length(cols)))
	colnames(df)<-cols
	df$Tree[1:1000]<-1:1000
	df$Tree[1001]<-"MCC"
	
	#then go through each tree
	for(x in 1:1001){
		bio<-NULL
		bio<-dat[,paste("Mean_Bio",biocno,sep="")]
		names(bio)<-rownames(dat)	

		sem<-NULL
		sem<-dat[,paste("SE_Bio",biocno,sep="")]
		names(sem)<-rownames(dat)
	
		
		#BM
		Fit<-NULL
		Fit<-fitContinuous(trees[[x]],dat=bio,SE=sem,model="BM",ncores=2)
		df[x,c("BM_lnL","BM_AICc","BM_sigsq","BM_z0","BM_rows","BM_convergence")]<-c(Fit$opt$lnL,Fit$opt$aicc,Fit$opt$sigsq,Fit$opt$z0,nrow(Fit$res),sum(Fit$res[,"convergence"]))

		#OU
		Fit<-NULL
		Fit<-fitContinuous(trees[[x]],dat=bio,SE=sem,model="OU",ncores=2)
		df[x,c("OU_lnL","OU_AICc","OU_alpha","OU_sigsq","OU_z0","OU_rows","OU_convergence")]<-c(Fit$opt$lnL,Fit$opt$aicc,Fit$opt$alpha,Fit$opt$sigsq,Fit$opt$z0,nrow(Fit$res),sum(Fit$res[,"convergence"]))

		#EB
		Fit<-NULL
		Fit<-fitContinuous(trees[[x]],dat=bio,SE=sem,model="EB",ncores=2)
		df[x,c("EB_lnL","EB_AICc","EB_a","EB_sigsq","EB_z0","EB_rows","EB_convergence")]<-c(Fit$opt$lnL,Fit$opt$aicc,Fit$opt$a,Fit$opt$sigsq,Fit$opt$z0,nrow(Fit$res),sum(Fit$res[,"convergence"]))

		#calculate aicc weights
		weights<-aicw(as.numeric(df[x,c("BM_AICc","OU_AICc","EB_AICc")]))
		df[x,c("BM_Delta_AICc","OU_Delta_AICc","EB_Delta_AICc")]<-weights[,"delta"]
		df[x,c("BM_AICc_Weight","OU_AICc_Weight","EB_AICc_Weight")]<-weights[,"w"]
		models<-c("BM","OU","EB")
		df[x,"Best"]<-models[weights$w==max(weights$w)]
		df[x,"Best_AICc_Weight"]<-max(weights$w)
	}
	write.csv(df,file=paste("/xxxx/bio",biocno,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),row.names=FALSE)
}

fitmodels(trees=trees,biocno=19,dat=tbc)


###################
#order according to ppca
#MCC
path<-"/xxxx/"
for(x in 1:19){
	toplot<-NULL
	toplot.mcc<-NULL
	toplot<-read.csv(file=paste(path,"bio",x,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
	toplot.mcc<-toplot[1001,c("BM_AICc_Weight","OU_AICc_Weight","EB_AICc_Weight")]
	rownames(toplot.mcc)<-"MCC"
	toplot.mcc<-t(toplot.mcc)
	toplot.mcc<-as.matrix(toplot.mcc[nrow(toplot.mcc):1,])
	if(x==1){
		plotd.mcc<-toplot.mcc
	}
	if(x>1){
		plotd.mcc<-cbind(plotd.mcc,toplot.mcc)
	}
	if(x==19){
		colnames(plotd.mcc)<-paste("BIO",1:19,sep="")
		rownames(plotd.mcc)<-c("Early Burst (EB)","Ornstein-Uhlenbeck (OU)","Brownian Motion (BM)")
	}
}


#Posterior
path<-"/xxxx/"
for(x in 1:19){
	toplot<-NULL
	toplot.mcc<-NULL
	toplot<-read.csv(file=paste(path,"bio",x,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
	toplot.post<-toplot[1:1000,c("BM_AICc_Weight","OU_AICc_Weight","EB_AICc_Weight")]
	tplpmn<-colMeans(toplot.post)
	rn<-names(tplpmn)
	tplpmn<-as.data.frame(matrix(tplpmn,nrow=1))
	rownames(tplpmn)<-"Posterior"
	colnames(tplpmn)<-rn
	tplpmn<-t(tplpmn)
	tplpmn<-as.matrix(tplpmn[nrow(tplpmn):1,])
	if(x==1){
		plotdp<-tplpmn
	}
	if(x>1){
		plotdp<-cbind(plotdp,tplpmn)
	}
	if(x==19){
		colnames(plotdp)<-paste("BIO",1:19,sep="")
		rownames(plotdp)<-c("Early Burst (EB)","Ornstein-Uhlenbeck (OU)","Brownian Motion (BM)")
	}
}

#rename and reorder
plotd.mcc<-plotd.mcc[,c(1,12,5,8,10,13,16,18,6,9,11,14,17,19,4,7,15,2,3)]
plotdp<-plotdp[,c(1,12,5,8,10,13,16,18,6,9,11,14,17,19,4,7,15,2,3)]

#plot
my.cbPalette <- c("darkgoldenrod3","yellowgreen","plum4")
pdf(file="~/model_fits_summary_better_cols_wSE_corrected_reordered_upd2.pdf",width=9,height=6)
par(mfrow=c(2,1))
barplot(plotd.mcc, main="AICc Weights for Evolutionary Models (MCC Tree)",xlab=NULL,ylab="AICc Weight",col=my.cbPalette,xlim=c(0,34),las=2,border=NA)
legend(x="topright",legend=rev(rownames(plotd.mcc)),col=rev(my.cbPalette),pch=16,pt.cex=2,box.lty=0,horiz=F)
barplot(plotdp, main="Average AICc Weights for Evolutionary Models (Posterior Distribution of Trees)",xlab=NULL,ylab="Average AICc Weight",col=my.cbPalette,xlim=c(0,34),las=2,border=NA)
legend(x="topright",legend=rev(rownames(plotd.mcc)),col=rev(my.cbPalette),pch=16,pt.cex=2,box.lty=0,horiz=F)
dev.off()
