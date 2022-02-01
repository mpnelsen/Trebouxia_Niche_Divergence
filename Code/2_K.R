#from lat_long_trebouxia21_MODIFIED_pmcUNIQUECLADESE_upd.R and plotting from K_genus_only_upd.R

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

#################
#K
#################
require(phytools)
#read sample of trees that have been renamed and pruned to only linclude those w bioclim data
trees<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees")

#read back in summary data w SE
taxa<-read.csv(file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES_wSE_corrected.csv",stringsAsFactors=FALSE,row.names=1)

#reduce to those that only have bioclim data
tbc<-taxa[!is.na(taxa$Mean_Bio1),]

kv<-as.data.frame(matrix(nrow=1000,ncol=39))
colnames(kv)<-c("Sample",paste("K","_","bio",1:19,sep=""),paste("P","_","bio",1:19,sep=""))

for(x in 1:nrow(kv)){
	for(b in 1:19){
		bio<-NULL
		bio<-tbc[,paste("Mean_Bio",b,sep="")]
		names(bio)<-rownames(tbc)
		
		sem<-NULL
		sem<-tbc[,paste("SE_Bio",b,sep="")]
		names(sem)<-rownames(tbc)

		k1<-NULL
		k1<-phylosig(trees[[x]],x=bio,se=sem,method='K',nsim=1000,test=TRUE)
		kv[x,paste("K","_","bio",b,sep="")]<-k1$K
		kv[x,paste("P","_","bio",b,sep="")]<-k1$P
	}
}

#add sample numbers
kv$Sample<-1:1000

#add in mcc tree
mcc<-kv[1,]
mcc[1,]<-NA
kvmcc<-rbind(kv,mcc)

kvmcc$Sample[1001]<-"MCC"

#read and clean mcc tree w outs and non-bioclim taxa dropped
tr<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre")

for(b in 1:19){
	bio<-NULL
	bio<-tbc[,paste("Mean_Bio",b,sep="")]
	names(bio)<-rownames(tbc)	

	sem<-NULL
	sem<-tbc[,paste("SE_Bio",b,sep="")]
	names(sem)<-rownames(tbc)

	k1<-NULL
	k1<-phylosig(tr,x=bio,se=sem,method='K',nsim=1000,test=TRUE)
	kvmcc[1001,paste("K","_","bio",b,sep="")]<-k1$K
	kvmcc[1001,paste("P","_","bio",b,sep="")]<-k1$P	
}

write.csv(kvmcc,file="/xxxx/K_mean_4nov2019_MODIFIED_UNIQUES_SE_corrected.csv",row.names=FALSE)


#make summary
vals<-read.csv(file="/xxxx/K_mean_4nov2019_MODIFIED_UNIQUES_SE_corrected.csv",stringsAsFactors=FALSE)


sums<-as.data.frame(matrix(nrow=19,ncol=8))
colnames(sums)<-c("K_MCC","P_MCC","Mean_K","Mean_P","Min_K","Min_P","Max_K","Max_P")
rownames(sums)<-paste("bio",1:19,sep="")

for(x in 1:19){
	sums[x,"K_MCC"]<-vals[1001,paste("K_","bio",x,sep="")]
	sums[x,"P_MCC"]<-vals[1001,paste("P_","bio",x,sep="")]
	sums[x,"Mean_K"]<-mean(vals[1:1000,paste("K_","bio",x,sep="")])
	sums[x,"Mean_P"]<-mean(vals[1:1000,paste("P_","bio",x,sep="")])
	sums[x,"Min_K"]<-min(vals[1:1000,paste("K_","bio",x,sep="")])
	sums[x,"Min_P"]<-min(vals[1:1000,paste("P_","bio",x,sep="")])
	sums[x,"Max_K"]<-max(vals[1:1000,paste("K_","bio",x,sep="")])
	sums[x,"Max_P"]<-max(vals[1:1000,paste("P_","bio",x,sep="")])
}

write.csv(sums,file="/Users/xxxx/summary_K_mean_4nov2019_MODIFIED_UNIQUES_SE_corrected.csv",row.names=TRUE)


#PLOT K and P values
require(ggplot2)
require(ggside)
vals<-read.csv(file="/xxxx/K_mean_4nov2019_MODIFIED_UNIQUES_SE_corrected.csv",stringsAsFactors=FALSE)
vals$Type<-NA
vals$Type[1:1000]<-"Posterior"
vals$Type[1001]<-"MCC"
vals$Type<-as.factor(vals$Type)


g1<-ggplot(vals[c(1:100),],aes(P_bio1,K_bio1,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio1[1001],vals$K_bio1[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO1 (Gr1, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g2<-ggplot(vals[1:1000,],aes(P_bio2,K_bio2,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio2[1001],vals$K_bio2[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO2 (Gr?, PC2/3)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g3<-ggplot(vals[1:1000,],aes(P_bio3,K_bio3,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio3[1001],vals$K_bio3[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO3 (Gr6, PC3)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g4<-ggplot(vals[1:1000,],aes(P_bio4,K_bio4,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio4[1001],vals$K_bio4[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO4 (Gr4, PC2)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g5<-ggplot(vals[1:1000,],aes(P_bio5,K_bio5,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio5[1001],vals$K_bio5[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO5 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g6<-ggplot(vals[1:1000,],aes(P_bio6,K_bio6,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio6[1001],vals$K_bio6[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO6 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g7<-ggplot(vals[1:1000,],aes(P_bio7,K_bio7,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio7[1001],vals$K_bio7[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO7 (Gr4, PC2)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g8<-ggplot(vals[1:1000,],aes(P_bio8,K_bio8,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio8[1001],vals$K_bio8[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO8 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g9<-ggplot(vals[1:1000,],aes(P_bio9,K_bio9,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio9[1001],vals$K_bio9[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO9 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g10<-ggplot(vals[1:1000,],aes(P_bio10,K_bio10,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio10[1001],vals$K_bio10[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO10 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g11<-ggplot(vals[1:1000,],aes(P_bio11,K_bio11,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio11[1001],vals$K_bio11[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO11 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g12<-ggplot(vals[1:1000,],aes(P_bio12,K_bio12,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio12[1001],vals$K_bio12[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO12 (Gr1, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g13<-ggplot(vals[1:1000,],aes(P_bio13,K_bio13,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio13[1001],vals$K_bio13[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO13 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g14<-ggplot(vals[1:1000,],aes(P_bio14,K_bio14,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio14[1001],vals$K_bio14[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO14 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g15<-ggplot(vals[1:1000,],aes(P_bio15,K_bio15,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio15[1001],vals$K_bio15[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO15 (Gr5, PC2)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g16<-ggplot(vals[1:1000,],aes(P_bio16,K_bio16,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio16[1001],vals$K_bio16[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO16 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g17<-ggplot(vals[1:1000,],aes(P_bio17,K_bio17,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio17[1001],vals$K_bio17[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO17 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g18<-ggplot(vals[1:1000,],aes(P_bio18,K_bio18,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio18[1001],vals$K_bio18[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO18 (Gr2, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)
g19<-ggplot(vals[1:1000,],aes(P_bio19,K_bio19,color=Type))+geom_point(shape=16,size=1.5,alpha=0.05,show.legend=FALSE)+annotate("point",vals$P_bio19[1001],vals$K_bio19[1001],shape=3,size=1)+scale_color_manual(values=c("#0072B2","black"))+geom_xsidedensity(aes(y=stat(density)),show.legend=FALSE)+geom_ysidedensity(aes(x=stat(density)),show.legend=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO19 (Gr3, PC1)",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+scale_ysidex_continuous(labels=NULL,breaks=NULL)+scale_xsidey_continuous(labels=NULL,breaks=NULL)

bp<-ggplot(vals,aes(P_bio1,K_bio1,color=Type,shape=Type))+geom_point(size=3,alpha=1,show.legend=TRUE)+annotate("point",vals$P_bio1[c(1001)],vals$K_bio1[1001],shape=3,size=3)+scale_color_manual(values=c("black","#0072B2"))+scale_shape_manual(values=c(3,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "grey"))+geom_vline(xintercept=0.05,linetype="dashed",color="grey",size=1)+labs(title="BIO1",x="p-value",y="K",hjust=0.5)+theme(plot.title = element_text(hjust = 0.5))+xlim(0,0.35)+ylim(0,2.0)+theme(legend.key=element_blank(),legend.title=element_blank())
#g_legend function from here: https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot and here http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(bp)
require(gridExtra)
require(grid)


pdf("/xxxx/K_plots_SE_corrected_genus_only_upd.pdf",width=8,height=10)
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,legend,nrow=5,top=textGrob("Phylogenetic Signal (K) in Climatic Variables", gp=gpar(fontsize=18,font=8,face="bold")))
dev.off()