#from rate_values4_bestfitmodel.R

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

require(coda)
require(ggpubr)
require(dplyr)
path<-"xxxx/"

toplot5<-read.csv(file=paste(path,"bio",5,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
sigsq5<-data.frame(matrix(nrow=nrow(toplot5),ncol=2))
colnames(sigsq5)<-c("Tree","Sigsq")
sigsq5$Tree<-toplot5[,"Tree"]
for(x in 1:nrow(sigsq5)){
	sigsq5$Sigsq[x]<-toplot5[x,paste(toplot5$Best[x],"_sigsq",sep="")]
}
sigsq5$Var<-"BIO5"
sigsq5$Dist<-NA
sigsq5$Dist[1:1000]<-"Posterior"
sigsq5$Dist[1001]<-"MCC"

toplot6<-read.csv(file=paste(path,"bio",6,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
sigsq6<-data.frame(matrix(nrow=nrow(toplot6),ncol=2))
colnames(sigsq6)<-c("Tree","Sigsq")
sigsq6$Tree<-toplot6[,"Tree"]
for(x in 1:nrow(sigsq6)){
	sigsq6$Sigsq[x]<-toplot6[x,paste(toplot6$Best[x],"_sigsq",sep="")]
}
sigsq6$Var<-"BIO6"
sigsq6$Dist<-NA
sigsq6$Dist[1:1000]<-"Posterior"
sigsq6$Dist[1001]<-"MCC"

HPDinterval(as.mcmc(sigsq5$Sigsq[1:1000]),0.95)
HPDinterval(as.mcmc(sigsq6$Sigsq[1:1000]),0.95)

sigsq56<-rbind(sigsq5,sigsq6)
sigsq56$Var<-as.factor(sigsq56$Var)
sigsq56$Dist<-as.factor(sigsq56$Dist)

p<-ggplot(sigsq56 %>% filter(Dist=="Posterior"),aes(x=Var,y=Sigsq,fill=Var))+geom_violin(trim=FALSE,color=NA)+scale_fill_manual(values=c("navajowhite4","royalblue"))+theme_minimal()+theme(legend.title=element_blank(),legend.position="none")+ggtitle("Thermal Variables")+ylab(bquote('Evolutionary Rate '(sigma^2)))+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank())+geom_segment(aes(x=1,xend=1,y=HPDinterval(as.mcmc(sigsq5$Sigsq[1:1000]),0.95)[[1]],yend=HPDinterval(as.mcmc(sigsq5$Sigsq[1:1000]),0.95)[[2]]),colour="gray77",size=0.25)+geom_segment(aes(x=2,xend=2,y=HPDinterval(as.mcmc(sigsq6$Sigsq[1:1000]),0.95)[[1]],yend=HPDinterval(as.mcmc(sigsq6$Sigsq[1:1000]),0.95)[[2]]),colour="gray77",size=0.25)+geom_point(data=sigsq56 %>% filter(Dist=="MCC"),aes(x=Var,y=Sigsq),size=5,color="black")
p

toplot16<-read.csv(file=paste(path,"bio",16,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
sigsq16<-data.frame(matrix(nrow=nrow(toplot16),ncol=2))
colnames(sigsq16)<-c("Tree","Sigsq")
sigsq16$Tree<-toplot16[,"Tree"]
for(x in 1:nrow(sigsq16)){
	sigsq16$Sigsq[x]<-toplot16[x,paste(toplot16$Best[x],"_sigsq",sep="")]
}
sigsq16$Var<-"BIO16"
sigsq16$Dist<-NA
sigsq16$Dist[1:1000]<-"Posterior"
sigsq16$Dist[1001]<-"MCC"

toplot17<-read.csv(file=paste(path,"bio",17,"_modelfitting_MODIFIED_UNIQUES_wSE_corrected.csv",sep=""),stringsAsFactors=FALSE)
sigsq17<-data.frame(matrix(nrow=nrow(toplot17),ncol=2))
colnames(sigsq17)<-c("Tree","Sigsq")
sigsq17$Tree<-toplot17[,"Tree"]
for(x in 1:nrow(sigsq17)){
	sigsq17$Sigsq[x]<-toplot17[x,paste(toplot17$Best[x],"_sigsq",sep="")]
}
sigsq17$Var<-"BIO17"
sigsq17$Dist<-NA
sigsq17$Dist[1:1000]<-"Posterior"
sigsq17$Dist[1001]<-"MCC"

HPDinterval(as.mcmc(sigsq16$Sigsq[1:1000]),0.95)
HPDinterval(as.mcmc(sigsq17$Sigsq[1:1000]),0.95)

sigsq1617<-rbind(sigsq16,sigsq17)
sigsq1617$Var<-as.factor(sigsq1617$Var)
sigsq1617$Dist<-as.factor(sigsq1617$Dist)

p2<-ggplot(sigsq1617 %>% filter(Dist=="Posterior"),aes(x=Var,y=Sigsq,fill=Var))+geom_violin(trim=FALSE,color=NA)+scale_fill_manual(values=c("lightblue4","coral4"))+theme_minimal()+theme(legend.title=element_blank(),legend.position="none")+ggtitle("Moisture Variables")+ylab(bquote('Evolutionary Rate '(sigma^2)))+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank())+geom_segment(aes(x=1,xend=1,y=HPDinterval(as.mcmc(sigsq16$Sigsq[1:1000]),0.95)[[1]],yend=HPDinterval(as.mcmc(sigsq16$Sigsq[1:1000]),0.95)[[2]]),colour="gray77",size=0.25)+geom_segment(aes(x=2,xend=2,y=HPDinterval(as.mcmc(sigsq17$Sigsq[1:1000]),0.95)[[1]],yend=HPDinterval(as.mcmc(sigsq17$Sigsq[1:1000]),0.95)[[2]]),colour="gray77",size=0.25)+geom_point(data=sigsq1617 %>% filter(Dist=="MCC"),aes(x=Var,y=Sigsq),size=5,color="black")
p2

p3<-ggarrange(p,p2,ncol=2,nrow=1)
p3

pdf("xxxx/rates_of_evolution_upd.pdf")
p3
dev.off()