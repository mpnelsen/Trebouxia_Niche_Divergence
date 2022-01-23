#from plot_multiple_dtts_temp_precip_facet_new._28aug2020.r

#The following data files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/csv_files_etc/USED_FOR_ANALYSES
#MASTER_accession_list_INFOADDED_MODIFIED.csv
#OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv

#The following tree files can be obtained from https://github.com/mpnelsen/Trebouxia_Macroecology/tree/main/analysis/files/algal_trees
#burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees
#burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED.tre

######################
#DTT
######################
mult.dtts<-function(data.frame=NULL,colnames=NULL,varnames=NULL,tree=NULL,nsim=NULL,path.to.save=NULL){
	for(p in 1:length(colnames)){
		var1<-data.frame[,colnames[p]]
		names(var1)<-rownames(data.frame)
		dtt.var1<-dtt(phy=tree,data=var1,index="avg.sq",plot=FALSE,calculateMDIp=TRUE,nsim=nsim)
		save(dtt.var1,file=paste(path.to.save,"dtt_",colnames[p],"_29aug2020.Rsave",sep=""))
		#ltt <- sort(branching.times(tree), decreasing = TRUE)
		#ltt <- c(0, (max(ltt) - ltt)/max(ltt))
		ylim = c(range(pretty(dtt.var1$dtt)))
		plot(dtt.var1$times, dtt.var1$dtt, xlab = "Time", ylab = "Mean Subclade Disparity", ylim = ylim, bty = "n", type = "n",axes=FALSE,main=paste("Mean Subclade Disparity Through Time (",coldefs[p],")",sep=""))
		axis(2,labels=TRUE)
		height<-round(max(nodeHeights(tree)))
		axis(1,at=c(0,1-125/height,1-100/height,1-75/height,1-50/height,1-25/height,1),labels=c(height,125,100,75,50,25,0))
		#for x in 1:length things in timescale, make a polygon...
		for(x in 1:nrow(timey)){
			rect(xleft=1-timey$Start[x]/height,ybottom=-0.2,xright=1-timey$End[x]/height,ytop=ylim[2]+0.25,col=timey$RGB[x],border=NA)
		}
		for(x in 1:nrow(timey)){
			text(x=1-timey$Midpoint[x]/height,y=-0.022*ylim[[2]],cex=0.75,col="grey35",labels=timey$Abbrev[x])
		}
		lines(x=c(-0.1,1),y=c(0,0))
		poly = geiger:::.dtt.polygon(dtt.var1$sim, dtt.var1$times, alpha = 1 - 0.95)
		polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgrey", 0.5), border = NA)
		median.sims <- apply(dtt.var1$sim, 1, median)
		lines(dtt.var1$times, median.sims, lty = 2)
		lines(dtt.var1$times, dtt.var1$dtt, type = "l", lwd = 2)
		text(x=0.8,y=ylim[2],labels=paste("MDI = ", round(dtt.var1$MDI[[1]],3),sep=""))
		if(dtt.var1$MDIpVal[[1]]<0.001){
			text(x=0.8,y=ylim[2]-(ylim[2]*0.05),labels=paste("p < ", round(1/ncol(dtt.var1$sim),3), sep=""))
		}
		if(dtt.var1$MDIpVal[[1]]>=0.001){
			text(x=0.8,y=ylim[2]-(ylim[2]*0.05),labels=paste("p = ", round(dtt.var1$MDIpVal[[1]],3),sep=""))
		}
	}	
}


#TEMP IN DEGREES

require(geiger)
require(phytools)
timescale<-read.csv(file="/xxxx/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)
for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<150,]
timey<-timey[timey$Start>0,]
timey$End[1]<-0
timey$Start[nrow(timey)]<-150
timey$Midpoint[nrow(timey)]<-timey$End[nrow(timey)]+(timey$Start[nrow(timey)]-timey$End[nrow(timey)])/2
tree<-read.tree(file="/xxxx/burnin25pcn_resampevery20k_mcc_medhts_OTUs_renamed_bioclim_retained_MODIFIED_ladderized.tre")
dat.in.tre<-read.csv(file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv",stringsAsFactors=FALSE,row.names=1)
dat.in.tre<-dat.in.tre[order(match(rownames(dat.in.tre),tree$tip.label)),]
tree
nrow(dat.in.tre)

###JUST DO ALL
cols<-c("Mean_Bio1","Mean_Bio2","Mean_Bio3","Mean_Bio4","Mean_Bio5","Mean_Bio6")
coldefs<-c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6")
path.to.save<-"/xxxx/"
png("/xxxx/dtt_mults_1-6.png",width=11,height=11,units="in",res=2400)
par(mfrow=c(3,2))
mult.dtts(data.frame=dat.in.tre,colnames=cols,varnames=coldefs,tree=tree,nsim=1000,path.to.save=path.to.save)
dev.off()

cols<-c("Mean_Bio7","Mean_Bio8","Mean_Bio9","Mean_Bio10","Mean_Bio11","Mean_Bio12")
coldefs<-c("BIO7","BIO8","BIO9","BIO10","BIO11","BIO12")
path.to.save<-"/xxxx/"
png("/xxxx/dtt_mults_7-12.png",width=11,height=11,units="in",res=2400)
par(mfrow=c(3,2))
mult.dtts(data.frame=dat.in.tre,colnames=cols,varnames=coldefs,tree=tree,nsim=1000,path.to.save=path.to.save)
dev.off()

cols<-c("Mean_Bio13","Mean_Bio14","Mean_Bio15","Mean_Bio16","Mean_Bio17","Mean_Bio18","Mean_Bio19")
coldefs<-c("BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
path.to.save<-"/xxxx/"
png("/xxxx/dtt_mults_13-19.png",width=11,height=11,units="in",res=2400)
par(mfrow=c(4,2))
mult.dtts(data.frame=dat.in.tre,colnames=cols,varnames=coldefs,tree=tree,nsim=1000,path.to.save=path.to.save)
dev.off()




###########
###########
###########
#DTT on multiple traits and over distribution of trees
###########
###########

mult.dtts<-function(data.frame=NULL,colnames=NULL,trees=NULL,nsim=NULL,path.to.save=NULL,df=NULL){
	for(t in 1:length(trees)){
		df$Tree[t]<-t
		for(p in 1:length(colnames)){
			var1<-data.frame[,colnames[p]]
			names(var1)<-rownames(data.frame)
			dtt.var1<-dtt(phy=tree[[t]],data=var1,index="avg.sq",plot=FALSE,calculateMDIp=TRUE,nsim=nsim)
			save(dtt.var1,file=paste(path.to.save,"dtt_tree_",t,"_",colnames[p],"_new.Rsave",sep=""))
			df[t,(p*2)]<-round(dtt.var1$MDI[[1]],3)
			df[t,((p*2)+1)]<-dtt.var1$MDIpVal[[1]]
		}
	}	
return(df)
}

require(geiger)
tree<-read.nexus(file="/xxxx/burnin25pcn_resampevery20k_OTUs_renamed_bioclim_retained_1krandomlyselected_MODIFIED.trees")
dat.in.tre<-read.csv(file="/xxxx/OTU_97.5_BIOCLIM_SUMMARY_MODIFIED_UNIQUES.csv",stringsAsFactors=FALSE,row.names=1)
#dat.in.tre<-dat.in.tre[order(match(rownames(dat.in.tre),tree$tip.label)),]
#tree
nrow(dat.in.tre)

sc<-c("Tree","MDI_1","p_1","MDI_2","p_2","MDI_3","p_3","MDI_4","p_4","MDI_5","p_5","MDI_6","p_6","MDI_7","p_7","MDI_8","p_8","MDI_9","p_9","MDI_10","p_10","MDI_11","p_11","MDI_12","p_12","MDI_13","p_13","MDI_14","p_14","MDI_15","p_15","MDI_16","p_16","MDI_17","p_17","MDI_18","p_18","MDI_19","p_19")
df<-as.data.frame(matrix(nrow=1000,ncol=length(sc)))
colnames(df)<-sc

##Regular of all variables
cols<-c("Mean_Bio1","Mean_Bio2","Mean_Bio3","Mean_Bio4","Mean_Bio5","Mean_Bio6","Mean_Bio7","Mean_Bio8","Mean_Bio9","Mean_Bio10","Mean_Bio11","Mean_Bio12","Mean_Bio13","Mean_Bio14","Mean_Bio15","Mean_Bio16","Mean_Bio17","Mean_Bio18","Mean_Bio19")
path.to.save<-"/xxxx/"
ndf<-mult.dtts(data.frame=dat.in.tre,colnames=cols,tree=tree,nsim=1000,path.to.save=path.to.save,df=df)
write.csv(ndf,file="/xxxx/dtt_on_pp_28aug2020.csv",row.names=FALSE)

#make value names and empty data frame
mdi<-c(paste("MDI_",1:19,sep=""))
pv<-c(paste("p_",1:19,sep=""))
cols<-c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
sc<-c("Bio","MDI_Min","MDI_Mean","MDI_Max","P_Min","P_Mean","P_Max")
df<-as.data.frame(matrix(nrow=length(cols),ncol=length(sc)))
colnames(df)<-sc
df$Bio<-cols

for(x in 1:19){
	#get MDI values
	df[x,"MDI_Min"]<-round(min(ndf[,mdi[x]]),4)
	df[x,"MDI_Mean"]<-round(mean(ndf[,mdi[x]]),4)
	df[x,"MDI_Max"]<-round(max(ndf[,mdi[x]]),4)
	
	#get p-values
	df[x,"P_Min"]<-round(min(ndf[,pv[x]]),4)
	df[x,"P_Mean"]<-round(mean(ndf[,pv[x]]),4)
	df[x,"P_Max"]<-round(max(ndf[,pv[x]]),4)
}

write.csv(df,file="/xxxx/dtt_on_pp_28aug2020_summary.csv",row.names=FALSE)

round(min(ndf$MDI_1),4)
round(min(ndf$MDI_2),4)
round(min(ndf$MDI_3),4)
round(min(ndf$MDI_4),4)
round(min(ndf$MDI_5),4)
round(min(ndf$MDI_6),4)
round(min(ndf$MDI_7),4)
round(min(ndf$MDI_8),4)
round(min(ndf$MDI_9),4)
round(min(ndf$MDI_10),4)
round(min(ndf$MDI_11),4)
round(min(ndf$MDI_12),4)
round(min(ndf$MDI_13),4)
round(min(ndf$MDI_14),4)
round(min(ndf$MDI_15),4)
round(min(ndf$MDI_16),4)
round(min(ndf$MDI_17),4)
round(min(ndf$MDI_18),4)
round(min(ndf$MDI_19),4)

round(mean(ndf$MDI_1),4)
round(mean(ndf$MDI_2),4)
round(mean(ndf$MDI_3),4)
round(mean(ndf$MDI_4),4)
round(mean(ndf$MDI_5),4)
round(mean(ndf$MDI_6),4)
round(mean(ndf$MDI_7),4)
round(mean(ndf$MDI_8),4)
round(mean(ndf$MDI_9),4)
round(mean(ndf$MDI_10),4)
round(mean(ndf$MDI_11),4)
round(mean(ndf$MDI_12),4)
round(mean(ndf$MDI_13),4)
round(mean(ndf$MDI_14),4)
round(mean(ndf$MDI_15),4)
round(mean(ndf$MDI_16),4)
round(mean(ndf$MDI_17),4)
round(mean(ndf$MDI_18),4)
round(mean(ndf$MDI_19),4)

round(max(ndf$MDI_1),4)
round(max(ndf$MDI_2),4)
round(max(ndf$MDI_3),4)
round(max(ndf$MDI_4),4)
round(max(ndf$MDI_5),4)
round(max(ndf$MDI_6),4)
round(max(ndf$MDI_7),4)
round(max(ndf$MDI_8),4)
round(max(ndf$MDI_9),4)
round(max(ndf$MDI_10),4)
round(max(ndf$MDI_11),4)
round(max(ndf$MDI_12),4)
round(max(ndf$MDI_13),4)
round(max(ndf$MDI_14),4)
round(max(ndf$MDI_15),4)
round(max(ndf$MDI_16),4)
round(max(ndf$MDI_17),4)
round(max(ndf$MDI_18),4)
round(max(ndf$MDI_19),4)

round(min(ndf$p_1),4)
round(min(ndf$p_2),4)
round(min(ndf$p_3),4)
round(min(ndf$p_4),4)
round(min(ndf$p_5),4)
round(min(ndf$p_6),4)
round(min(ndf$p_7),4)
round(min(ndf$p_8),4)
round(min(ndf$p_9),4)
round(min(ndf$p_10),4)
round(min(ndf$p_11),4)
round(min(ndf$p_12),4)
round(min(ndf$p_13),4)
round(min(ndf$p_14),4)
round(min(ndf$p_15),4)
round(min(ndf$p_16),4)
round(min(ndf$p_17),4)
round(min(ndf$p_18),4)
round(min(ndf$p_19),4)

round(mean(ndf$p_1),4)
round(mean(ndf$p_2),4)
round(mean(ndf$p_3),4)
round(mean(ndf$p_4),4)
round(mean(ndf$p_5),4)
round(mean(ndf$p_6),4)
round(mean(ndf$p_7),4)
round(mean(ndf$p_8),4)
round(mean(ndf$p_9),4)
round(mean(ndf$p_10),4)
round(mean(ndf$p_11),4)
round(mean(ndf$p_12),4)
round(mean(ndf$p_13),4)
round(mean(ndf$p_14),4)
round(mean(ndf$p_15),4)
round(mean(ndf$p_16),4)
round(mean(ndf$p_17),4)
round(mean(ndf$p_18),4)
round(mean(ndf$p_19),4)

round(max(ndf$p_1),4)
round(max(ndf$p_2),4)
round(max(ndf$p_3),4)
round(max(ndf$p_4),4)
round(max(ndf$p_5),4)
round(max(ndf$p_6),4)
round(max(ndf$p_7),4)
round(max(ndf$p_8),4)
round(max(ndf$p_9),4)
round(max(ndf$p_10),4)
round(max(ndf$p_11),4)
round(max(ndf$p_12),4)
round(max(ndf$p_13),4)
round(max(ndf$p_14),4)
round(max(ndf$p_15),4)
round(max(ndf$p_16),4)
round(max(ndf$p_17),4)
round(max(ndf$p_18),4)
round(max(ndf$p_19),4)


