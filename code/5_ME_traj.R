#plot module eigenprotein trajectory with diagnosis and traits

#rename MSSM diagnosis for plotting
targets.MSSM$Diagnosis=targets.MSSM$DxByPath
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==0] <- 'Control'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==1] <- 'AD.poss'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==2] <- 'AD.pro'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==3] <- 'AD.def'

#make ApoE4 as dosage
targets.BLSA$ApoE4 = targets.BLSA$ApoE
targets.BLSA$ApoE4[targets.BLSA$ApoE4==-1]=0
targets.Mayo$ApoE4=0
targets.Mayo$ApoE4[targets.Mayo$ApoE=='e24' | targets.Mayo$ApoE=='e34']=1
targets.Mayo$ApoE4[targets.Mayo$ApoE=='e44']=2
targets.Banner$ApoE4=0
targets.Banner$ApoE4[targets.Banner$ApoE=='e24' | targets.Banner$ApoE=='e34']=1
targets.Banner$ApoE4[targets.Banner$ApoE=='e44']=2

#order plotting so largest module will be plotted first
modorder=names(sort(table(geneInfo.cons$moduleColor.cons),decreasing=T))
memodorder=paste0('ME',modorder)
MEs_MSSM=MEs_MSSM[,memodorder]
MEs_BLSA=MEs_BLSA[,memodorder]
MEs_Banner=MEs_Banner[,memodorder]
MEs_ACT=MEs_ACT[,memodorder]
MEs_Mayo=MEs_Mayo[,paste0('ME',1:length(memodorder)-1)]
names(MEs_Mayo)=memodorder #change so Mayo column names are colors not numbers

#setup dataframe. These are the metadata of interest
corvars=c("MSSM.PlaqueMean",'MSSM.CDR','MSSM.NP1','BLSA.ApoE4','BLSA.CERAD','BLSA.BRAAK','ACT.CASI','ACT.CERAD','ACT.BRAAK','Banner.ApoE4','Banner.Cerad.NP','Banner.Braak.score','Mayo.ApoE4','Mayo.AgeAtDeath','Mayo.Gender')
corDF=data.frame(matrix(NA,nrow=ncol(MEs_Banner),ncol=length(corvars)))
names(corDF)=corvars
rownames(corDF)=colnames(MEs_Banner)

#function used within WGCNA to calculate correlation and p-value within bicorAndPvalue. Trajectory plot will match moduleTraitPvalue
WGCNAcorp=function(x,y){
    cor=signif(cor(x,y,use='pairwise'),2)
    corp = signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))),2)
    return(corp)
}

pdf(paste0(resdir,'ME_trajectory_Consensus.pdf'),width=16,height=12,useDingbats=F)

par(mfrow=c(5,4))
par(mar=c(5,6,4,2))

for (i in 1:ncol(MEs_Banner)) {
    #MSSM  
    boxplot(MEs_MSSM[,i]~factor(as.vector(as.factor(targets.MSSM$Diagnosis)),c('Control','AD.poss','AD.pro','AD.def')),col=substring(colnames(MEs_BLSA)[i],3,40),ylab="Eigenprotein Value",main=paste(colnames(MEs_BLSA)[i],'MSSM',sep='.'),xlab=NULL,las=2)
    verboseScatterplot(x=as.numeric(targets.MSSM[,"PlaqueMean"]),y=MEs_MSSM[,i],xlab="PlaqueMean",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19)
    corDF[i,'MSSM.PlaqueMean']=WGCNAcorp(as.numeric(targets.MSSM[,"PlaqueMean"]), MEs_MSSM[,i])
    verboseScatterplot(x=as.numeric(targets.MSSM[,"CDR"]),y=MEs_MSSM[,i],xlab="CDR",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19)
    corDF[i,'MSSM.CDR']=WGCNAcorp(as.numeric(targets.MSSM[,"CDR"]), MEs_MSSM[,i])
    verboseScatterplot(x=as.numeric(targets.MSSM[,"NP1"]),y=MEs_MSSM[,i],xlab="NP1",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19)
    corDF[i,'MSSM.NP1']=WGCNAcorp(as.numeric(targets.MSSM[,"NP1"]), MEs_MSSM[,i])
    #BLSA
    boxplot(MEs_BLSA[,i]~factor(as.vector(as.factor(targets.BLSA$Diagnosis)),c('Control','AsymAD','AD')),col=substring(colnames(MEs_BLSA)[i],3,40),ylab="Eigenprotein Value",main=paste(colnames(MEs_BLSA)[i],'BLSA',sep='.'),xlab=NULL,las=2)
    verboseScatterplot(x=targets.BLSA[,"ApoE4"],y=MEs_BLSA[,i],xlab="ApoE4 allele",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19);axis(side=1,at=0:2)
    corDF[i,'BLSA.ApoE4']=WGCNAcorp(targets.BLSA[,"ApoE4"], MEs_BLSA[,i])
    verboseScatterplot(x=as.numeric(factor(targets.BLSA[,"CERAD"])),y=MEs_BLSA[,i],xlab="CERAD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19)
    corDF[i,'BLSA.CERAD']=WGCNAcorp(as.numeric(factor(targets.BLSA[,"CERAD"])), MEs_BLSA[,i])
    verboseScatterplot(x=as.numeric(factor(targets.BLSA[,"BRAAK"])),y=MEs_BLSA[,i],xlab="BRAAK",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_BLSA)[i],3,40),pch=19)
    corDF[i,'BLSA.BRAAK']=WGCNAcorp(as.numeric(factor(targets.BLSA[,"BRAAK"])), MEs_BLSA[,i])
    #ACT
    boxplot(MEs_ACT[,i]~factor(as.vector(as.factor(targets.ACT$Diagnosis)),c('Control','AsymAD','AD')),col=substring(colnames(MEs_ACT)[i],3,40),ylab="Eigenprotein Value",main=paste(colnames(MEs_ACT)[i],'ACT',sep='.'),xlab=NULL,las=2)
    verboseScatterplot(x=as.numeric(factor(targets.ACT[,"CASI"])),y=MEs_ACT[,i],xlab="CASI",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_ACT)[i],3,40),pch=19)
    corDF[i,'ACT.CASI']=WGCNAcorp(as.numeric(factor(targets.ACT[,"CASI"])), MEs_ACT[,i])
    verboseScatterplot(x=as.numeric(factor(targets.ACT[,"CERAD"])),y=MEs_ACT[,i],xlab="CERAD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_ACT)[i],3,40),pch=19)
    corDF[i,'ACT.CERAD']=WGCNAcorp(as.numeric(factor(targets.ACT[,"CERAD"])), MEs_ACT[,i])
    verboseScatterplot(x=as.numeric(factor(targets.ACT[,"BRAAK"])),y=MEs_ACT[,i],xlab="BRAAK",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_ACT)[i],3,40),pch=19)
    corDF[i,'ACT.BRAAK']=WGCNAcorp(as.numeric(factor(targets.ACT[,"BRAAK"])), MEs_ACT[,i])
    #Banner
    boxplot(MEs_Banner[,i]~factor(as.vector(as.factor(targets.Banner$Diagnosis)),c('CTL','ASYM','MCI','AD')),col=substring(colnames(MEs_Banner)[i],3,40),ylab="Eigenprotein Value",main=paste(colnames(MEs_Banner)[i],'Banner',sep='.'),xlab=NULL,las=2)
    verboseScatterplot(x=targets.Banner[,"ApoE4"],y=MEs_Banner[,i],xlab="ApoE4 allele",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Banner)[i],3,40),pch=19);axis(side=1,at=0:2)
    corDF[i,'Banner.ApoE4']=WGCNAcorp(targets.Banner[,"ApoE4"], MEs_Banner[,i])
    verboseScatterplot(x=as.numeric(factor(targets.Banner[,"Cerad.NP"],c('Criteria not met','not AD','possible AD','probable AD','definite AD'))),y=MEs_Banner[,i],xlab="CERAD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Banner)[i],3,40),pch=19)
    corDF[i,'Banner.Cerad.NP']=WGCNAcorp(as.numeric(factor(targets.Banner[,"Cerad.NP"],,c('Criteria not met','not AD','possible AD','probable AD','definite AD'))), MEs_Banner[,i])
    verboseScatterplot(x=as.numeric(factor(targets.Banner[,"Braak.score"],c('I','II','III','IV','V','VI'))),y=MEs_Banner[,i],xlab="BRAAK",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Banner)[i],3,40),pch=19)
    corDF[i,'Banner.Braak.score']=WGCNAcorp(as.numeric(factor(targets.Banner[,"Braak.score"],c('I','II','III','IV','V','VI'))), MEs_Banner[,i])
    #May
    boxplot(MEs_Mayo[,i]~factor(as.vector(as.factor(targets.Mayo$Diagnosis)),c('Control','AD','PSP')),col=substring(colnames(MEs_Mayo)[i],3,40),ylab="Eigenprotein Value",main=paste(colnames(MEs_Mayo)[i],'Mayo',sep='.'),xlab=NULL,las=2)
    verboseScatterplot(x=targets.Mayo[,"ApoE4"],y=MEs_Mayo[,i],xlab="ApoE4 allele",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Mayo)[i],3,40),pch=19);axis(side=1,at=0:2)
    corDF[i,'Mayo.ApoE4']=WGCNAcorp(targets.Mayo[,"ApoE4"], MEs_Mayo[,i])
    verboseScatterplot(x=as.numeric(factor(targets.Mayo[,"AgeAtDeath"])),y=MEs_Mayo[,i],xlab="AgeAtDeath",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Mayo)[i],3,40),pch=19)
    corDF[i,'Mayo.AgeAtDeath']=WGCNAcorp(as.numeric(factor(targets.Mayo[,"AgeAtDeath"])), MEs_Mayo[,i])
    verboseScatterplot(x=as.numeric(factor(targets.Mayo[,"Gender"])),y=MEs_Mayo[,i],xlab="Gender",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=substring(colnames(MEs_Mayo)[i],3,40),pch=19)
    corDF[i,'Mayo.Gender']=WGCNAcorp(as.numeric(factor(targets.Mayo[,"Gender"])), MEs_Mayo[,i])

}

dev.off()

write.csv(corDF,file=paste0(resdir,'ME_trajectory_Consensus.csv'))

