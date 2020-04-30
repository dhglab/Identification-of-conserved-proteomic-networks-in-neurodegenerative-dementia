###correlate module eigenprotein with trait

#save module eigenprotein to ned variable
MEs_MSSM=merged$newMEs[[1]]$data
MEs_BLSA=merged$newMEs[[2]]$data
MEs_ACT=merged$newMEs[[3]]$data
MEs_Banner=merged$newMEs[[4]]$data
MEs_Mayo=merged$newMEs[[5]]$data

#label modules numbers as colors
tmp=as.data.frame(cbind(labels2colors(c(0,1,2,3,4,5,6,7,8,9,10,11,12)),paste('ME',c(0,1,2,3,4,5,6,7,8,9,10,11,12),sep='')))
tmp=tmp[match(colnames(MEs_MSSM),tmp$V2),]
colnames(MEs_MSSM)<- colnames(MEs_BLSA) <- colnames(MEs_ACT) <- colnames(MEs_Banner) <- paste('ME',tmp$V1,sep='')

#create KME
consKME1=consensusKME(multiExpr=multiExpr, moduleColor.cons,
    multiEigengenes = NULL,
    consensusQuantile = 0.2,
    signed = TRUE)

consensus.KMEs=consKME1[,regexpr('consensus.kME',names(consKME1))>0]
geneInfo.cons=as.data.frame(cbind(colnames(multiExpr[[1]]$data),colnames(multiExpr[[1]]$data),moduleColor.cons,consensus.KMEs))
colnames(geneInfo.cons)[c(1,2)] <- c('ID','Symbol')
geneInfo.cons$Symbol=do.call('rbind',strsplit(as.character(geneInfo.cons$Symbol),'[|]'))[,1]
#write out module membership and KME to file
write.csv(geneInfo.cons,paste0(resdir,'geneInfo.cons.csv'))


##set up metadata of cohorts and perform ME correlation

#Mayo#
numericMeta=targets.Mayo[,c(3,7,4,5,2)]
numericMeta$Diagnosis=as.numeric(factor(numericMeta$Diagnosis,c('Control','AD','PSP')))-1
numericMeta$Gender=as.numeric(factor(numericMeta$Gender))-1
numericMeta$ApoE=as.numeric(factor(numericMeta$ApoE))-1
numericMeta$batch=as.numeric(factor(numericMeta$batch))-1
#PSP vs Control
numericMeta.PSP=subset(numericMeta,Diagnosis!='1') ## 0=con, 1=ad, 2=psp
MEs.PSP=MEs_Mayo[targets.Mayo$Diagnosis!='AD',]
MEcors.PSP <- bicorAndPvalue(MEs.PSP,numericMeta.PSP)
moduleTraitCor.PSP <- MEcors.PSP$bicor
moduleTraitPvalue.PSP <- MEcors.PSP$p
colnames(moduleTraitPvalue.PSP) <- colnames(moduleTraitCor.PSP) <- paste('PSPvsControl',colnames(moduleTraitCor.PSP),sep='\n')
#AD vs Control
numericMeta.AD=subset(numericMeta,Diagnosis!='2')
MEs.AD=MEs_Mayo[targets.Mayo$Diagnosis!='PSP',]
MEcors.AD <- bicorAndPvalue(MEs.AD,numericMeta.AD)
moduleTraitCor.AD <- MEcors.AD$bicor
moduleTraitPvalue.AD <- MEcors.AD$p
colnames(moduleTraitPvalue.AD) <- colnames(moduleTraitCor.AD) <- paste('ADvsControl',colnames(moduleTraitCor.AD),sep='\n')
#PSP vs AD
numericMeta.AD.PSP=subset(numericMeta,Diagnosis!='0') ## no control
numericMeta.AD.PSP$Diagnosis=numericMeta.AD.PSP$Diagnosis-1 ## AD=0, PSP=1
MEs.AD.PSP=MEs_Mayo[targets.Mayo$Diagnosis!='Control',]
MEcors.AD.PSP <- bicorAndPvalue(MEs.AD.PSP,numericMeta.AD.PSP)
moduleTraitCor.AD.PSP <- MEcors.AD.PSP$bicor
moduleTraitPvalue.AD.PSP <- MEcors.AD.PSP$p
colnames(moduleTraitPvalue.AD.PSP) <- colnames(moduleTraitCor.AD.PSP) <- paste('PSPvsAD',colnames(moduleTraitCor.AD.PSP),sep='\n')
#combine together
moduleTraitCor.Mayo <- cbind(moduleTraitCor.PSP,moduleTraitCor.AD,moduleTraitCor.AD.PSP)
moduleTraitPvalue.Mayo <- cbind(moduleTraitPvalue.PSP,moduleTraitPvalue.AD,moduleTraitPvalue.AD.PSP)

#MSSM#
#CDR, Age, Sex
MEcors.MSSM <- bicorAndPvalue(MEs_MSSM,targets.MSSM[,c('CDR','AOD','SEX')])
moduleTraitCor.MSSM <- MEcors.MSSM$bicor
moduleTraitPvalue.MSSM <- MEcors.MSSM$p
colnames(moduleTraitCor.MSSM) <- colnames(moduleTraitPvalue.MSSM) <- paste(colnames(moduleTraitCor.MSSM),'MSSM',sep='.')
#diagnosis AD possible
adMat=MEs_MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==1,]
cond1=as.numeric(targets.MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==1,'DxByPath'])
MEcors.AD.poss_MSSM=bicorAndPvalue(adMat,cond1)
#Diagnosis AD probable
adMat=MEs_MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==2,]
cond1=as.numeric(targets.MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==2,'DxByPath'])
MEcors.AD.prob_MSSM=bicorAndPvalue(adMat,cond1)
#Diagnosis AD definite
adMat=MEs_MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==3,]
cond1=as.numeric(targets.MSSM[targets.MSSM$DxByPath==0|targets.MSSM$DxByPath==3,'DxByPath'])
MEcors.AD.def_MSSM=bicorAndPvalue(adMat,cond1)
#combine together
moduleTraitCor.MSSM <- cbind(as.numeric(MEcors.AD.poss_MSSM$bicor[,1]),as.numeric(MEcors.AD.prob_MSSM$bicor[,1]),as.numeric(MEcors.AD.def_MSSM$bicor[,1]),moduleTraitCor.MSSM)
moduleTraitPvalue.MSSM  <- cbind(as.numeric(MEcors.AD.poss_MSSM$p[,1]),as.numeric(MEcors.AD.prob_MSSM$p[,1]),as.numeric(MEcors.AD.def_MSSM$p[,1]),moduleTraitPvalue.MSSM)
colnames(moduleTraitCor.MSSM)[c(1:3)] <- colnames(moduleTraitPvalue.MSSM)[c(1:3)] <- c('AD.possible','AD.prob','AD.def')

#BLSA#
#set up diagnosis
group=targets.BLSA
group[group$AsymAD==1,'AsymAD'] <- 'AsymAD'
group[group$AD==1,'AD'] <- 'AD'
group[group$CT==1,'CT'] <- 'Control'
targets.BLSA$Diagnosis=targets.BLSA$CT
ind=which(group$CT=='Control')
targets.BLSA$Diagnosis[ind] <- 'Control'
ind=which(group$AsymAD=='AsymAD')
targets.BLSA$Diagnosis[ind] <- 'AsymAD'
ind=which(group$AD=='AD')
targets.BLSA$Diagnosis[ind] <- 'AD'
#age, sex
MEcors.BLSA <- bicorAndPvalue(MEs_BLSA,targets.BLSA[,c('AGE','SEX')])
moduleTraitCor.BLSA <- MEcors.BLSA$bicor
moduleTraitPvalue.BLSA <- MEcors.BLSA$p
colnames(moduleTraitCor.BLSA) <- colnames(moduleTraitPvalue.BLSA) <- paste(colnames(moduleTraitCor.BLSA),'BLSA',sep='.')
#Control vs AD
adMat=MEs_BLSA[targets.BLSA$Diagnosis=='Control'|targets.BLSA$Diagnosis=='AD',]
cond1=factor(targets.BLSA[targets.BLSA$Diagnosis=='Control'|targets.BLSA$Diagnosis=='AD','Diagnosis'],c('Control','AD'))
cond1=as.numeric(cond1)-1
MEcors.AD_BLSA=bicorAndPvalue(adMat,cond1)
#Control vs AsymAD
asyMat=MEs_BLSA[targets.BLSA$Diagnosis=='Control'|targets.BLSA$Diagnosis=='AsymAD',]
cond1=factor(targets.BLSA[targets.BLSA$Diagnosis=='Control'|targets.BLSA$Diagnosis=='AsymAD','Diagnosis'],c('Control','AsymAD'))
cond1=as.numeric(cond1)-1
MEcors.ASYM=bicorAndPvalue(asyMat,cond1)
#combine together
moduleTraitCor.BLSA <- cbind(as.numeric(MEcors.ASYM$bicor[,1]),as.numeric(MEcors.AD_BLSA$bicor[,1]),moduleTraitCor.BLSA)
moduleTraitPvalue.BLSA  <- cbind(as.numeric(MEcors.ASYM$p[,1]),as.numeric(MEcors.AD_BLSA$p[,1]),moduleTraitPvalue.BLSA)
colnames(moduleTraitCor.BLSA)[c(1,2)] <- colnames(moduleTraitPvalue.BLSA)[c(1,2)] <- c('AsymAD.BLSA','AD.BLSA')

#ACT#
#setup diagnosis
targets.ACT$Diagnosis=targets.ACT$CT.Asym.AD
targets.ACT$Diagnosis[targets.ACT$Diagnosis==0] <- 'Control'
targets.ACT$Diagnosis[targets.ACT$Diagnosis==1] <- 'AsymAD'
targets.ACT$Diagnosis[targets.ACT$Diagnosis==2] <- 'AD'

#age, sex
MEcors.ACT <- bicorAndPvalue(MEs_ACT,targets.ACT[,c('Age','Gender')])
moduleTraitCor.ACT <- MEcors.ACT$bicor
moduleTraitPvalue.ACT <- MEcors.ACT$p
colnames(moduleTraitCor.ACT) <- colnames(moduleTraitPvalue.ACT) <- paste(colnames(moduleTraitCor.ACT),'ACT',sep='.')
#AD vs Control
adMat=MEs_ACT[targets.ACT$Diagnosis=='Control'|targets.ACT$Diagnosis=='AD',]
cond1=factor(targets.ACT[targets.ACT$Diagnosis=='Control'|targets.ACT$Diagnosis=='AD','Diagnosis'],c('Control','AD'))
cond1=as.numeric(cond1)-1
MEcors.AD_ACT=bicorAndPvalue(adMat,cond1)
#AsymAD vs Control
asyMat=MEs_ACT[targets.ACT$Diagnosis=='Control'|targets.ACT$Diagnosis=='AsymAD',]
cond1=factor(targets.ACT[targets.ACT$Diagnosis=='Control'|targets.ACT$Diagnosis=='AsymAD','Diagnosis'],c('Control','AsymAD'))
cond1=as.numeric(cond1)-1
MEcors.ASYM=bicorAndPvalue(asyMat,cond1)
#combine together
moduleTraitCor.ACT <- cbind(as.numeric(MEcors.ASYM$bicor[,1]),as.numeric(MEcors.AD_ACT$bicor[,1]),moduleTraitCor.ACT)
moduleTraitPvalue.ACT  <- cbind(as.numeric(MEcors.ASYM$p[,1]),as.numeric(MEcors.AD_ACT$p[,1]),moduleTraitPvalue.ACT)
colnames(moduleTraitCor.ACT)[c(1,2)] <- colnames(moduleTraitPvalue.ACT)[c(1,2)] <- c('AsymAD.ACT','AD.ACT')

#Banner#
#match sample IDs
rownames(MEs_Banner)=targets.Banner$SampleID
targets.Banner=targets.Banner[complete.cases(targets.Banner$Diagnosis),]
MEs_Banner=MEs_Banner[na.omit(match(targets.Banner$SampleID,rownames(MEs_Banner))),]

#age, gender
MEcors.Banner <- bicorAndPvalue(MEs_Banner,targets.Banner[,c('age','gender')])
moduleTraitCor.Banner <- MEcors.Banner$bicor
moduleTraitPvalue.Banner <- MEcors.Banner$p
colnames(moduleTraitCor.Banner) <- colnames(moduleTraitPvalue.Banner) <- paste(colnames(moduleTraitCor.Banner),'Banner',sep='.')

#AD vs control
adMat=MEs_Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='AD',]
cond1=factor(targets.Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='AD','Diagnosis'],c('CTL','AD'))
cond1=as.numeric(cond1)-1
MEcors.AD_Banner=bicorAndPvalue(adMat,cond1)
#AsymAD vs Control
asyMat=MEs_Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='ASYM',]
cond1=factor(targets.Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='ASYM','Diagnosis'],c('CTL','ASYM'))
cond1=as.numeric(cond1)-1
MEcors.ASYM=bicorAndPvalue(asyMat,cond1)
#MCI vs Control
MCIMat=MEs_Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='MCI',]
cond1=factor(targets.Banner[targets.Banner$Diagnosis=='CTL'|targets.Banner$Diagnosis=='MCI','Diagnosis'],c('CTL','MCI'))
cond1=as.numeric(cond1)-1
MEcors.MCI=bicorAndPvalue(MCIMat,cond1)
#combine together
moduleTraitCor.Banner <- cbind(as.numeric(MEcors.ASYM$bicor[,1]),as.numeric(MEcors.MCI$bicor[,1]),as.numeric(MEcors.AD_Banner$bicor[,1]),moduleTraitCor.Banner)
moduleTraitPvalue.Banner  <- cbind(as.numeric(MEcors.ASYM$p[,1]),as.numeric(MEcors.MCI$p[,1]),as.numeric(MEcors.AD_Banner$p[,1]),moduleTraitPvalue.Banner)
colnames(moduleTraitCor.Banner)[c(1:3)] <- colnames(moduleTraitPvalue.Banner)[c(1:3)] <- c('AsymAD.Banner','MCI.Banner','AD.Banner')

##Combine all cohorts together
moduleTraitCor=cbind(moduleTraitCor.MSSM,moduleTraitCor.BLSA,moduleTraitCor.ACT,moduleTraitCor.Banner,moduleTraitCor.Mayo)
moduleTraitPvalue=cbind(moduleTraitPvalue.MSSM,moduleTraitPvalue.BLSA,moduleTraitPvalue.ACT,moduleTraitPvalue.Banner,moduleTraitPvalue.Mayo)

moduleTraitCor=moduleTraitCor[-11,] ##remove grey
moduleTraitPvalue=moduleTraitPvalue[-11,]
#add MSSM name
colnames(moduleTraitCor)[c(1:3)]=paste(colnames(moduleTraitCor)[c(1:3)],'MSSM',sep='\n')
colnames(moduleTraitPvalue)[c(1:3)]=paste(colnames(moduleTraitPvalue)[c(1:3)],'MSSM',sep='\n')

#plot module eigengene and trait heatmap
pdf(paste0(resdir,'ME_TraitCor_Consensus.pdf'),width=26,height=18)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
colvec <- rep("white",100)
colvec[11:50] <- "pink"
colvec[1:10] <- "red"

labeledHeatmap(Matrix = moduleTraitPvalue,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = colvec,
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.8,
    zlim = c(0,0.1), #was 0,1
    main = paste("Module-trait relationships\n bicor r-value \n (p-value)"),
    cex.main=0.8)
dev.off()

