###3. plot dendrogram and correlation of gene with diagnosis/covariates 

####create Diagnosis for targets that don't have Diagnosis###
#BLSA
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

#ACT
targets.ACT$Diagnosis=targets.ACT$CT.Asym.AD
targets.ACT$Diagnosis[targets.ACT$Diagnosis==0] <- 'Control'
targets.ACT$Diagnosis[targets.ACT$Diagnosis==1] <- 'AsymAD'
targets.ACT$Diagnosis[targets.ACT$Diagnosis==2] <- 'AD'

#MSSM
targets.MSSM$Diagnosis=targets.MSSM$DxByPath
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==0] <- 'Control'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==1] <- 'AD.poss'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==2] <- 'AD.pro'
targets.MSSM$Diagnosis[targets.MSSM$Diagnosis==3] <- 'AD.def'

#targets.Mayo have appropriate Diagnosis
#targets.Banner have appropriate Diagnosis

###setting up data to plot dendrogram###

#Mayo#
Diagnosis=factor(as.character(targets.Mayo$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=model.matrix(~Diagnosis)
Diagnosis=Diagnosis[,-1]
Age=as.numeric(targets.Mayo$AgeAtDeath)
Gender=targets.Mayo$Gender

geneSigs=matrix(NA,nrow=4,ncol=ncol(datExpr.Mayo)) # create a vector to hold the data

#calculate correlations
for(i in 1:ncol(geneSigs)) {

    exprvec=as.numeric(datExpr.Mayo[,i]) # get the expression vector for ith gene
    condition.ADr=bicor(exprvec, Diagnosis[,1],use="pairwise.complete.obs")
    condition.PSPr=bicor(exprvec, Diagnosis[,2],use="pairwise.complete.obs")
    ager=bicor(Age,exprvec,use="pairwise.complete.obs")
    sexr=sqrt(max(summary(lm(exprvec~as.factor(Gender)))$adj.r.squared,0)) # calculate adjusted R^2s square-root for categorical variables
    geneSigs[,i]=c(condition.ADr, condition.PSPr, ager, sexr)
}
#assign color
geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) # For categorical variables  we do not want values, thus lim=c(0,1), and signed and centered=F

#create rownames
rownames(geneSigs)=c("Mayo.AD","Mayo.PSP","Age","Sex")
#reassign variable name
geneSigs.Mayo=geneSigs
rm(geneSigs)


#MSMM#
cat('MSSM\n')
Diagnosis=factor(as.character(targets.MSSM$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=model.matrix(~Diagnosis)
Diagnosis=Diagnosis[,-1]
Age=as.numeric(targets.MSSM$AOD)
Gender=as.numeric(targets.MSSM$SEX)

geneSigs=matrix(NA,nrow=5,ncol=ncol(datExpr.MSSM)) # create a vector to hold the data

#calculate correlation
for(i in 1:ncol(geneSigs)) {

	exprvec=as.numeric(datExpr.MSSM[,i]) # get the expression vector for ith gene
	condition.ADdefr=bicor(exprvec, Diagnosis[,1],use="pairwise.complete.obs")
    condition.ADpossr=bicor(exprvec, Diagnosis[,2],use="pairwise.complete.obs")
    condition.ADprobr=bicor(exprvec, Diagnosis[,3],use="pairwise.complete.obs")
    ager=bicor(Age,exprvec,use="pairwise.complete.obs")
    sexr=bicor(Gender,exprvec,use="pairwise.complete.obs")

    geneSigs[,i]=c(condition.ADpossr, condition.ADprobr,condition.ADdefr, ager, sexr)
}
#assign to color
geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[5,] =numbers2colors(as.numeric(geneSigs[5,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

#create row names
rownames(geneSigs)=c("MSSM.AD.poss","MSSM.AD.prob","MSSM.AD.def","Age","Sex")
#reassign variable names
geneSigs.MSSM=geneSigs
rm(geneSigs)


# Banner#
cat('Banner\n')
#match Banner datExpr and targets
gnS=intersect(rownames(datExpr.Banner),targets.Banner$SampleID)
datExpr.Banner=datExpr.Banner[match(gnS,rownames(datExpr.Banner)),]
targets.Banner=targets.Banner[match(gnS,targets.Banner$SampleID),]

Diagnosis=as.character(targets.Banner$Diagnosis) #TC
Diagnosis[is.na(Diagnosis)]='Unknown' #set the NA to unknown. has to match datExpr #TC
Diagnosis=factor(Diagnosis) #TC
Diagnosis<-relevel(Diagnosis,'CTL')
Diagnosis=model.matrix(~Diagnosis) #DiagnosisUnknown is last. So doesn't change anything
Diagnosis=Diagnosis[,-1]
Age=as.numeric(targets.Banner$age)
Gender=as.numeric(targets.Banner$gender)

geneSigs=matrix(NA,nrow=5,ncol=ncol(datExpr.Banner)) # create a vector to hold the data
#correlation
for(i in 1:ncol(geneSigs)) {

  exprvec=as.numeric(datExpr.Banner[,i]) # get the expression vector for ith gene
  condition.ASYMr=bicor(exprvec, Diagnosis[,2],use="pairwise.complete.obs")
  condition.MCIr=bicor(exprvec, Diagnosis[,3],use="pairwise.complete.obs")
  condition.ADr=bicor(exprvec, Diagnosis[,1],use="pairwise.complete.obs")
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(Gender,exprvec,use="pairwise.complete.obs")
  geneSigs[,i]=c(condition.MCIr,condition.ASYMr,condition.ADr,ager, sexr)
}

#assign to color
geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[5,] =numbers2colors(as.numeric(geneSigs[5,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

#create row name
rownames(geneSigs)=c("Banner.MCI","Banner.AsymAD","Banner.AD","Age","Sex")
#reassign variable name
geneSigs.Banner=geneSigs
rm(geneSigs)

#ACT#
Diagnosis=factor(as.character(targets.ACT$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=model.matrix(~Diagnosis)
Diagnosis=Diagnosis[,-1]
Age=as.numeric(targets.ACT$Age)
Gender=as.numeric(targets.ACT$Gender)

geneSigs=matrix(NA,nrow=4,ncol=ncol(datExpr.ACT)) # create a vector to hold the data

#correlation
for(i in 1:ncol(geneSigs)) {
 	exprvec=as.numeric(datExpr.ACT[,i]) # get the expression vector for ith gene
 	condition.AsymADr=bicor(exprvec, Diagnosis[,2],use="pairwise.complete.obs")
    condition.ADr=bicor(exprvec, Diagnosis[,1],use="pairwise.complete.obs")
    ager=bicor(Age,exprvec,use="pairwise.complete.obs")
    sexr=bicor(Gender,exprvec,use="pairwise.complete.obs")

    geneSigs[,i]=c(condition.AsymADr, condition.ADr, ager, sexr)
}

#assign to color
geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

#create row name
rownames(geneSigs)=c("ACT.AsymAD","ACT.AD","Age","Sex")
#reassign variable name
geneSigs.ACT=geneSigs
rm(geneSigs)

#BLSA#
Diagnosis=factor(as.character(targets.BLSA$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=model.matrix(~Diagnosis)
Diagnosis=Diagnosis[,-1]
Age=as.numeric(targets.BLSA$AGE)
Gender=as.numeric(targets.BLSA$SEX)

geneSigs=matrix(NA,nrow=4,ncol=ncol(datExpr.BLSA)) # create a vector to hold the data
#correlation
for(i in 1:ncol(geneSigs)) {

    exprvec=as.numeric(datExpr.BLSA[,i]) # get the expression vector for ith gene
    condition.AsymADr=bicor(exprvec, Diagnosis[,2],use="pairwise.complete.obs")
    condition.ADr=bicor(exprvec, Diagnosis[,1],use="pairwise.complete.obs")
    ager=bicor(Age,exprvec,use="pairwise.complete.obs")
    sexr=bicor(Gender,exprvec,use="pairwise.complete.obs")

    geneSigs[,i]=c(condition.AsymADr, condition.ADr, ager, sexr)
}

#assign to color
geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

#create row name
rownames(geneSigs)=c("BLSA.AsymAD","BLSA.AD","Age","Sex")
#reassign variable name
geneSigs.BLSA=geneSigs
rm(geneSigs)


###plot dendrogram###
mColorh1=cbind(mColorh,geneSigs.Mayo[1,],geneSigs.Mayo[2,],geneSigs.Mayo[3,],geneSigs.Mayo[4,],
    mColorh, geneSigs.MSSM[1,],geneSigs.MSSM[2,],geneSigs.MSSM[3,],geneSigs.MSSM[4,],geneSigs.MSSM[5,],
    geneSigs.Banner[1,],geneSigs.Banner[2,],geneSigs.Banner[3,],geneSigs.Banner[4,],geneSigs.Banner[5,],
    geneSigs.ACT[1,],geneSigs.ACT[2,],geneSigs.ACT[3,],geneSigs.ACT[4,],
    geneSigs.BLSA[1,],geneSigs.BLSA[2,],geneSigs.BLSA[3,],geneSigs.BLSA[4,])

mLabelh1=c(mLabelh,rownames(geneSigs.Mayo),rownames(geneSigs.MSSM),rownames(geneSigs.Banner),rownames(geneSigs.ACT),rownames(geneSigs.BLSA))

pdf(paste0(resdir,"ConsensusTOM_FinalDendro_geneSigs.pdf"),height=8,width=15)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 12, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()


