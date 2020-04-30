#1. prepare data for WGCNA

library(WGCNA)

#results directory
#resdir

#load R data files so you have the following datExpr and targets variables for Mayo, MSSM (Mt Sinait), BLSA, ACT, Banner
#datExpr.*cohort* is expression data - sample (row) x expression (columns)
#targets.*cohort* is metadata - sample (row) x metadata (columns)

#Mayo
#datExpr.Mayo
#targets.Mayo

#MSMM
#datExpr.MSSM
#targets.MSSM

#BLSA
#datExpr.BLSA
#targets.BLSA

#ACT
#datExpr.ACT
#targets.ACT

#Banner
#datExpr.Banner
#targets.Banner

#get proteins that are in all cohorts
gnS=intersect(intersect(colnames(datExpr.MSSM),colnames(datExpr.BLSA)),intersect(colnames(datExpr.ACT),colnames(datExpr.Banner))) #2142
gnS=intersect(gnS,colnames(datExpr.Mayo)) ##2005

#keep proteins from each cohort that are presentin all cohorts
datExpr.Mayo=datExpr.Mayo[,match(gnS,colnames(datExpr.Mayo))]
datExpr.MSSM=datExpr.MSSM[,match(gnS,colnames(datExpr.MSSM))]
datExpr.BLSA=datExpr.BLSA[,match(gnS,colnames(datExpr.BLSA))]
datExpr.ACT=datExpr.ACT[,match(gnS,colnames(datExpr.ACT))]
datExpr.Banner=datExpr.Banner[,match(gnS,colnames(datExpr.Banner))]

# Read the data and make a multiexpression set
nSets=5
setLabels=c("MSSM","BLSA","ACT","Banner","Mayo")
shortLabels=c("MSSM","BLSA","ACT","Banner","Mayo")

multiExpr=vector(mode="list",length=nSets)

multiExpr[[1]] = list(data=as.data.frame(datExpr.MSSM)) # MSSM
names(multiExpr[[1]]$data)=colnames(datExpr.MSSM)
rownames(multiExpr[[1]]$data)=rownames(datExpr.MSSM)

multiExpr[[2]] = list(data=as.data.frame(datExpr.BLSA)) # BLSA
names(multiExpr[[2]]$data)=colnames(datExpr.BLSA)
rownames(multiExpr[[2]]$data)=rownames(datExpr.BLSA)

multiExpr[[3]] = list(data=as.data.frame(datExpr.ACT)) #ACT
names(multiExpr[[3]]$data)=colnames(datExpr.ACT)
rownames(multiExpr[[3]]$data)=rownames(datExpr.ACT)

multiExpr[[4]] = list(data=as.data.frame(datExpr.Banner)) #Banner
names(multiExpr[[4]]$data)=colnames(datExpr.Banner)
rownames(multiExpr[[4]]$data)=rownames(datExpr.Banner)

multiExpr[[5]] = list(data=as.data.frame(datExpr.Mayo)) # Mayo
names(multiExpr[[5]]$data)=colnames(datExpr.Mayo)
rownames(multiExpr[[5]]$data)=rownames(datExpr.Mayo)

checkSets(multiExpr) # check data size

#create multi meta-data
multiMeta=list(MSSM =list(data=targets.MSSM),BLSA =list(data=targets.BLSA),ACT=list(data=targets.ACT),Banner=list(data=targets.Banner),Mayo=list(data=targets.Mayo))

#######

