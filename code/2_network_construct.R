## 2. Network Construction

#resdir from 1_data_prep.R

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,40, by=2));

# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);

# Call the network topology analysis function for each set in turn
for (set in 1:nSets){
    powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
        verbose = 5,networkType="signed",corFnc="bicor")[[2]]);
}

# Plot the results:
pdf(paste0(resdir,"1_Power.pdf"), height=10, width=18)

colors = c("blue", "red","green","black",'magenta')
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets) {
    for (col in 1:length(plotCols)){
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) { 
    for (set in 1:nSets) {
    
        if (set==1) {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                main = colNames[col]);
            addGrid();
        }
    
        if (col==1) {
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                labels=powers,cex=cex1,col=colors[set]);
        } else {
            text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                labels=powers,cex=cex1,col=colors[set]);
        }
    
        if (col==1){
            legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
        } else {
            legend("topright", legend = setLabels, col = colors, pch = 20) ;
        }
    }
}
dev.off()


#choose soft power 
softPower=16

#create consensus network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                                         maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                                         randomSeed = 12345,
                                         corType = "pearson", ## no use for bicor
                                         power = softPower,
                                         consensusQuantile = 0.2,
                                         networkType = "signed",
                                         TOMType = "unsigned",
                                         TOMDenom = "min",
                                         scaleTOMs = TRUE, scaleQuantile = 0.2,
                                         sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                         useDiskCache = TRUE, chunkSize = NULL,
                                         deepSplit = 2,
                                         detectCutHeight = 0.99999999, minModuleSize = 20,
                                         mergeCutHeight = 0.07,
                                         saveConsensusTOMs = TRUE,
                                         consensusTOMFilePattern = paste0(resdir,"ConsensusTOM-block.%b.rda"))


#assign consensus module eigengene, colors and dendrogram to new variables
consMEs = net$multiMEs;
moduleColors = net$colors;
consTree = net$dendrograms[[1]];

# Various Tree Cutting Params
consTree= hclust(1-consTomDS,method="average");

mColorh <- mLabelh <- colorLabels <- NULL
for (minModSize in c(20,40,100)) {
    for (dthresh in c(0.07,0.1,0.2,0.25)) {
        for (ds in c(2,4)) {
            print("Trying parameters:")
            print(c(minModSize,dthresh,ds))
            tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
                minClusterSize = minModSize, cutHeight = 0.99999999,
                deepSplit = ds, distM = as.matrix(1-consTomDS))

            merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                cutHeight = dthresh)
            mColorh <- cbind(mColorh,labels2colors(merged$colors))
            mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
        }
    }
}


pdf(paste0(resdir,"ConsensusTOM_MultiDendro.pdf"),height=30,width=35)
plotDendroAndColors(consTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed consensus network with power = 12"));
dev.off()


### Final network 

#choose the following parameters
mms=20
ds=4
dthresh=0.07


tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
    minClusterSize = mms, cutHeight = 0.99999999,
    deepSplit = ds, distM = as.matrix(1-consTomDS))

merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels, cutHeight = dthresh)
moduleColor.cons <- labels2colors(merged$colors)

mColorh <- cbind(labels2colors(merged$colors))
mLabelh <- c("Merged Colors")



pdf(paste0(resdir,"ConsensusTOM_FinalDendro.pdf"),height=15,width=25)
plotDendroAndColors(consTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,ylim=c(0.8,1),main= paste("Signed bicor network with power = 16, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()



