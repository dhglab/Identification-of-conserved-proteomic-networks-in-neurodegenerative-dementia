###create plots of each module
library(igraph);
library(RColorBrewer);

#KMEs
KME.vis= consensus.KMEs
#get gene names
annot=as.data.frame(do.call("rbind",strsplit(as.character(colnames(datExpr.Banner)),"[|]")))
KME.vis$Symbol=annot$V1
geneInfo.vis=as.data.frame(cbind(KME.vis$Symbol,moduleColor.cons, KME.vis))

geneInfo.vis=geneInfo.vis[,-ncol(geneInfo.vis)] # check if last column is Ensembl gene id
colnames(geneInfo.vis)[1]= "Symbol"
colnames(geneInfo.vis)[2]= "Module.Color"

TOM.matrix = as.matrix(consTomDS);
#Get the top connected genes in the module
uniquemodcolors = unique(moduleColor.cons); #moduleColors
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

pdf(paste0(resdir,'iGraphsCons_Dense.pdf'),height=9,width=10,useDingbats=F)

for (mod in uniquemodcolors)  {
    numgenesingraph = 100;
    numconnections2keep = 2000; 
    cat('module:',mod,'\n');
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="NA",]
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="",]
    colind = which(colnames(geneInfo.vis)== paste("consensus.kME",mod,sep=""));
    rowind = which(geneInfo.vis[,2]==mod);
    cat(' ',length(rowind),'probes in module\n');
    submatrix = geneInfo.vis[rowind,];
    orderind = order(submatrix[,colind],decreasing=TRUE);
    if (length(rowind) < numgenesingraph) {
      numgenesingraph = length(rowind);
      numconnections2keep = numgenesingraph/2 * (numgenesingraph/6 - 1); #added /2 and /6 9/14/2015
    }
    cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
    submatrix = submatrix[orderind[1:numgenesingraph],];

    #Identify the columns in the TOM that correspond to these hub probes
    matchind = match(submatrix$Symbol,annot$Symbol);
    reducedTOM = TOM.matrix[matchind,matchind];

    orderind = order(reducedTOM,decreasing=TRUE);
    connections2keep = orderind[1:numconnections2keep];
    reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] = 1;

    g0 <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMata <- layout.circle(g0)

    if (ncol(reducedTOM) < 51) {
        g0 <- graph.adjacency(as.matrix(reducedTOM[11:ncol(reducedTOM),11:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMatb <- layout.circle(g0)

        g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75)
    } else {

        g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMatb <- layout.circle(g0)

        g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMatc <- layout.circle(g0)
        g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc)
    }

    plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(submatrix$Symbol),vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*8,main=paste(mod,"module"))

}
dev.off();




