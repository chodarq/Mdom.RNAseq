#### 0. Load libraries ####
library("limma")
library("edgeR")
library("RColorBrewer")
library("Heatplus")
library("clValid")
library("ggplot2")
library("gplots")
library("affycoretools")
library("apcluster")
library("cluster")
library("Mfuzz")
library("biomaRt")
library("topGO")
library("data.table")
library("plyr")
#### 1. Read counts from file and transform to DGElist object ####

mdcounts<-read.csv("counts.csv",sep="\t",header = TRUE,row.names = 1) #read total counts
groups<-c("SB","SB","CB","CB","G","G","GBE","GBE","DC","DC") # grouped by stage
mdomDGE<-DGEList(counts=mdcounts,group = groups)

#### 2. Calc for Normalization factors by TMM (Trimmed Mean of M-values) ####
normMdomDGE<-calcNormFactors(mdomDGE,method = c("TMM"))
plotMDS(normMdomDGE,main="Mdom")
cpmMdom<-cpm(normMdomDGE, normalized.lib.sizes=TRUE) # stored CPM counts

#### 3. Global data clustering ####
cpmMdom.mean <-  cbind(apply(cpmMdom[,1:2],1,mean),apply(cpmMdom[,3:4],1,mean),
                       apply(cpmMdom[,5:6],1,mean),apply(cpmMdom[,7:8],1,mean),
                       apply(cpmMdom[,9:10],1,mean))
colnames(cpmMdom.mean)<-c("SB","CB","G","GBE","DC")

## Graph complete cluster ##
dist1<-as.dist(1-cor(t(cpmMdom.mean)))
hc1<-hclust(dist1,method="ward.D2")
figura<-heatmap.2(cpmMdom.mean,dendrogram = c("row"),trace="none",density.info = "none",
          Rowv=as.dendrogram(hc1),labRow=NA,Colv="none",scale = c("row"),
          col=colorRampPalette(brewer.pal(11, "RdBu"))(100))
plot(as.dendrogram(hc1))
rect.hclust(hc1, k = 4, border = c("red","green","blue","black"))
dendro<-cutree(hc1,k=4)
clusterCols <- c("darkred","red","indianred","tomato")
myClusterSideBar <- clusterCols[dendro]
figura2<-heatmap.2(cpmMdom.mean,dendrogram = c("row"),trace="none",density.info = "none",
                   Rowv=as.dendrogram(hc1),labRow=NA,Colv="none",scale = c("row"),
                   col=colorRampPalette(brewer.pal(11, "RdBu"))(100),RowSideColors=myClusterSideBar)

cpmMdom.mean.ord=(cpmMdom.mean)[row.hc1$order,] #reorder your d table according to the heatmap 
grupo1<-cpmMdom.mean[dendro==1,]
grupo2<-cpmMdom.mean[dendro==2,]
grupo3<-cpmMdom.mean[dendro==3,]
grupo4<-cpmMdom.mean[dendro==4,]

g1<-melt(grupo1)
g2<-melt(grupo2)
g3<-melt(grupo3)
g4<-melt(grupo4)

d1<-ggplot(g1,aes(x=Var2,y=log(value,10)))+
  stat_summary(fun.y='mean',geom='line',colour="black",aes(group=1),size=1)+
  labs(title="Gene expression profile",x="Developmental stage",y="log10(cpm)")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("1","2","3","4","5"),labels=c("SB","CB","G","GBE","DC"))+
  geom_boxplot(aes(group=Var2),alpha=0.3,size=0.2,width=0.3,fill="pink")
d2<-ggplot(g2,aes(x=Var2,y=log(value,10)))+
  stat_summary(fun.y='mean',geom='line',colour="black",aes(group=1),size=1)+
  labs(title="Gene expression profile",x="Developmental stage",y="log10(cpm)")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("1","2","3","4","5"),labels=c("SB","CB","G","GBE","DC"))+
  geom_boxplot(aes(group=Var2),alpha=0.3,size=0.2,width=0.3,fill="pink")
d3<-ggplot(g3,aes(x=Var2,y=log(value,10)))+
  stat_summary(fun.y='mean',geom='line',colour="black",aes(group=1),size=1)+
  labs(title="Gene expression profile",x="Developmental stage",y="log10(cpm)")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("1","2","3","4","5"),labels=c("SB","CB","G","GBE","DC"))+
  geom_boxplot(aes(group=Var2),alpha=0.3,size=0.2,width=0.3,fill="pink")
d4<-ggplot(g4,aes(x=Var2,y=log(value,10)))+
  stat_summary(fun.y='mean',geom='line',colour="black",aes(group=1),size=1)+
  labs(title="Gene expression profile",x="Developmental stage",y="log10(cpm)")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("1","2","3","4","5"),labels=c("SB","CB","G","GBE","DC"))+
  geom_boxplot(aes(group=Var2),alpha=0.3,size=0.2,width=0.3,fill="pink")

multiplot(d1, d2, d3, d4, cols=2)
write.csv(grupo1,file="cluster1.csv",)
write.csv(grupo2,file="cluster2.csv")
write.csv(grupo3,file="cluster3.csv")
write.csv(grupo4,file="cluster4.csv")

#### 4. voom  transformation of data ####
design<-model.matrix(~0+groups,data=normMdomDGE$samples)
colnames(design)<-levels(normMdomDGE$samples$group)
v<-voomWithQualityWeights(normMdomDGE$counts,design=design,plot=TRUE)

meanNormMdom <-  cbind(apply(v$E[,1:2],1,mean),apply(v$E[,3:4],1,mean),
                       apply(v$E[,5:6],1,mean),apply(v$E[,7:8],1,mean),
                       apply(v$E[,9:10],1,mean))
colnames(meanNormMdom) <- c("SB","CB","G","GBE","DC")

#### 5. Contrasts testing consecutive stages ####
cont <- makeContrasts(contrasts=c("-SB+CB", "-CB+G", "-G+GBE", "-GBE+DC")
                         ,levels=design) #contrastes por estado
fit<-lmFit(v,design)
fit_contrast<-contrasts.fit(fit,cont)
fit_contrast<-eBayes(fit_contrast)
list<-topTable(fit_contrast,adjust="fdr",sort.by = "none",number = Inf)

results.list<-decideTests(fit_contrast,method="global",p.value = 0.05,adjust.method = "fdr")

par(mfrow=c(1,2))
vennDiagram(results.list,circle.col = c("red", "blue", "green3","yellow"),include="up",main="Up regulated DEs",cex=c(0.8,0.8,0.8))
vennDiagram(results.list,circle.col = c("red", "blue", "green3","yellow"),include="down",main="Down regulated DEs",cex=c(0.8,0.8,0.8))

global.adjusted.p <- fit_contrast$p.value # get individual p-values per contrast
global.adjusted.p[] <- p.adjust(global.adjusted.p, method="BH")
lista.ids.DE<-unique(c(rownames(global.adjusted.p[global.adjusted.p[,1]<0.05,]),
                       rownames(global.adjusted.p[global.adjusted.p[,2]<0.05,]),
                       rownames(global.adjusted.p[global.adjusted.p[,3]<0.05,]),
                       rownames(global.adjusted.p[global.adjusted.p[,4]<0.05,])))

meanNormMdom.DE<-meanNormMdom[rownames(meanNormMdom) %in% lista.ids.DE,] # 2414 cited in manuscript
meanNormMdom.DEexpset<-ExpressionSet(meanNormMdom.DE)
meanNormMdom.DEs<-standardise(meanNormMdom.DEexpset)


#### 6. Graphics ####

cpmMdom.mean <-  cbind(apply(cpmMdom[,1:2],1,mean),apply(cpmMdom[,3:4],1,mean),
                       apply(cpmMdom[,5:6],1,mean),apply(cpmMdom[,7:8],1,mean),
                       apply(cpmMdom[,9:10],1,mean))
colnames(cpmMdom.mean)<-c("SB","CB","G","GBE","DC")

# Graph complete cluster #
d<-as.dist(1-cor(t(cpmMdom.mean)))
hc<-hclust(d,method="ward.D2")
full.cluster<-heatmap.2(cpmMdom.mean,dendrogram = c("row"),trace="none",density.info = "none",
          Rowv=as.dendrogram(hc),labRow=NA,Colv="none",scale = c("row"),key=TRUE,
          col=colorRampPalette(brewer.pal(11, "RdBu"))(100))

two.groups<-cutree(hc,k=2)
plot(hc)
rect.hclust(hc, k = 2, border = 2:5)
group1<-names(two.groups[two.groups==1])
group2<-names(two.groups[two.groups==2])

ids.bs.bc<-rownames(global.adjusted.p[global.adjusted.p[,1]<0.05,])
values.bs.bc<-as.data.frame(meanNormMdom[rownames(meanNormMdom) %in% ids.bs.bc,])
values.bs.bc$GeneId<-rownames(values.bs.bc)
bs.bc.m<-melt(values.bs.bc,id.vars="GeneId")

ids.bc.g<-rownames(global.adjusted.p[global.adjusted.p[,2]<0.05,])
values.bc.g<-as.data.frame(meanNormMdom[rownames(meanNormMdom) %in% ids.bc.g,])
values.bc.g$GeneId<-rownames(values.bc.g)
bc.g.m<-melt(values.bc.g,id.vars="GeneId")

ids.g.gbe<-rownames(global.adjusted.p[global.adjusted.p[,3]<0.05,])
values.g.gbe<-as.data.frame(meanNormMdom[rownames(meanNormMdom) %in% ids.g.gbe,])
values.g.gbe$GeneId<-rownames(values.g.gbe)
g.gbe.m<-melt(values.g.gbe,id.vars="GeneId")

ids.gbe.dc<-rownames(global.adjusted.p[global.adjusted.p[,4]<0.05,])
values.gbe.dc<-as.data.frame(meanNormMdom[rownames(meanNormMdom) %in% ids.gbe.dc,])
values.gbe.dc$GeneId<-rownames(values.gbe.dc)
gbe.dc.m<-melt(values.gbe.dc,id.vars="GeneId")

ggplot(bs.bc.m,aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
ggplot(bc.g.m,aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
ggplot(g.gbe.m,aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
ggplot(gbe.dc.m,aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")

upDE.bs.bc<-values.bs.bc[values.bs.bc$CB>values.bs.bc$SB,]
downDE.bs.bc<-values.bs.bc[values.bs.bc$SB>values.bs.bc$CB,]
p1<-ggplot(melt(upDE.bs.bc,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
p2<-ggplot(melt(downDE.bs.bc,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
multiplot(p1,p2,cols=2)

upDE.bc.g<-values.bc.g[values.bc.g$G>values.bc.g$CB,]
downDE.bc.g<-values.bc.g[values.bc.g$CB>values.bc.g$G,]
p1<-ggplot(melt(upDE.bc.g,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
p2<-ggplot(melt(downDE.bc.g,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
multiplot(p1,p2,cols=2)

upDE.bc.g<-values.bc.g[values.bc.g$G>values.bc.g$CB,]
downDE.bc.g<-values.bc.g[values.bc.g$CB>values.bc.g$G,]
p1<-ggplot(melt(upDE.bc.g,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
p2<-ggplot(melt(downDE.bc.g,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
multiplot(p1,p2,cols=2)

upDE.g.gbe<-values.g.gbe[values.g.gbe$GBE>values.g.gbe$G,]
downDE.g.gbe<-values.g.gbe[values.g.gbe$G>values.g.gbe$GBE,]
p1<-ggplot(melt(upDE.g.gbe,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
p2<-ggplot(melt(downDE.g.gbe,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
multiplot(p1,p2,cols=2)

upDE.gbe.dc<-values.gbe.dc[values.gbe.dc$DC>values.gbe.dc$GBE,]
downDE.gbe.dc<-values.gbe.dc[values.gbe.dc$GBE>values.gbe.dc$DC,]
p1<-ggplot(melt(upDE.gbe.dc,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
p2<-ggplot(melt(downDE.gbe.dc,id.vars="GeneId"),aes(x=variable,y=value,group=GeneId,color=GeneId))+geom_line()+theme(legend.position="none")
multiplot(p1,p2,cols=2)

#### 7. Flymine data analysis ####
md.dm<-read.csv("Md.Dmel.Orthagogue.csv",sep=",",header = TRUE)
dmel.conv<-read.csv("fbgn_fbtr_fbpp.txt.csv",sep=",",header=TRUE)

# Bs.vs.BC
up<-md.dm[md.dm$MdomID %in% rownames(upDE.bs.bc),]
down<-md.dm[md.dm$MdomID %in% rownames(downDE.bs.bc),]
up.dm<-unique(dmel.conv[dmel.conv[,3] %in% up$DmelID,1])
down.dm<-unique(dmel.conv[dmel.conv[,3] %in% down$DmelID,1])
# Bc.vs.G
up<-md.dm[md.dm$MdomID %in% rownames(upDE.bc.g),]
down<-md.dm[md.dm$MdomID %in% rownames(downDE.bc.g),]
up.dm<-unique(dmel.conv[dmel.conv[,3] %in% up$DmelID,1])
down.dm<-unique(dmel.conv[dmel.conv[,3] %in% down$DmelID,1])
# G.vs.GBE
up<-md.dm[md.dm$MdomID %in% rownames(upDE.g.gbe),]
up.dm<-unique(dmel.conv[dmel.conv[,3] %in% up$DmelID,1])
# GBE.vs.DC
up<-md.dm[md.dm$MdomID %in% rownames(upDE.gbe.dc),]
down<-md.dm[md.dm$MdomID %in% rownames(downDE.gbe.dc),]
up.dm<-unique(dmel.conv[dmel.conv[,3] %in% up$DmelID,1])
down.dm<-unique(dmel.conv[dmel.conv[,3] %in% down$DmelID,1])

#### topGO ####
geneNames<-gos.mdom$ensembl_gene_id
geneID2GO <- readMappings(file = "md.geneid2go.map")

geneList<-factor(as.integer(geneNames %in% lista.ids.DE))
names(geneList)<-geneNames

geneListbsbc<-factor(as.integer(geneNames %in% ids.bs.bc))
names(geneListbsbc)<-geneNames
geneListbcg<-factor(as.integer(geneNames %in% ids.bc.g))
names(geneListbcg)<-geneNames
geneListggbe<-factor(as.integer(geneNames %in% ids.g.gbe))
names(geneListggbe)<-geneNames
geneListgbedc<-factor(as.integer(geneNames %in% ids.gbe.dc))
names(geneListgbedc)<-geneNames

GOdata <- new("topGOdata", ontology = "BP", description="DE genes",nodeSize=20,allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.1 <- new("topGOdata", ontology = "BP", description="Bs vs Bc",allGenes = geneListbsbc,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.2 <- new("topGOdata", ontology = "BP", description="Bc vs G",allGenes = geneListbcg,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.3 <- new("topGOdata", ontology = "BP", description="G vs GBE",allGenes = geneListggbe,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.4 <- new("topGOdata", ontology = "BP", description="GBE vs DC",allGenes = geneListgbedc,annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFisher<-runTest(GOdata,algorithm = "elim",statistic="fisher")
resultFisher1<-runTest(GOdata.1,algorithm = "elim",statistic="fisher")
resultFisher2<-runTest(GOdata.2,algorithm = "elim",statistic="fisher")
resultFisher3<-runTest(GOdata.3,algorithm = "elim",statistic="fisher")
resultFisher4<-runTest(GOdata.4,algorithm = "elim",statistic="fisher")

#### Nuevos cruces ####
commons.dm.md<-meanNormMdom.DE[rownames(meanNormMdom.DE) %in% md.dm$MdomID,] #1770 common proteins between 2414 DEs of Mdom and Orthagogue result
subset.orthagogue<-md.dm[md.dm[,1] %in% rownames(commons.dm.md),c(1,2)] # cross between 1770 commons and fbgn_fbpp drosophila conversion. n=3862

subsetOrt.cruce.FBids<-(dmel.conv[dmel.conv[,3] %in% subset.orthagogue$DmelID,c(1,2,3)]) # 3862 Dmel proteins transformed to genes
MDs.a.genesDM<-(unique(subsetOrt.cruce.FBids)) # final set of 1841 genes


#### Functions ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}