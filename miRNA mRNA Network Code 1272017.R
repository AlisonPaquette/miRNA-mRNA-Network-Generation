##########################################
##title: "Improved Code Correlations"
##author: "Alison Paquette"
##date: "10/23/2017"
##########################################

#Step 1: Visit Targetscan Databae (http://www.targetscan.org/vert_71/) and download .txt files for miRNAs of Interest

#Step 2: Load .txt files into R as a List
directory<-"~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/TargetScanFiles/"
setwd("~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/TargetScanFiles/")
Files <- list.files(directory)
length(Files) #Does this match up with the number of files you would expect?
Files2<-gsub(".txt","",Files)
MeanContexts<-data.frame(matrix(NA,nrow=length(Files),ncol=5))
rownames(MeanContexts)<-Files2
colnames(MeanContexts)<-c("MedianContextScore","initial","cutoff")
miRNA_mRNA<-list()
for (i in 1:length(Files)){
  tmp = read.table(file=Files[i],header=T,fill=T, sep="\t",stringsAsFactors = FALSE,na.strings="N/A")
  tmp<-as.data.frame(cbind(as.character(tmp$Ortholog.of.target.gene),as.numeric(tmp$Cumulative.weighted.context...score)))
  colnames(tmp)<-c("Gene","CWContextScore")
  tmp<-na.omit(tmp)
  tmp<-subset(tmp,duplicated(tmp$Gene)==F)
  tmp$CWContextScore<-as.numeric(as.character(tmp$CWContextScore))
  rownames(tmp)<-as.character(tmp$Gene)
  hist(tmp$CWContextScore)
  cutoff<-median(tmp$CWContextScore)
  MeanContexts[i,2]<-length(rownames(tmp))
  MeanContexts[i,1]<-cutoff
  tmp2<-subset(tmp,CWContextScore<cutoff)
  MeanContexts[i,3]<-length(rownames(tmp2))
  tmp3<-as.character(tmp2$Gene)
  miRNA_mRNA[[Files2[i]]] <- tmp3
}

###Load miRNA Expression Matrix
load("~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/Pax_miRNA_NormalizeData.RData")
DE_miRNAs<-Files2
miRNA_mat<-NormData[Files2,]

###Load mRNA expression Matrix
load("~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/Pax_MRNAFORDIRAC.RData")

mRNA_Mat<-Pax_MRNA

#Ensure miRNA and mRNA names Match
Samples<-intersect(colnames(mRNA_Mat),colnames(miRNA_mat))
miRNA<-miRNA_mat[,Samples]
mRNA<-as.matrix(mRNA_Mat[,Samples])

#Load Correlation Functions
MakeCorrelations <- function(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal,CorValue){


  #First, Ensure MiRNAs and target mRNAs Match
  if (any((colnames(mRNA)==colnames(miRNA))==F)){
    stop("Sample IDs from mRNA and miRNA Matrixes do not match!")
  }


  d <- data.frame()
  for (i in 1:length(DE_miRNAs)){
    targets<-miRNA_mRNA[i]
    microRNA<-names(targets)
    targets<-targets[[1]]
    #Isolte mRNA expression of interest
    mRNA.Target<-mRNA[intersect(rownames(mRNA),targets),]
    if(is.na(mRNA.Target[1])==T) {break}
    targets2<-rownames(mRNA.Target)
    #miRNA.rep<-(miRNA[colnames(mRNA_miRNA)[i],])
    miRNA.rep<-(miRNA[microRNA,])
    Cor<-as.data.frame(matrix(NA,nrow=length(targets2),ncol=2))

    for (j in 1:length(targets2)){
      test<-cor.test(as.numeric(mRNA.Target[j,]),miRNA.rep,method="spearman",paired=T,exact=F)
      Cor[j,1]<-(test$p.value)
      Cor[j,2]<-(test$estimate)
    }
    colnames(Cor)<-c("P.Unadjust","Corr")
    rownames(Cor)<-targets2
    Cor$P.FDR<-p.adjust(Cor$P.Unadjust, method = "fdr", n = length(Cor$P.Unadjust))
    Cor$miRNA<-rep(microRNA,times=length(rownames(Cor)))

    Cor2<-subset(Cor,Corr<(CorValue))
    # Cor2<-subset(Cor2,Cor$P.Unadjust<PVal)
    Cor2<-subset(Cor2,Cor2$P.Unadjust<PVal)
    d<-rbind(d,Cor2)
  }
  d
}
EvalPredictions <- function(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal,CorValue){

  #First, Ensure MiRNAs and target mRNAs Match
  if (any((colnames(mRNA)==colnames(miRNA))==F)){
    stop("Sample IDs from mRNA and miRNA Matrixes do not match!")
  }


  d <- as.data.frame(matrix(NA,nrow=length(DE_miRNAs),ncol=3))
  rownames(d)<-DE_miRNAs
  colnames(d)<-c("Npredictions","nMRNAsonArray","nSig")

  for (i in 1:length(DE_miRNAs)){
    targets<-miRNA_mRNA[i]
    microRNA<-names(targets)
    targets<-targets[[1]]
    d[i,1]<-  length(targets)
    #Isolte mRNA expression of interest
    mRNA.Target<-mRNA[intersect(rownames(mRNA),targets),]
    if(is.na(mRNA.Target[1])==T) {break}
    targets2<-rownames(mRNA.Target)
    d[i,2]<-  length(targets2)
    miRNA.rep<-(miRNA[microRNA,])


    Cor<-as.data.frame(matrix(NA,nrow=length(targets2),ncol=2))

    for (j in 1:length(targets2)){
      test<-cor.test(as.numeric(mRNA.Target[j,]),miRNA.rep,method="spearman",paired=T,exact=F)
      Cor[j,1]<-(test$p.value)
      Cor[j,2]<-(test$estimate)
    }
    colnames(Cor)<-c("P.Unadjust","Corr")
    rownames(Cor)<-targets2
    Cor$P.FDR<-p.adjust(Cor$P.Unadjust, method = "fdr", n = length(Cor$P.Unadjust))
    Cor$miRNA<-rep(microRNA,times=length(rownames(Cor)))

    Cor2<-subset(Cor,Corr<(CorValue))
    Cor2<-subset(Cor2,Cor2$P.Unadjust<PVal)
    d[i,3]<-  length(rownames(Cor2))
  }
  d
}

#Perform Correlations###
#Series of Spearman Correlation Tests

#This Generates the Correlation Matrix
WBCorMatrix<-MakeCorrelations(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal=0.05,CorValue=-0.3)

#This generates The Data Needed for Cytoscape Visualization
PredictionMatrix<-EvalPredictions(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal=0.05,CorValue=-0.3)
```

#Make Figures
Results<-merge(MeanContexts,PredictionMatrix,by='row.names',all=T)

#Plot Context ++ Score
BarData<-Results[order(Results[,2],decreasing=T),]
#barplot(sort(Results[,2]),ylim=c(-0.8,0.2),ylab="Mean Context ++ Score",col="dodgerblue",main="Whole Blood miRNAs")
par("mar"=c(8,4,1.5,1))
m<-barplot(sort(Results[,2]),ylim=c(-0.8,0.3),ylab="Mean Context ++ Score",col="dodgerblue",main="Whole Blood miRNAs")
m
axis(1,las=3,at=m,labels=as.character(BarData[,1]))


#Load List of DE mRNAs from LIMMA/EdgeR
DE_mRNAs<-read.csv("~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/DE_Pax_TMMNORM.csv")
rownames(DE_mRNAs)<-DE_mRNAs$X
sigNetworks<-merge(DE_mRNAs,WBCorMatrix,by='row.names',all=F)
#write.csv(sigNetworks,file="sigmRNAmiRNAnetworks.csv")

#Create Files For Cytoscape####
flag<-cbind(sigNetworks$Row.names,(rep("PTL_mRNA",times=length(sigNetworks$Row.names))))
colnames(flag)<-c("GENE","Sig")
rownames(flag)<-flag[,1]

forcyto<-merge(WBCorMatrix,flag,by='row.names',all=T)
forcyto$CytoP<-(-(log10(forcyto$P.FDR)))

CytoEdges<-forcyto[,c(1,3,4,5,8)]

CytoNodes1<-cbind(forcyto[,1],rep("Gene",times=length(rownames(CytoEdges))),as.character(forcyto[,7]))
CytoNodes2<-cbind(unique(forcyto$miRNA),rep("miRNA",times=length(unique(forcyto$miRNA))),rep(NA,times=length(unique(forcyto$miRNA))))
CytoNodes<-rbind(CytoNodes1,CytoNodes2)

colnames(CytoNodes)<-c("Node","Type","DEG")

write.csv(CytoNodes,file="~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/CytoNodes.csv")
write.csv(CytoEdges,file="~/Dropbox/OBC_Project/PaxMIRNA/IntegratedAnalysis/CytoEdges.csv")


