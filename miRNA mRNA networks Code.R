##############################################################
#######     miRNA & mRNA Correlation Networks           #####
########       Alison G Paquette, January 3rd 2017  ##########
##############################################################


#Part 0A: Load Packages and Functions To Complete Analysis

#Part 0B: Load & Prepare All Necessary Data
#1.  Need Normalized miRNA and mRNA Gene expression Martixes (from microarray or RNAseq), and Genes of Interest derived from EdgeR, DESEQ2, or LIMMA

#miRNA Matrix
load("/Users/alisonpaquette/Documents/Sickkidsptb/MONOCYTE_INTEGRATED/Monocyte_TMMNormalizedData.RData")
M_miRNA<-NormData
#miRNA Genes of Interest
DE_miRNAS1<-read.csv("~/Documents/Sickkidsptb/M_Integrated/DE_mIRNAs.csv")
DE_miRNAS<-as.character(DE_miRNAS1[,1]) #Vector of N miRNAs significantly assocatied with Phenotype

#mRNA Matrix
load("/Users/alisonpaquette/Documents/Sickkidsptb/MONOCYTE_INTEGRATED/mRNA_M_TMMNorm.RData")

#mRNA Genes of Interest
DE_mRNAs1<-read.csv("~/Documents/Sickkidsptb/MONOCYTE_INTEGRATED/DE_M_MRNA_TMM.csv")
rownames(DE_mRNAs1)<-DE_mRNAs1[,1]
DE_mRNAs<-DE_mRNAs1[,1]


#Load miRNA prediction Database: Here we use Targetscan
mRNA_miRNA<-read.csv("~/Documents/Sickkidsptb/MONOCYTE_INTEGRATED/miRNAMRNAMonocytesTargetScan.csv",na.strings=c("", " "))
DE_miRNAS<-t(mRNA_miRNA[1,])
DE_miRNAS<-as.character(DE_miRNAS[,1])
colnames(mRNA_miRNA)<-DE_miRNAS
mRNA_miRNA<-mRNA_miRNA[-1,]

#Ensure miRNA and mRNA names Match
Samples<-intersect(colnames(M_MRNA),colnames(M_miRNA))
miRNA<-M_miRNA[,Samples]
mRNA<-as.matrix(M_MRNA[,Samples])

#Generate Correlation Network from mRNA data
MakeCorrelations <- function(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal,CorValue){

#First, Ensure MiRNAs and target mRNAs Match
if (any((colnames(mRNA)==colnames(miRNA))==F)){
    stop("Sample IDs from mRNA and miRNA Matrixes do not match!")
  }


d <- data.frame()
for (i in 1:length(DE_miRNAS)){
  targets<-unique(as.character(na.omit(mRNA_miRNA[,i])))
  microRNA<-colnames(mRNA_miRNA)[i]
  length(targets)
  #Isolte mRNA expression of interest
  mRNA.Target<-mRNA[intersect(rownames(mRNA),targets),]
  targets2<-rownames(mRNA.Target)
  miRNA.rep<-(miRNA[colnames(mRNA_miRNA)[i],])

Cor<-as.data.frame(matrix(NA,nrow=length(targets2),ncol=2))

  for (j in 1:length(targets2)){
    test<-cor.test(as.numeric(mRNA.Target[j,]),miRNA.rep,method="pearson",paired=T,exact=F)
    Cor[j,1]<-(test$p.value)
    Cor[j,2]<-(test$estimate)
  }
  colnames(Cor)<-c("P.Unadjust","Corr")
  rownames(Cor)<-targets2
  Cor$P.FDR<-p.adjust(Cor$P.Unadjust, method = "fdr", n = length(Cor$P.Unadjust))
  Cor$miRNA<-rep(microRNA,times=length(rownames(Cor)))

  Cor2<-subset(Cor,Corr<(CorValue))
  Cor2<-subset(Cor2,P.FDR<PVal)
  d<-rbind(d,Cor2)
}
d
}
EvalPredictions <- function(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal,CorValue){

  #First, Ensure MiRNAs and target mRNAs Match
  if (any((colnames(mRNA)==colnames(miRNA))==F)){
    stop("Sample IDs from mRNA and miRNA Matrixes do not match!")
  }


  d <- as.data.frame(matrix(NA,nrow=length(DE_miRNAS),ncol=3))
  rownames(d)<-DE_miRNAS
  colnames(d)<-c("Npredictions","nMRNAsonArray","nSig")

  for (i in 1:length(DE_miRNAS)){
    targets<-unique(as.character(na.omit(mRNA_miRNA[,i])))
    microRNA<-colnames(mRNA_miRNA)[i]
    d[i,1]<-  length(targets)
    #Isolte mRNA expression of interest
    mRNA.Target<-mRNA[intersect(rownames(mRNA),targets),]
    targets2<-rownames(mRNA.Target)
    d[i,2]<-  length(targets2)
    miRNA.rep<-(miRNA[colnames(mRNA_miRNA)[i],])


    Cor<-as.data.frame(matrix(NA,nrow=length(targets2),ncol=2))

    for (j in 1:length(targets2)){
      test<-cor.test(as.numeric(mRNA.Target[j,]),miRNA.rep,method="pearson",paired=T,exact=F)
      Cor[j,1]<-(test$p.value)
      Cor[j,2]<-(test$estimate)
    }
    colnames(Cor)<-c("P.Unadjust","Corr")
    rownames(Cor)<-targets2
    Cor$P.FDR<-p.adjust(Cor$P.Unadjust, method = "fdr", n = length(Cor$P.Unadjust))
    Cor$miRNA<-rep(microRNA,times=length(rownames(Cor)))

    Cor2<-subset(Cor,Corr<(CorValue))
    Cor2<-subset(Cor2,P.FDR<PVal)
    d[i,3]<-  length(rownames(Cor2))
  }
  d
}

MonocytesCorMatrix<-MakeCorrelations(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal=0.1,CorValue=-0.3)
PredictionMatrix<-EvalPredictions(mRNA_miRNA,MRNA,miRNA,DE_miRNAS,PVal=0.1,CorValue=-0.3)

write.csv(PredictionMatrix,file="MonocytePredictionMatrix.csv")
#Identify mRNAS associated with PTL
sigNetworks<-merge(DE_mRNAs1,MonocytesCorMatrix,by='row.names',all=F)
write.csv(sigNetworks,file="sigmRNAmiRNAnetworks.csv")

#CreateCytoscapeNodeFile
flag<-cbind(sigNetworks$Row.names,(rep("PTL_mRNA",times=length(sigNetworks$Row.names))))
colnames(flag)<-c("GENE","Sig")
rownames(flag)<-flag[,1]

forcyto<-merge(MonocytesCorMatrix,flag,by='row.names',all=T)
forcyto$CytoP<-(-(log10(forcyto$P.FDR)))

CytoEdges<-forcyto[,c(3,4,5,8)]

CytoNodes1<-cbind(forcyto[,1],rep("Gene",times=length(rownames(CytoEdges))),as.character(forcyto[,7]))
CytoNodes2<-cbind(unique(forcyto$miRNA),rep("miRNA",times=length(unique(forcyto$miRNA))),rep(NA,times=length(unique(forcyto$miRNA))))
CytoNodes<-rbind(CytoNodes1,CytoNodes2)

colnames(CytoNodes)<-c("Node","Type","DEG")
