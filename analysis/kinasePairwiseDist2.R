setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)
library(plotrix)

source("/home/colinge/net-r/expandGraph.R")

# Reference PPI network
alldb <- read.csv("ppi-20140211.txt",sep="\t",check.names=FALSE,stringsAsFactors=F)
l.no.self <- as.vector(alldb[["source"]]) != as.vector(alldb[["target"]])
from <- as.vector(alldb[["source"]][l.no.self])
to <- as.vector(alldb[["target"]][l.no.self])
ft <- cbind(from,to)
alldb.g <- ugraph(ftM2graphNEL(ft, edgemode="directed"))
n.alldb.g <- nodes(alldb.g)

# Reference KSI network
ksi <- read.table("ksi-20140211.txt",sep="\t",stringsAsFactors=F,header=T)
k.no.self <- as.vector(ksi[["source"]]) != as.vector(ksi[["target"]])
from <- as.vector(ksi[["source"]][k.no.self])
to <- as.vector(ksi[["target"]][k.no.self])
ft <- cbind(from,to)
ksi.g <- ftM2graphNEL(ft, edgemode="directed")
n.ksi.g <- nodes(ksi.g)

# Gene symbols
bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=F,sep="\t",stringsAsFactors=F)
sp2name.1 <- sp2name[sp2name[["Accession"]] %in% n.alldb.g,]
already.symb <- unique(sp2name.1[["Gene name"]])
sp2name.2 <- sp2name[!(sp2name[["Gene name"]]%in%already.symb),]
sp2name.3 <- rbind(sp2name.1,sp2name.2)
sp.names <- sp2name.3[["Gene name"]]
names(sp.names) <- sp2name.3[["Accession"]]
sp.names[sp.names==""] <- sp2name.3[["Accession"]][sp.names==""]
sp.names[bcrabl] <- "BCR-ABL"
names.sp <- sp2name.3[["Accession"]]
names(names.sp) <- sp2name.3[["Gene name"]]
names.sp["BCR-ABL"] <- bcrabl
names.sp["SPFH2"] <- "O94905"
sp.names["Q9Y3C1"] <- "NOP16"
id.sp <- sp2name.3[["Accession"]]
names(id.sp) <- sp2name.3[["Entry name"]]

# Human kinases --------------------------------------------
kinases <- read.csv("UniProt_human_protein_kinases_formated-atypical-1class.txt",sep="\t",stringsAsFactors=F)
kin.in.net <- intersect(kinases[["UniProt_ID"]],n.alldb.g)
rownames(kinases) <- kinases[["UniProt_ID"]]

# Jaccard distance -----------------

getJaccardDistances <- function(nodes,min.connect=5,alpha=0.5){
  ppis <- adj(alldb.g,nodes)
  num.ppis <- unlist(lapply(ppis,length))
  ksis <- adj(ksi.g,nodes)
  num.ksis <- unlist(lapply(ksis,length))
  s <- nodes[num.ppis+num.ksis>=min.connect]
  size <- length(s)
  ppi <- matrix(rep(0,size*length(n.alldb.g)),nrow=size,dimnames=list(s,n.alldb.g))
  ksi <- matrix(rep(0,size*length(n.ksi.g)),nrow=size,dimnames=list(s,n.ksi.g))
  for (kin in s){
    ppi[kin,ppis[[kin]]] <- 1
    ksi[kin,ksis[[kin]]] <- 1
  }
  jdist <- matrix(rep(0,size*size),nrow=size,dimnames=list(s,s))
  for (i in 1:(size-1)){
    cat(i," ")
    for (j in (i+1):size){
      M10p <- sum(ppi[i,]>ppi[j,])
      M01p <- sum(ppi[i,]<ppi[j,])
      M11p <- sum(ppi[i,]==ppi[j,])
      M10k <- sum(ksi[i,]>ksi[j,])
      M01k <- sum(ksi[i,]<ksi[j,])
      M11k <- sum(ksi[i,]==ksi[j,])
      jdist[i,j] <- alpha*(M10p+M01p)/(M10p+M01p+M11p)+(1-alpha)*(M10k+M01k)/(M10k+M01k+M11k)
    }
  }
  cat("\n")
  jdist
} # getJaccardDistances

w.net.dist <- getJaccardDistances(kin.in.net)
net.dist <- as.dist(t(w.net.dist))
save(net.dist,file="network-distance-jaccard.data")
hist(net.dist,breaks=100)
hist(net.dist^0.1,breaks=100)
hist(exp(net.dist),breaks=100)
clust <- hclust(net.dist^0.1)
plot(clust)
heatmap(as.matrix(net.dist^0.1))

classes <- unique(kinases[["Type"]])
cind <- 1:length(classes)
names(cind) <- classes
rain <- c("gray","blue","pink","cyan","orange","yellow","red","green","black","violet","salmon")
cols <- rain[cind[kinases[labels(net.dist),"Type"]]]
heatmap(as.matrix(net.dist^0.1),RowSideColors=cols,ColSideColors=cols)



