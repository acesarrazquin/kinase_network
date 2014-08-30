setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)

source("/home/colinge/net-r/expandGraph.R")

# Full PPI & KSI network
full <- read.csv("full-network-20140211.txt",sep="\t",check.names=FALSE,stringsAsFactors=F)
l.no.self <- as.vector(full[["source"]]) != as.vector(full[["target"]])
from.ppi <- as.vector(full[l.no.self & full[["type"]]=="PPI","source"])
to.ppi <- as.vector(full[l.no.self & full[["type"]]=="PPI","target"])
from.ksi <- as.vector(full[l.no.self & (full[["type"]]=="KSI" | full[["type"]]=="computeKSI"),"source"])
to.ksi <- as.vector(full[l.no.self & (full[["type"]]=="KSI" | full[["type"]]=="computeKSI"),"target"])
sinks <- setdiff(to.ksi,c(from.ppi,to.ppi,from.ksi))
length(to.ppi)
length(to.ksi)
length(as.vector(full[l.no.self & full[["type"]]=="KSI","source"]))
length(as.vector(full[l.no.self & full[["type"]]=="computeKSI","source"]))
# PPI bidirectional and KSI unidirectional, a self-connected sink added for outDegree==0 KSI targets
ft <- cbind(c(from.ppi,to.ppi,from.ksi,sinks,"sink"),c(to.ppi,from.ppi,to.ksi,rep("sink",length(sinks)+1)))
full.g <- ftM2graphNEL(ft, edgemode="directed")
full.cc <- connectedComp(full.g)
full.n.nodes <- length(full.cc[[1]])
full.G <- subGraph(full.cc[[1]],full.g)
n.full.G <- nodes(full.G)
#full.am <- graph2SparseM(full.G,useweights=F)
#full.P <- adjmat2transprob(full.am)
#save(full.P,file="full-network-20140211-transmat.data")
load("full-network-20140211-transmat.data")

#shortest.G <- floyd.warshall.all.pairs.sp(full.G)
#save(shortest.G,file="full-network-20140211-shortest-path.data")
load("full-network-20140211-shortest-path.data")

# ID mapping ===============================

bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=F,sep="\t",stringsAsFactors=F)
sp2name.1 <- sp2name[sp2name[["Accession"]] %in% n.full.G,]
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
kin.in.net <- intersect(kinases[["UniProt_ID"]],n.full.G)
rownames(kinases) <- kinases[["UniProt_ID"]]

# Distance computation ----------------

getDistances <- function(nodes,min.connect=5){
  targets <- adj(full.G,nodes)
  num.targets <- unlist(lapply(targets,length))
  s <- nodes[num.targets>=min.connect]
  size <- length(s)
  w.net.dist <- matrix(rep(0,size*size),nrow=size,dimnames=list(s,s))
  for (i in 1:(size-1)){
    cat(i," ")
    a <- s[i]
    a.targets <- targets[[a]]
    for (j in (i+1):size){
      b <- s[j]
      b.targets <- targets[[b]]
      a.spe <- setdiff(a.targets,b.targets)
      b.spe <- setdiff(b.targets,a.targets)
      tot <- 0
      wn.dist <- 0
      for (t in a.spe){
        mt <- min(shortest.G[t,b.targets])
        if (!is.infinite(mt)){
          wn.dist <- wn.dist+mt
          tot <- tot+1
        }
        else{
          #cat(i,j,t,"\n")
        }
      }
      for (t in b.spe){
        tot <- tot+1
        mt <- min(shortest.G[t,a.targets])
        if (!is.infinite(mt)){
          wn.dist <- wn.dist+mt
          tot <- tot+1
        }
        else{
          #cat(i,j,t,"\n")
        }
      }
      m <- length(a.spe)+length(b.spe)
      n <- length(intersect(union(a.targets,b.targets),n.full.G))
      if (tot+m>=min.connect)
        w.net.dist[a,b] <- wn.dist/tot*m/n
      else
        w.net.dist[a,b] <- NA
    }
  }
  cat("\n")
  w.net.dist
} # getDistances


w.net.dist <- getDistances(kin.in.net)
net.dist <- as.dist(t(w.net.dist))
save(net.dist,file="network-distance.data")

load("network-distance.data")

clust <- hclust(net.dist,method="ward")
plot(clust)
heatmap(as.matrix(net.dist))
hist(net.dist)

classes <- unique(kinases[["Type"]])
# cind <- 1:length(classes)
# names(cind) <- classes
# rain <- rainbow(length(classes))
#cols <- rain[cind[kinases[labels(net.dist),"Type"]]]
rain <- c(rgb(0,1,1,1),rgb(1,0.6,0,1),rgb(0,0,1,1),rgb(1,0.6,0.6,1),rgb(0,1,0,1),
          rgb(0.6,0.8,0,1),rgb(0.8,0.8,0.8,1),rgb(1,1,0,1),rgb(1,0,1,1),rgb(1,0,0,1),rgb(0.4,0.4,1,1))
names(rain) <- c("AGC","Atypical","CAMK","CK1","CMGC","NEK","Other","STE","TKL","Tyr","RGC")
cols <- rain[kinases[labels(net.dist),"Type"]]
lab <- sp.names[labels(net.dist)]
hcols <- heat.colors(50)
heatmap(as.matrix(net.dist),RowSideColors=cols,ColSideColors=cols,col=hcols)

pdf("figures/kinase-network-distance.pdf",width=20,height=16,useDingbats=F,pointsize=6)
hm <- heatmap(as.matrix(net.dist),RowSideColors=cols,ColSideColors=cols,scale="none",keep.dendro=T,labRow=lab,labCol=lab,col=hcols)
dev.off()
