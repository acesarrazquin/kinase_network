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
alldb.cc <- connectedComp(alldb.g)
alldb.n.nodes <- length(alldb.cc[[1]])
alldb.G <- subGraph(alldb.cc[[1]],alldb.g)
n.alldb.G <- nodes(alldb.G)

# ID mapping ===============================

bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=F,sep="\t",stringsAsFactors=F)
sp2name.1 <- sp2name[sp2name[["Accession"]] %in% n.alldb.G,]
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

# Mann's 11 cell line core proteomes ------------------
core <- read.csv("/home/colinge/haplogene/core-prots.csv",sep="\t",stringsAsFactors=F)
present <- NULL
names <- NULL
for (j in seq(11,41,by=3)){
  pres <- (core[[j]] > 0) | (core[[j+1]] > 0) | (core[[j+2]] > 0)
  present <- cbind(present,pres)
  names <- c(names,strsplit(names(core)[j],split="[._]",perl=T)[[1]][2])
}
colnames(present) <- names
present[is.na(present)] <- FALSE
num.cl <- apply(present,1,sum)
sp.acs <- read.csv("/home/colinge/haplogene/human-sp-acs.txt",sep="\t",stringsAsFactors=F)[[1]]
acs <- strsplit(core[["Uniprot"]]," ")
acs.sp <- lapply(acs,function(x){intersect(sp.acs,gsub("-.+","",x))[1]})
num.names <- unlist(acs.sp)
good <- !is.na(num.names)
common <- num.cl[good]
names(common) <- num.names[good]

# Human kinases --------------------------------------------
kinases <- read.csv("UniProt_human_protein_kinases_formated-atypical-1class.txt",sep="\t",stringsAsFactors=F)
kin.in.net <- intersect(kinases[["UniProt_ID"]],n.alldb.G)
rownames(kinases) <- kinases[["UniProt_ID"]]
high.prot <- intersect(n.alldb.G,names(common)[common==11])

# KSI network ----------------------
ksi <- read.table("ksi-20140211.txt",sep="\t",stringsAsFactors=F,header=T)
k.no.self <- as.vector(ksi[["source"]]) != as.vector(ksi[["target"]])
from <- as.vector(ksi[["source"]][k.no.self])
to <- as.vector(ksi[["target"]][k.no.self])
ft <- cbind(from,to)
ksi.g <- ftM2graphNEL(ft, edgemode="directed")
n.ksi.g <- nodes(ksi.g)
ksi.cc <- connectedComp(ksi.g)
ksi.n.nodes <- length(ksi.cc[[1]])
ksi.G <- subGraph(ksi.cc[[1]],ksi.g)
n.ksi.G <- nodes(ksi.G)


# Global analysis =========================================================================

# PPI network topological analysis ----------
ppi.deg <- degree(alldb.G)
h.ppi <- hist(log10(ppi.deg),breaks=100)
pos.p <- h.ppi$counts>0
plot(x=h.ppi$mid[pos.p],y=log10(h.ppi$counts[pos.p]))

# KSI network topological analysis ---------------
ksi.deg <- degree(ksi.G)
h.ksi.o <- hist(log10(ksi.deg$outDegree),breaks=50)
h.ksi.i <- hist(log10(ksi.deg$inDegree),breaks=50)
pos.o <- h.ksi.o$counts>0
pos.i <- h.ksi.i$counts>0
plot(x=h.ksi.o$mid[pos.o],y=log10(h.ksi.o$counts[pos.o]),ylim=c(0,3))
lines(x=h.ksi.i$mid[pos.i],y=log10(h.ksi.i$counts[pos.i]),type="p",col="blue")


# Kinase-kinase analysis ==================================================================

# PPI network
is.kinase <- alldb[["source"]]%in%kinases[["UniProt_ID"]] & alldb[["target"]]%in%kinases[["UniProt_ID"]]
from <- as.vector(alldb[["source"]][l.no.self & is.kinase])
to <- as.vector(alldb[["target"]][l.no.self & is.kinase])
ft <- cbind(from,to)
kin.ppi.g <- ugraph(ftM2graphNEL(ft, edgemode="directed"))
kin.ppi.cc <- connectedComp(kin.ppi.g)
kin.ppi.G <- subGraph(kin.ppi.cc[[1]],kin.ppi.g)
n.kin.ppi.G <- nodes(kin.ppi.G)

sum((alldb[["source"]]%in%kinases[["UniProt_ID"]] | alldb[["target"]]%in%kinases[["UniProt_ID"]]) & l.no.self)

kin.ppi.deg <- degree(kin.ppi.G)
h.ppi <- hist(log10(kin.ppi.deg),breaks=100)
pos.p <- h.ppi$counts>0
plot(x=h.ppi$mid[pos.p],y=log10(h.ppi$counts[pos.p]))

top <- kin.ppi.deg[order(kin.ppi.deg,decreasing=T)][20]
collab <- kin.ppi.deg>=top
sp.names[n.kin.ppi.G[collab]]
partner.g <- subGraph(n.kin.ppi.G[collab],kin.ppi.G)
plot(partner.g)

writeDirectedNetwork(ctrl.g,"figures/controling",kin.ksi.deg,ksi.deg)


# KSI network
is.kinase <- ksi[["source"]]%in%kinases[["UniProt_ID"]] & ksi[["target"]]%in%kinases[["UniProt_ID"]]
from <- as.vector(ksi[["source"]][k.no.self & is.kinase])
to <- as.vector(ksi[["target"]][k.no.self & is.kinase])
ft <- cbind(from,to)
kin.ksi.g <- ftM2graphNEL(ft, edgemode="directed")
kin.ksi.cc <- connectedComp(kin.ksi.g)
kin.ksi.G <- subGraph(kin.ksi.cc[[1]],kin.ksi.g)
n.kin.ksi.G <- nodes(kin.ksi.G)

kin.ksi.deg <- degree(kin.ksi.G)
h.ksi.io <- hist(log10(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree),breaks=50)
h.ksi.o <- hist(log10(kin.ksi.deg$outDegree),breaks=50)
h.ksi.i <- hist(log10(kin.ksi.deg$inDegree),breaks=50)
pos.io <- h.ksi.io$counts>0
pos.o <- h.ksi.o$counts>0
pos.i <- h.ksi.i$counts>0
pdf("figures/fig-4-power-law.pdf",width=2,height=2,pointsize=8,useDingbats=F)
plot(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io]),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
symbols(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io]),
        circles=rep(0.03,sum(pos.io)),inches=F,bg="black",fg=NULL,xlim=c(0,3.5),add=T)
symbols(x=h.ksi.o$mid[pos.o],y=log10(h.ksi.o$counts[pos.o]),
        circles=rep(0.03,sum(pos.o)),inches=F,bg="red",fg=NULL,xlim=c(0,3.5),add=T)
symbols(x=h.ksi.i$mid[pos.i],y=log10(h.ksi.i$counts[pos.i]),
        circles=rep(0.03,sum(pos.i)),inches=F,bg="blue",fg=NULL,xlim=c(0,3.5),add=T)
reg.io <- lm(y~x,data.frame(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io])))
abline(a=reg.io$coefficients[1],b=reg.io$coefficients[2],col="black",lwd=1)
reg.o <- lm(y~x,data.frame(x=h.ksi.o$mid[pos.o],y=log10(h.ksi.o$counts[pos.o])))
abline(a=reg.o$coefficients[1],b=reg.o$coefficients[2],col="red",lwd=1)
reg.i <- lm(y~x,data.frame(x=h.ksi.i$mid[pos.i],y=log10(h.ksi.i$counts[pos.i])))
abline(a=reg.i$coefficients[1],b=reg.i$coefficients[2],col="blue",lwd=1)
dev.off()




plot(y=kin.ksi.deg$outDegree,x=kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree))
plot(y=kin.ksi.deg$outDegree,x=kin.ksi.deg$inDegree)

plot(x=kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree),y=ksi.deg$outDegree[n.kin.ksi.G]/(ksi.deg$outDegree[n.kin.ksi.G]+ksi.deg$inDegree[n.kin.ksi.G]))

plot(x=jitter(kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)),y=kin.ksi.deg$outDegree/ksi.deg$outDegree[n.kin.ksi.G])
abline(h=0.12,col="red")

controlers <- (kin.ksi.deg$outDegree>=10 & kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.8) |
  (kin.ksi.deg$outDegree>=4 & kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)==1)
sp.names[n.kin.ksi.G[controlers]]
ctrl.g <- subGraph(n.kin.ksi.G[controlers],kin.ksi.G)
plot(ctrl.g)

writeDirectedNetwork(ctrl.g,"figures/controling3",kin.ksi.deg,ksi.deg)

# controlers <- kin.ksi.deg$outDegree>=10 & kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.8 & kin.ksi.deg$outDegree>0.12*ksi.deg$outDegree[n.kin.ksi.G]
# sp.names[n.kin.ksi.G[controlers]]
# ctrl.g <- subGraph(n.kin.ksi.G[controlers],kin.ksi.G)
# plot(ctrl.g)
# 
# controlers <- kin.ksi.deg$outDegree>0.25*ksi.deg$outDegree[n.kin.ksi.G]
# sp.names[n.kin.ksi.G[controlers]]
# ctrl.g <- subGraph(n.kin.ksi.G[controlers],kin.ksi.G)
# plot(ctrl.g)
# 

plot(y=kin.ksi.deg$inDegree,x=kin.ksi.deg$inDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree))
plot(y=kin.ksi.deg$inDegree,x=kin.ksi.deg$outDegree)

controled <- (kin.ksi.deg$inDegree>=5 & kin.ksi.deg$inDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.75)
sp.names[n.kin.ksi.G[controled]]
ctrled.g <- subGraph(n.kin.ksi.G[controled],kin.ksi.G)


plot(ctrled.g)
writeDirectedNetwork(ctrled.g,"figures/controled",kin.ksi.deg,ksi.deg)




writeDirectedNetwork <- function(g,base.name,deg,deg2){
  edge.file <- paste(base.name,"-edges.txt",sep="")
  e <- edges(g)
  cat("source\ttarget\tedge.type\n",sep="",file=edge.file)
  for (i in 1:length(e)){
    source <- names(e)[i]
    for (target in e[[i]])
      cat(source,"\t",target,"\tPPI\n",sep="",file=edge.file,append=T)
  }
  node.file <- paste(base.name,"-nodes.txt",sep="")
  n <- nodes(g)
  write.table(data.frame(node=n,name=sp.names[n],kin.type=kinases[n,"Type"],out.deg=deg$outDegree[n],in.deg=deg$inDegree[n],
                         deg.ratio=deg$outDegree[n]/(deg$outDegree[n]+deg$inDegree[n]),
                         full.out.deg=deg2$outDegree[n],stringsAsFactors=F),
              quote=F,sep="\t",row.names=F,file=node.file)
}


# Hierarchical topology test ---------------------
# Exclude kinases with 0 or 1 substrates only!
good <- ksi.deg$outDegree[n.kin.ksi.G]>=5
kin.ksi.g2 <- subGraph(n.kin.ksi.G[good],kin.ksi.G)
kin.ksi.cc2 <- connectedComp(kin.ksi.g2)
kin.ksi.G2 <- subGraph(kin.ksi.cc2[[1]],kin.ksi.g2)
n.kin.ksi.G2 <- nodes(kin.ksi.G2)
kin.ksi.deg2 <- degree(kin.ksi.G2)


deg.ratios <- kin.ksi.deg2$outDegree/(kin.ksi.deg2$outDegree+kin.ksi.deg2$inDegree)
ratios <- deg.ratios[order(deg.ratios,decreasing=T)]
linked.ratios <- vapply(names(ratios),function(x){sum(ratios[unlist(adj(kin.ksi.G2,x))]<ratios[x])},0)
plot(linked.ratios)
sp.names[names(linked.ratios)[linked.ratios>35]]
score <- sum(linked.ratios)

r.scores <- NULL
for (r in 1:200){
  rg <- randomNodeGraph(kin.ksi.deg2$outDegree+kin.ksi.deg2$inDegree)
  e <- edges(rg)
  for (n in nodes(rg))
    if (n %in% e[[n]])
      removeEdge(n,n,rg)
  r.deg <- degree(rg)
  r.deg.ratios <- r.deg$outDegree/(r.deg$outDegree+r.deg$inDegree)
  r.ratios <- r.deg.ratios[order(r.deg.ratios,decreasing=T)]
  r.linked.ratios <- vapply(names(r.ratios),function(x){sum(ratios[unlist(adj(rg,x))]<r.ratios[x])},0)
  r.scores <- c(r.scores,sum(r.linked.ratios)/numEdges(rg)*numEdges(kin.ksi.G))  
}
write.table(r.scores,file="hierarchical-r-scores-ksi.txt",row.names=F,col.names=F)
r.scores <- read.table("hierarchical-r-scores-ksi.txt")[[1]]
score/mean(r.scores)

pdf("figures/fig-4-hierarchical-ksi.pdf",width=2.5,height=2,pointsize=8,useDingbats=F)
plot(density(r.scores),xlim=c(min(c(score,r.scores)),max(c(score,r.scores))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(r.scores))
abline(v=score,col="red")
dev.off()


rg <- randomNodeGraph(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)
e <- edges(rg)
for (n in nodes(rg))
  if (n %in% e[[n]])
    removeEdge(n,n,rg)
r.deg <- degree(rg)
h.ksi.io <- hist(log10(r.deg$outDegree+r.deg$inDegree),breaks=50)
h.ksi.o <- hist(log10(r.deg$outDegree),breaks=50)
h.ksi.i <- hist(log10(r.deg$inDegree),breaks=50)
pos.io <- h.ksi.io$counts>0
pos.o <- h.ksi.o$counts>0
pos.i <- h.ksi.i$counts>0
pdf("figures/fig-4-power-law-random-true.pdf",width=2,height=2,pointsize=8,useDingbats=F)
plot(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io]),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
symbols(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io]),
        circles=rep(0.03,sum(pos.io)),inches=F,bg="black",fg=NULL,xlim=c(0,3.5),add=T)
symbols(x=h.ksi.o$mid[pos.o],y=log10(h.ksi.o$counts[pos.o]),
        circles=rep(0.03,sum(pos.o)),inches=F,bg="red",fg=NULL,xlim=c(0,3.5),add=T)
symbols(x=h.ksi.i$mid[pos.i],y=log10(h.ksi.i$counts[pos.i]),
        circles=rep(0.03,sum(pos.i)),inches=F,bg="blue",fg=NULL,xlim=c(0,3.5),add=T)
reg.io <- lm(y~x,data.frame(x=h.ksi.io$mid[pos.io],y=log10(h.ksi.io$counts[pos.io])))
abline(a=reg.io$coefficients[1],b=reg.io$coefficients[2],col="black",lwd=1)
reg.o <- lm(y~x,data.frame(x=h.ksi.o$mid[pos.o],y=log10(h.ksi.o$counts[pos.o])))
abline(a=reg.o$coefficients[1],b=reg.o$coefficients[2],col="red",lwd=1)
reg.i <- lm(y~x,data.frame(x=h.ksi.i$mid[pos.i],y=log10(h.ksi.i$counts[pos.i])))
abline(a=reg.i$coefficients[1],b=reg.i$coefficients[2],col="blue",lwd=1)
dev.off()


# Experimental KSI only
is.kinase <- ksi[["source"]]%in%kinases[["UniProt_ID"]] & ksi[["target"]]%in%kinases[["UniProt_ID"]] & ksi[["type"]]=="KSI"
from <- as.vector(ksi[["source"]][k.no.self & is.kinase])
to <- as.vector(ksi[["target"]][k.no.self & is.kinase])
ft <- cbind(from,to)
kin.eksi.g <- ftM2graphNEL(ft, edgemode="directed")
kin.eksi.cc <- connectedComp(kin.eksi.g)
kin.eksi.G <- subGraph(kin.eksi.cc[[1]],kin.eksi.g)
n.kin.eksi.G <- nodes(kin.eksi.G)
good <- ksi.deg$outDegree[n.kin.eksi.G]>=5
kin.eksi.g <- subGraph(n.kin.eksi.G[good],kin.eksi.G)
kin.eksi.cc <- connectedComp(kin.eksi.g)
kin.eksi.G <- subGraph(kin.eksi.cc[[1]],kin.eksi.g)
n.kin.eksi.G <- nodes(kin.eksi.G)
kin.eksi.deg <- degree(kin.eksi.G)

deg.ratios <- kin.eksi.deg$outDegree/(kin.eksi.deg$outDegree+kin.eksi.deg$inDegree)
ratios <- deg.ratios[order(deg.ratios,decreasing=T)]
linked.ratios <- vapply(names(ratios),function(x){sum(ratios[unlist(adj(kin.eksi.G,x))]<ratios[x])},0)
plot(y=linked.ratios,x=jitter(ratios))
plot(linked.ratios)
e.ksi.score <- sum(linked.ratios)

e.r.scores <- NULL
for (r in 1:200){
  rg <- randomNodeGraph(kin.eksi.deg$outDegree+kin.eksi.deg$inDegree)
  e <- edges(rg)
  for (n in nodes(rg))
    if (n %in% e[[n]])
      removeEdge(n,n,rg)
  r.deg <- degree(rg)
  r.deg.ratios <- r.deg$outDegree/(r.deg$outDegree+r.deg$inDegree)
  r.ratios <- r.deg.ratios[order(r.deg.ratios,decreasing=T)]
  r.linked.ratios <- vapply(names(r.ratios),function(x){sum(ratios[unlist(adj(rg,x))]<r.ratios[x])},0)
  e.r.scores <- c(e.r.scores,sum(r.linked.ratios)/numEdges(rg)*numEdges(kin.eksi.G))  
}
write.table(e.r.scores,file="hierarchical-r-scores-exp-ksi.txt",quote=F,row.names=F)
e.r.scores <- read.table("hierarchical-r-scores-exp-ksi.txt",header=T)[[1]]
e.ksi.score/mean(e.r.scores)

pdf("figures/fig-4-hierarchical-exp-ksi.pdf",width=2.5,height=2,pointsize=8,useDingbats=F)
plot(density(e.r.scores),xlim=c(min(c(e.ksi.score,e.r.scores)),max(c(e.ksi.score,e.r.scores))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(e.r.scores))
abline(v=e.ksi.score,col="red")
dev.off()


# Exact random model (not really random any more...)
e.r.scores <- NULL
for (r in 1:20){
  cat(r," ")
  rg <- kin.eksi.G
  rr <- 0
  while (rr < 400){
    ed <- edges(rg)
    src <- sample(n.kin.eksi.G,2)
    ts1 <- ed[[src[1]]]
    ts2 <- ed[[src[2]]]
    if ((length(ts1)>0) && (length(ts2)>1)){
      t1 <- sample(ts1,1)
      t2 <- sample(setdiff(ts2,t1),1)
      rg <- removeEdge(src[1],t1,rg)
      rg <- removeEdge(src[2],t2,rg)
      rg <- addEdge(src[1],t2,rg)
      rg <- addEdge(src[2],t1,rg)
      #cat(src[1],t1,src[2],t2,"\n")
      rr <- rr+1
    }    
  }
  r.deg <- degree(rg)
  r.deg.ratios <- r.deg$outDegree/(r.deg$outDegree+r.deg$inDegree)
  r.ratios <- r.deg.ratios[order(r.deg.ratios,decreasing=T)]
  r.linked.ratios <- vapply(names(r.ratios),function(x){sum(ratios[unlist(adj(rg,x))]<r.ratios[x])},0)
  e.r.scores <- c(e.r.scores,sum(r.linked.ratios)/numEdges(rg)*numEdges(kin.eksi.G))  
}

# PPIs
deg.ratios <- kin.ppi.deg
ratios <- deg.ratios[order(deg.ratios,decreasing=T)]
linked.ratios <- vapply(names(ratios),function(x){sum(ratios[unlist(adj(kin.ppi.G,x))]<ratios[x])},0)
plot(linked.ratios)
p.score <- sum(linked.ratios)

p.r.scores <- NULL
for (r in 1:200){
  rg <- ugraph(randomNodeGraph(kin.ppi.deg))
  e <- edges(rg)
  for (n in nodes(rg))
    if (n %in% e[[n]])
      removeEdge(n,n,rg)
  r.deg <- degree(rg)
  r.deg.ratios <- r.deg
  r.ratios <- r.deg.ratios[order(r.deg.ratios,decreasing=T)]
  r.linked.ratios <- vapply(names(r.ratios),function(x){sum(ratios[unlist(adj(rg,x))]<r.ratios[x])},0)
  p.r.scores <- c(p.r.scores,sum(r.linked.ratios)/numEdges(rg)*numEdges(kin.ppi.G))  
}
write.table(p.r.scores,file="hierarchical-r-scores-ppi.txt",row.names=F,quote=F)
p.r.scores <- read.table("hierarchical-r-scores-ppi.txt")[[1]]
p.score/mean(p.r.scores)

pdf("figures/fig-4-hierarchical-ppi.pdf",width=2.5,height=2,pointsize=8,useDingbats=F)
plot(density(p.r.scores),xlim=c(min(c(p.score,p.r.scores)),max(c(p.score,p.r.scores))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(p.r.scores))
abline(v=p.score,col="red")
dev.off()
sum(p.r.scores>=p.score)/length(p.r.scores) # 0.02

# Full PPI network
deg.ratios <- degree(alldb.G)
ratios <- deg.ratios[order(deg.ratios,decreasing=T)]
linked.ratios <- vapply(names(ratios),function(x){sum(ratios[unlist(adj(alldb.G,x))]<ratios[x])},0)
plot(linked.ratios)
a.score <- sum(linked.ratios)
alldb.deg <- degree(alldb.G)

a.r.scores <- NULL
for (r in 1:10){
  rg <- ugraph(randomNodeGraph(alldb.deg))
  e <- edges(rg)
  for (n in nodes(rg))
    if (n %in% e[[n]])
      removeEdge(n,n,rg)
  r.deg <- degree(rg)
  r.deg.ratios <- r.deg
  r.ratios <- r.deg.ratios[order(r.deg.ratios,decreasing=T)]
  r.linked.ratios <- vapply(names(r.ratios),function(x){sum(ratios[unlist(adj(rg,x))]<r.ratios[x])},0)
  a.r.scores <- c(a.r.scores,sum(r.linked.ratios)/numEdges(rg)*numEdges(alldb.G))  
}
write.table(a.r.scores,file="hierarchical-r-scores-all-ppi.txt",row.names=F,quote=F)
a.r.scores <- read.table("hierarchical-r-scores-all-ppi.txt")[[1]]

a.r.scores <- c(103177.9,103245.2,103112.5,103154.7,103179.4,103138.3,103123.8,103276.8,103102.0,103149.5,103184.0,103221.2,103247.1,103196.3,
               103193.0,103258.2,103171.7,103162.0,103206.4,103190.2,103156.7,103084.8,103301.6,103146.3,103183.8,103238.0,103162.1,103184.8)
a.score <- 103737

a.score/mean(a.r.scores)

pdf("figures/fig-4-hierarchical-all-ppi.pdf",width=2.5,height=2,pointsize=8,useDingbats=F)
plot(density(a.r.scores),xlim=c(min(c(a.score,a.r.scores)),max(c(a.score,a.r.scores))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(a.r.scores))
abline(v=a.score,col="red")
dev.off()
sum(a.r.scores>=a.score)/length(a.r.scores)

# ==========================================================================
# Selection of specific node sets ==========================================

# Controlers
#kin.controlers.i <- (kin.ksi.deg$outDegree>=10 & kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.8) 
kin.controlers.i <- (kin.ksi.deg$outDegree>=10 & kin.ksi.deg$outDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.8) | (kin.ksi.deg$outDegree>=5 & kin.ksi.deg$inDegree==0)
kin.controlers <- sp.names[n.kin.ksi.G[kin.controlers.i]]
write.table(kin.controlers,file="prot-lists/kin-controlers.txt",quote=F,row.names=F,col.names=F)

controlers.i <- (ksi.deg$outDegree>=10 & ksi.deg$outDegree/(ksi.deg$outDegree+ksi.deg$inDegree)>0.8) 
controlers <- sp.names[n.ksi.G[controlers.i]]
write.table(controlers,file="prot-lists/controlers.txt",quote=F,row.names=F,col.names=F)

setdiff(kin.controlers,controlers)
setdiff(controlers,kin.controlers)

kin.controled.i <- (kin.ksi.deg$inDegree>=5 & kin.ksi.deg$inDegree/(kin.ksi.deg$outDegree+kin.ksi.deg$inDegree)>0.8)
kin.controled <- sp.names[n.kin.ksi.G[kin.controled.i]]
write.table(kin.controled,file="prot-lists/kin-controled.txt",quote=F,row.names=F,col.names=F)

controled.i <- (ksi.deg$inDegree>=10 & ksi.deg$inDegree/(ksi.deg$outDegree+ksi.deg$inDegree)>0.8)
controled <- sp.names[n.ksi.G[controled.i]]
write.table(controled,file="prot-lists/controled.txt",quote=F,row.names=F,col.names=F)
deg.controled <- ksi.deg$inDegree[controled.i]
hist(deg.controled)
top.controled <- sp.names[names(deg.controled[deg.controled>quantile(deg.controled,prob=0.75)])]
write.table(top.controled,file="prot-lists/top-controled.txt",quote=F,row.names=F,col.names=F)
setdiff(top.controled,sp.names[kinases[["Name"]]])

kdg <- kin.ksi.deg$outDegree+kin.ksi.deg$inDegree
integrators.i <- kdg>=10 & kin.ksi.deg$outDegree/kdg>0.4 & kin.ksi.deg$outDegree/kdg<0.6
integrators <- sp.names[n.kin.ksi.G[integrators.i]]
write.table(integrators,file="prot-lists/integrators.txt",quote=F,row.names=F,col.names=F)

cdg <- ppi.deg[intersect(n.alldb.G,n.ksi.G[controlers.i])]
hist(cdg)
bound.controlers <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.controlers,file="prot-lists/bound-controlers.txt",quote=F,row.names=F,col.names=F)

cdg <- kin.ppi.deg[intersect(n.kin.ppi.G,n.kin.ksi.G[kin.controlers.i])]
hist(cdg)
bound.kin.controlers <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.kin.controlers,file="prot-lists/bound-kin.controlers.txt",quote=F,row.names=F,col.names=F)

cdg <- ppi.deg[intersect(n.alldb.G,n.kin.ksi.G[kin.controlers.i])]
hist(cdg)
bound.kin.controlers2 <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.kin.controlers2,file="prot-lists/bound-kin-controlers2.txt",quote=F,row.names=F,col.names=F)

cdg <- kin.ppi.deg[intersect(n.kin.ppi.G,n.kin.ksi.G[integrators.i])]
hist(cdg)
bound.integrators <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.integrators,file="prot-lists/bound-integrators.txt",quote=F,row.names=F,col.names=F)

cdg <- ppi.deg[intersect(n.alldb.G,n.kin.ksi.G[integrators.i])]
hist(cdg)
bound.integrators2 <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.integrators2,file="prot-lists/bound-integrators2.txt",quote=F,row.names=F,col.names=F)

cdg <- kin.ppi.deg[intersect(n.kin.ppi.G,n.kin.ksi.G[kin.controled.i])]
hist(cdg)
bound.kin.controled <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.kin.controled,file="prot-lists/bound-kin-controled.txt",quote=F,row.names=F,col.names=F)

cdg <- ppi.deg[intersect(n.alldb.G,n.ksi.G[controled.i])]
hist(cdg)
bound.controled <- sp.names[names(cdg)[cdg>quantile(cdg,prob=0.75)]]
write.table(bound.controled,file="prot-lists/bound-controled.txt",quote=F,row.names=F,col.names=F)

hubs <- names(ppi.deg[ppi.deg>quantile(ppi.deg,prob=0.99)])
no.control.hubs <- setdiff(sp.names[hubs],controled)

plot(density(ksi.deg$outDegree[n.ksi.G[controlers.i]]+ksi.deg$inDegree[n.ksi.G[controlers.i]]))
lines(density(ksi.deg$outDegree[n.ksi.G[integrators.i]]+ksi.deg$inDegree[n.ksi.G[integrators.i]]),col="blue")
lines(density(ksi.deg$outDegree[n.ksi.G[controled.i]]+ksi.deg$inDegree[n.ksi.G[controled.i]]),col="green")

plot(density(ppi.deg[intersect(n.alldb.G,n.ksi.G[controlers.i])]))
lines(density(ppi.deg[intersect(n.alldb.G,n.ksi.G[integrators.i])]),col="blue")
lines(density(ppi.deg[intersect(n.alldb.G,n.ksi.G[controled.i])]),col="green")

strong.targets <- read.csv("strong-dpi.txt",sep="\t",stringsAsFactors=F)[[1]]
sp.names[intersect(n.kin.ksi.G[kin.controlers.i],strong.targets)]
sp.names[setdiff(n.kin.ksi.G[kin.controlers.i],strong.targets)]

sp.names[intersect(n.ksi.G[controlers.i],strong.targets)]
sp.names[setdiff(n.ksi.G[controlers.i],strong.targets)]

sp.names[intersect(n.kin.ksi.G[integrators.i],strong.targets)]
sp.names[setdiff(n.kin.ksi.G[integrators.i],strong.targets)]

sp.names[intersect(n.kin.ksi.G[kin.controled.i],strong.targets)]
sp.names[setdiff(n.kin.ksi.G[kin.controled.i],strong.targets)]

# Cancer
census <- read.csv("diseases//cancer_gene_census.txt",stringsAsFactors=F,check.names=F,sep="\t")
in.cen <- kin.names[["Kinase"]] %in% census[[1]]
sp.kegg <- read.csv("/home/colinge/databases/uniprot_to_kegg_human-2012-08-23.csv",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c('character','character'))
kegg2names <- read.csv("/home/colinge/databases/kegg_names-2012-08-23.txt",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c("character","character"))
kegg.names <- kegg2names[[2]]
names(kegg.names) <- kegg2names[[1]]
kegg.cancer <- NULL
for (term in unique(sp.kegg[[2]]))
  if (length(intersect(c('cancer','leukemia','carcinoma','Glioma','Melanoma'),unlist(strsplit(kegg.names[term],split=" ")))) >= 1){
    cat(term,kegg.names[term],"\n")
    mem <- sp.kegg[sp.kegg[[2]]==term,1]
    #inpw <- intersect(kin.names[["AC"]],mem)    
    kegg.cancer <- c(kegg.cancer,mem)
  }
kegg.cancer <- unique(kegg.cancer)
kin.cancer <- union(kinases[kinases[["Name"]]%in%census[[1]],"UniProt_ID"],intersect(kinases[["UniProt_ID"]],kegg.cancer))

sp.names[intersect(n.kin.ksi.G[kin.controlers.i],kin.cancer)]
sp.names[setdiff(n.kin.ksi.G[kin.controlers.i],kin.cancer)]

sp.names[intersect(n.ksi.G[controlers.i],kin.cancer)]
sp.names[setdiff(n.ksi.G[controlers.i],kin.cancer)]

sp.names[intersect(n.kin.ksi.G[integrators.i],kin.cancer)]
sp.names[setdiff(n.kin.ksi.G[integrators.i],kin.cancer)]

sp.names[intersect(n.kin.ksi.G[kin.controled.i],kin.cancer)]
sp.names[setdiff(n.kin.ksi.G[kin.controled.i],kin.cancer)]

