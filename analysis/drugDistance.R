setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)
library(multtest)

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
#save(full.P,file="full-transmat.data")
load("full-transmat.data")

# Kinases ---------------------------
kinases <- read.csv("UniProt_human_protein_kinases_formated-atypical-1class.txt",sep="\t",stringsAsFactors=F)
kin.class <- kinases[["Type"]]
names(kin.class) <- kinases[["UniProt_ID"]]
classes <- unique(kinases[["Type"]])
length(intersect(nodes(full.g),kinases[["UniProt_ID"]]))

# Network of kinases and proteins linked to at least n kinases by PPIs, all KSIs kept
prots <- unique(c(from.ppi,to.ppi,from.ksi,to.ksi))
counts <- rep(0,length(prots))
names(counts) <- prots
ik <- from.ppi%in%kinases[["UniProt_ID"]]
counts[to.ppi[ik]] <- counts[to.ppi[ik]]+1
ik <- to.ppi%in%kinases[["UniProt_ID"]]
counts[from.ppi[ik]] <- counts[from.ppi[ik]]+1
min.n.kinases <- 2
counts[to.ksi] <- counts[to.ksi]+min.n.kinases
counts[kinases[["UniProt_ID"]]] <- counts[kinases[["UniProt_ID"]]]+min.n.kinases
connected <- names(counts)[counts>=min.n.kinases]
kin.part <- full[full[["source"]]%in%connected & full[["target"]]%in%connected,]
kl.no.self <- as.vector(kin.part[["source"]]) != as.vector(kin.part[["target"]])
kfrom.ppi <- as.vector(kin.part[kl.no.self & kin.part[["type"]]=="PPI","source"])
kto.ppi <- as.vector(kin.part[kl.no.self & kin.part[["type"]]=="PPI","target"])
kfrom.ksi <- as.vector(kin.part[kl.no.self & (kin.part[["type"]]=="KSI" | kin.part[["type"]]=="computeKSI"),"source"])
kto.ksi <- as.vector(kin.part[kl.no.self & (kin.part[["type"]]=="KSI" | kin.part[["type"]]=="computeKSI"),"target"])
sinks <- setdiff(kto.ksi,c(kfrom.ppi,kto.ppi,kfrom.ksi))
ft <- cbind(c(kfrom.ppi,kto.ppi,kfrom.ksi,sinks,"sink"),c(kto.ppi,kfrom.ppi,kto.ksi,rep("sink",length(sinks)+1)))
kin.part.g <- ftM2graphNEL(ft, edgemode="directed")
kin.part.cc <- connectedComp(kin.part.g)
kin.part.n.nodes <- length(kin.part.cc[[1]])
kin.part.G <- subGraph(kin.part.cc[[1]],kin.part.g)
n.kin.part.G <- nodes(kin.part.G)
#kin.part.am <- graph2SparseM(kin.part.G,useweights=F)
#kin.part.P <- adjmat2transprob(kin.part.am)
#save(kin.part.P,file="kin-part-transmat.data")
load("kin-part-transmat.data")


# KEGG mapping ====================================================================
sp.kegg <- read.csv("/home/colinge/databases/uniprot_to_kegg_human-2012-08-23.csv",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c('character','character'))
kegg2names <- read.csv("/home/colinge/databases/kegg_names-2012-08-23.txt",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c("character","character"))
kegg.names <- kegg2names[[2]]
names(kegg.names) <- kegg2names[[1]]
sp.kegg <- sp.kegg[!is.na(sp.kegg[[2]]),]
for (i in 1:dim(sp.kegg)[1]){
  l <- nchar(sp.kegg[i,2])
  sp.kegg[i,2] <- paste(paste(rep("0",5-l),collapse=""),sp.kegg[i,2],sep="")
}

sp2k <- sp.kegg[sp.kegg[["AC"]] %in% n.kin.part.G,]
k.pws <- unique(sp2k[["KEGG.ID"]])


# Impact on KEGG pathways ---------------------------
ambit <- read.csv("davis2011_drugs2kinases.txt",sep="\t",check.names=FALSE,stringsAsFactors=F)
ambit <- ambit[nchar(ambit[["Target_UniProt_ID"]])>1,]
drugs <- unique(ambit[["Source_name"]])
proc <- "BH"
n.top <- 500
nperm <- 1000
kegg.res <- NULL
for (drug in drugs){
  cat(drug," ")
  d.data <- ambit[ambit[["Source_name"]]==drug,c("Target_UniProt_ID","Strength")]
  d.tar <- unique(d.data[[1]])
  d.s <- NULL
  for (t in d.tar)
    d.s <- c(d.s,min(d.data[d.data[[1]]==t,2]))
  d.s <- ambit[ambit[["Source_name"]]==drug,"Strength"]
  d.w <- 1/d.s
  d.w <- d.w/max(d.w)
  d.good <- d.tar%in%n.kin.part.G
  d.tar <- d.tar[d.good]
  d.w <- d.w[d.good]
  lx <- limitProb(P=kin.part.P,nodes=n.kin.part.G,seeds=d.tar,weights=d.w)
  p.sink <- lx["sink"]
  lx["sink"] <- 0
  lx <- lx/(1-p.sink)
  names(lx) <- n.kin.part.G
  lx.2 <- lx
  o <- order(lx.2,decreasing=T)
  thres <- lx.2[o][n.top]
  lx.2[lx.2<thres] <- 0
  
  res <- NULL
  for (pw in k.pws){
    prots <- sp2k[sp2k[["KEGG.ID"]] %in% pw,"AC"]
    if (length(prots) > 0){
      score <- sum(lx[prots])
      score.2 <- sum(lx.2[prots])
      rscores <- NULL
      for (r in 1:nperm)
        rscores <- c(rscores,sum(sample(lx,length(prots))))
      b <- var(rscores)/mean(rscores)
      c <- (mean(rscores)/sd(rscores))^2
      pval <- pgamma(score,shape=c,scale=b,lower.tail=F)
      
      rscores.2 <- NULL
      for (r in 1:nperm)
        rscores.2 <- c(rscores.2,sum(sample(lx.2,length(prots))))
      b <- var(rscores.2)/mean(rscores.2)
      c <- (mean(rscores.2)/sd(rscores.2))^2
      pval.2 <- pgamma(score.2,shape=c,scale=b,lower.tail=F)
      
      pwname <- kegg.names[pw]
      if (is.na(pwname))
        pwname <- ""
      res <- rbind(res,data.frame(drug=drug,pw=pw,descr=pwname,pwsize=length(prots),score=score,pval=pval,score.2=score.2,pval.2=pval.2,stringsAsFactors=F))
    }
  }
  rawp <- res[["pval"]]
  adj <- mt.rawp2adjp(rawp, proc)
  qval <- adj$adjp[order(adj$index),proc]
  rawp.2 <- res[["pval.2"]]
  adj.2 <- mt.rawp2adjp(rawp.2, proc)
  qval.2 <- adj.2$adjp[order(adj.2$index),proc]
  res <- cbind(res,data.frame(qval=qval,qval.2=qval.2))
  kegg.res <- rbind(kegg.res,res)
}
write.table(kegg.res,file="kegg-ambit-scores.txt",row.names=F,sep="\t")
kegg.res <- read.csv("kegg-ambit-scores.txt",sep="\t",stringsAsFactors=F,as.is=T,colClasses=c("character","character","character","integer","numeric","numeric","numeric","numeric","numeric","numeric"))

kegg.res[kegg.res[["qval.2"]]<0.001,]
kegg.res[kegg.res[["qval"]]<0.001,]

kegg.mat <- matrix(rep(1,length(drugs)*length(k.pws)),nrow=length(drugs),dimnames=list(drugs,k.pws))
for (i in 1:dim(kegg.res)[1])
  kegg.mat[kegg.res[i,"drug"],kegg.res[i,"pw"]] <- kegg.res[i,"qval"]
mq <- apply(kegg.mat,2,min)
good <- mq<0.05 & colnames(kegg.mat)!="05200"
mat <- kegg.mat[,good]
mat[mat<1E-10] <- 1E-10

pdf("figures/fig5-hm-kegg-pval-kin-part.pdf",width=12,height=8,pointsize=8,useDingbats=F)
heatmap(-log10(t(mat)),scale="none",labRow=paste(colnames(mat),kegg.names[colnames(mat)]),col=rev(heat.colors(50)))
dev.off()

kegg.mat <- matrix(rep(1,length(drugs)*length(k.pws)),nrow=length(drugs),dimnames=list(drugs,k.pws))
for (i in 1:dim(kegg.res)[1])
  kegg.mat[kegg.res[i,"drug"],kegg.res[i,"pw"]] <- kegg.res[i,"score.2"]
mq <- apply(kegg.mat,2,max)
good <- mq>quantile(kegg.mat,prob=0.95) & colnames(kegg.mat)!="05200"
mat <- kegg.mat[,good]

pdf("figures/fig5-hm-kegg-scores-kin-part.pdf",width=12,height=8,pointsize=8,useDingbats=F)
heatmap(t(mat),scale="none",labRow=paste(colnames(mat),kegg.names[colnames(mat)]),col=rev(heat.colors(50)))
dev.off()

d <- dist(sqrt(mat))
h <- hclust(d)
pdf("figures/fig-5-dendro-drug-kegg-kin-part.pdf",pointsize=6)
plot(h)
dev.off()


# Effect on proteins ===================================================

net.effect <- NULL
for (drug in drugs){
  d.data <- ambit[ambit[["Source_name"]]==drug,c("Target_UniProt_ID","Strength")]
  d.tar <- unique(d.data[[1]])
  d.s <- NULL
  for (t in d.tar)
    d.s <- c(d.s,min(d.data[d.data[[1]]==t,2]))
  d.s <- ambit[ambit[["Source_name"]]==drug,"Strength"]
  d.w <- 1/d.s
  d.w <- d.w/max(d.w)
  d.good <- d.tar%in%n.kin.part.G
  d.tar <- d.tar[d.good]
  d.w <- d.w[d.good]
  lx <- limitProb(P=kin.part.P,nodes=n.kin.part.G,seeds=d.tar,weights=d.w)
  p.sink <- lx["sink"]
  lx["sink"] <- 0
  lx <- lx/(1-p.sink)
  names(lx) <- n.kin.part.G
  net.effect <- rbind(net.effect,lx)
}
rownames(net.effect) <- drugs
thres <- apply(net.effect,1,quantile,prob=0.95)
mat <- net.effect
for (i in 1:dim(mat)[1])
  mat[i,mat[i,]<thres[i]] <- 0
no0 <- apply(mat,2,sum) > 0
mat <- mat[,no0]
sum(no0)

cl.col <- rainbow(n=length(classes))
names(cl.col) <- classes

pdf("figures/fig5-hm-prot-scores-kin-part.pdf",width=8,height=8,pointsize=6)
colors <- cl.col[kin.class[colnames(mat)]]
colors[is.na(colors)] <- "grey"
h <- heatmap(sqrt(mat),scale="none",ColSideColors=colors,col=rev(heat.colors(50)))
dev.off()

d <- dist(sqrt(mat))
h <- hclust(d)
plot(h)
pdf("figures/fig-5-dendro-drug-prot-kin-part.pdf",pointsize=6)
plot(h)
dev.off()

# kinases only
mat <- net.effect[,colnames(net.effect)%in%kinases[["UniProt_ID"]]]

cl.col <- rainbow(n=length(classes))
names(cl.col) <- classes

pdf("figures/fig5-hm-prot-scores-kin-part-kinases.pdf",width=8,height=8,pointsize=6)
colors <- cl.col[kin.class[colnames(mat)]]
colors[is.na(colors)] <- "grey"
h <- heatmap(sqrt(mat),scale="none",ColSideColors=colors,col=rev(heat.colors(50)))
dev.off()

d <- dist(sqrt(mat))
h <- hclust(d)
plot(h)
pdf("figures/fig-5-dendro-drug-prot-kin-part-kinases.pdf",pointsize=6)
plot(h)
dev.off()


# unfiltered
mat <- net.effect

cl.col <- rainbow(n=length(classes))
names(cl.col) <- classes

pdf("figures/fig5-hm-prot-scores-kin-part-unfiltered.pdf",width=8,height=8,pointsize=6)
colors <- cl.col[kin.class[colnames(mat)]]
colors[is.na(colors)] <- "grey"
h <- heatmap(sqrt(mat),scale="none",ColSideColors=colors,col=rev(heat.colors(50)))
dev.off()

d <- dist(sqrt(mat))
h <- hclust(d)
plot(h)
pdf("figures/fig-5-dendro-drug-prot-kin-part-unfiltered.pdf",pointsize=6)
plot(h)
dev.off()



# Loads GO annotations ===========================
library(GO.db)
xx <- as.list(GOTERM)
ids <- lapply(xx,GOID)
go.rel.date <- "GO-2012-08-21"
load(file=paste("/home/colinge/goa/",go.rel.date,"/sprot2GO.gen.data",sep=""))
load(file=paste("/home/colinge/goa/",go.rel.date,"/sprot2GO.split.gen.data",sep=""))
anc.bp <- as.list(GOBPANCESTOR)
anc.bp <- anc.bp[!is.na(xx)]
load(paste("/home/colinge/goa/",go.rel.date,"/all.gen.terms.data",sep=""))
all.gen.descr <- c()
for (s in all.gen.terms){
  all.gen.descr <- c(all.gen.descr, Term(xx[s][[1]]))
}
row.labels.gen <- paste(all.gen.terms,all.gen.descr)
names(row.labels.gen) <- all.gen.terms
load(paste("/home/colinge/goa/",go.rel.date,"/nr.sprot2GO.bp.data",sep=""))
load(paste("/home/colinge/goa/",go.rel.date,"/num.nr.noanc.terms.data",sep=""))

sp2g <- nr.sprot2GO.bp[intersect(n.kin.part.G,names(nr.sprot2GO.bp))]
go.bps <- setdiff(unique(unlist(sp2g)),c("all"))

proc <- "BH"
n.top <- 500
nperm <- 1000
go.res <- NULL
for (drug in drugs){
  cat(drug," ")
  d.data <- ambit[ambit[["Source_name"]]==drug,c("Target_UniProt_ID","Strength")]
  d.tar <- unique(d.data[[1]])
  d.s <- NULL
  for (t in d.tar)
    d.s <- c(d.s,min(d.data[d.data[[1]]==t,2]))
  d.s <- ambit[ambit[["Source_name"]]==drug,"Strength"]
  d.w <- 1/d.s
  d.w <- d.w/max(d.w)
  d.good <- d.tar%in%n.kin.part.G
  d.tar <- d.tar[d.good]
  d.w <- d.w[d.good]
  lx <- limitProb(P=kin.part.P,nodes=n.kin.part.G,seeds=d.tar,weights=d.w)
  p.sink <- lx["sink"]
  lx["sink"] <- 0
  lx <- lx/(1-p.sink)
  names(lx) <- n.kin.part.G
  lx.2 <- lx
  o <- order(lx.2,decreasing=T)
  thres <- lx.2[o][n.top]
  lx.2[lx.2<thres] <- 0
  
  res <- NULL
  for (bp in go.bps){
    match <- unlist(lapply(sp2g,function(x){sum(x %in% bp)}))
    prots <- names(match[match>0])
    if (length(prots) > 0){
      score <- sum(lx[prots])
      score.2 <- sum(lx.2[prots])      
      pwname <- row.labels.gen[bp]
      if (is.na(pwname))
        pwname <- ""
      res <- rbind(res,data.frame(drug=drug,pw=bp,descr=pwname,pwsize=length(prots),score=score,score.2=score.2,stringsAsFactors=F))
    }
  }
  go.res <- rbind(go.res,res)
}
write.table(go.res,file="GOBP-ambit-scores.txt",row.names=F,sep="\t")
go.res <- read.csv("GOBP-ambit-scores.txt",sep="\t",stringsAsFactors=F,as.is=T,colClasses=c("character","character","character","integer","numeric","numeric"))

go.mat <- matrix(rep(1,length(drugs)*length(go.bps)),nrow=length(drugs),dimnames=list(drugs,go.bps))
for (i in 1:dim(go.res)[1])
  go.mat[go.res[i,"drug"],go.res[i,"pw"]] <- go.res[i,"score.2"]
mq <- apply(go.mat,2,max)
good <- mq>0.2 & colnames(go.mat)!="GO:0008150"
mat <- go.mat[,good]

pdf("figures/fig5-hm-gobp-scores-kin-part.pdf",width=12,height=8,pointsize=8,useDingbats=F)
heatmap(sqrt(t(mat)),scale="none",labRow=paste(colnames(mat),row.labels.gen[colnames(mat)]),col=rev(heat.colors(50)))
dev.off()

d <- dist(sqrt(mat))
h <- hclust(d)
pdf("figures/fig-5-dendro-drug-gobp-kin-part.pdf",pointsize=6)
plot(h)
dev.off()

