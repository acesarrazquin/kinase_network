setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)

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
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2014-01-30.txt",check.names=F,sep="\t",stringsAsFactors=F)
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

# Rix, Colinge, PLoS ONE, PMID 24130846, 2013 Oct =========================================
drugModelALL <- function(drug,sc.table,CL,cid=NA,plot=F){
  score.col <- paste(drug," score",sep="")
  good <- sc.table[,score.col] > 0
  targets <- sc.table[good,"AC"]
  weights <- sc.table[good,score.col]
  strong <- rep('no',length(weights))
  if (plot)
    boxplot(weights,main=drug)
  bs <- boxplot.stats(weights)
  strong[weights%in%bs$out] <- 'yes'
  strong[weights>=0.5*max(weights)] <- 'yes'
  len <- length(targets)
  data.frame(source=rep(drug,len),target=targets,source_name=rep(drug,len),target_name=sp.names[targets],PMIDs=rep('24130846',len),dates=rep(2013,len),
             sources=rep('24130846',len),type=rep('DPI',len),source_IDs=rep(cid,len),target_form=rep("",len),strong_target=strong,strength=weights,strength_type=rep('empirical.score',len),strength_direction=rep("high",len),
             technology=rep('chemprot',len),cell_type=rep(CL,len),stringsAsFactors=F)
}

bv.173 <- read.csv("/home/colinge/CP22-ALL/CP22_Scores_BV173_for JC.csv",check.names=F,sep="\t",stringsAsFactors=F)
bv.173[bv.173[["AC"]]=="bcrabl","AC"] <- bcrabl
bv.173[["AC"]] <- sapply(strsplit(bv.173[["AC"]],split="-"),function(x) x[1])
cp22 <- drugModelALL("Dasatinib",bv.173,"BV-173","CID:3062316,CHEMBL1421,CHEBI:49375,DB01254",plot=T)
cp22 <- rbind(cp22,drugModelALL("Bosutinib",bv.173,"BV-173","CID:5328940,CHEMBL288441,CHEBI:39112",plot=T))
cp22 <- rbind(cp22,drugModelALL("Bafetinib",bv.173,"BV-173","CID:24853523",plot=T))
cp22 <- rbind(cp22,drugModelALL("Nilotinib",bv.173,"BV-173","CID:644241,CHEBI:52172",plot=T))

z.119 <- read.csv("/home/colinge/CP22-ALL/CP22_Scores_Z119_for JC.csv",check.names=F,sep="\t",stringsAsFactors=F)
z.119[z.119[["AC"]]=="bcrabl","AC"] <- bcrabl
z.119[["AC"]] <- sapply(strsplit(z.119[["AC"]],split="-"),function(x) x[1])
cp22 <- rbind(cp22,drugModelALL("Dasatinib",z.119,"Z119","CID:3062316,CHEMBL1421,CHEBI:49375,DB01254"))
cp22 <- rbind(cp22,drugModelALL("Bosutinib",z.119,"Z119","CID:5328940,CHEMBL288441,CHEBI:39112"))
cp22 <- rbind(cp22,drugModelALL("Bafetinib",z.119,"Z119","CID:24853523"))
cp22 <- rbind(cp22,drugModelALL("Nilotinib",z.119,"Z119","CID:644241,CHEBI:52172"))

write.table(cp22,file="dpi/dpi-20130846-plosone-rix-colinge.txt",quote=F,sep="\t",row.names=F)

# Old timers ===============================================================================

# First TKI paper by Uwe, 17720881, 2007 Dec
drugModel <- function(targets,drug,CL,pmid,date,cid=NA,plot=F){
  targets[targets[["ID"]]=='bcrabl',"AC"] <- bcrabl
  t <- targets[["AC"]]
  weights <- targets[["sc"]]
  strong <- rep('no',length(weights))
  if (plot)
    boxplot(weights,main=drug)
  bs <- boxplot.stats(weights)
  strong[weights%in%bs$out] <- 'yes'
  strong[weights>=0.5*max(weights)] <- 'yes'
  len <- length(t)
  data.frame(source=rep(drug,len),target=t,source_name=rep(drug,len),target_name=sp.names[t],PMIDs=rep(pmid,len),
             dates=rep(date,len),sources=rep(pmid,len),type=rep('DPI',len),source_IDs=rep(cid,len),target_form=rep("",len),strong_target=strong,
             strength=weights,strength_type=rep('spectral.count',len),strength_direction=rep("high",len),
             technology=rep('chemprot',len),cell_type=rep(CL,len),stringsAsFactors=F) 
}

ima <- read.csv("/home/colinge/comp-chem-prot/smallCML/pulldowns/ima-edited-targets-kinases.txt",sep="\t",check.names=F,stringsAsFactors=F)
old.tki <- drugModel(ima,"Imatinib","K562",'17720881',2007,"CID:5291,CHEMBL941,CHEBI:45783,DB00619,DB03261",plot=T)
dasa <-  read.csv("/home/colinge/comp-chem-prot/smallCML/pulldowns/dasa-edited-targets-kinases.txt",sep="\t",check.names=F,stringsAsFactors=F)
old.tki <- rbind(old.tki,drugModel(dasa,"Dasatinib","K562",'17720881',2007,"CID:3062316,CHEMBL1421,CHEBI:49375,DB01254",plot=T))
nilo <- read.csv("nilotinib.txt",sep="\t",check.names=F,stringsAsFactors=F)
targets0 <- nilo[[1]]
t <- names.sp[targets0]
names(t) <- NULL
t[targets0=='ABL1'] <- bcrabl
weights <- nilo[[3]]
boxplot(weights,main="Nilotinib")
bs <- boxplot.stats(weights)
strong <- rep('no',length(weights))
strong[weights%in%bs$out] <- 'yes'
strong[weights>=0.5*max(weights)] <- 'yes'
len <- length(t)
old.tki <- rbind(old.tki,
                 data.frame(source=rep("Nilotinib",len),target=t,source_name=rep("Nilotinib",len),target_name=sp.names[t],PMIDs=rep('17720881',len),dates=rep(2007,len),sources=rep('17720881',len),
                            type=rep('DPI',len),source_IDs=rep("CID:644241,CHEBI:52172",len),target_form=rep("",len),strong_target=strong,strength=weights,strength_type=rep('peptide.count',len),strength_direction=rep("high",len),
                            technology=rep('chemprot',len),cell_type=rep("K562",len),stringsAsFactors=F)
)

# Bafe, Uwe, 19890374, 2010 Jan
inno <-  read.csv("/home/colinge/comp-chem-prot/smallCML/pulldowns/inno-edited-targets-kinases.txt",sep="\t",check.names=F,stringsAsFactors=F)
old.tki <- rbind(old.tki,drugModel(inno,"Bafetinib","K562",'19890374',2010,"CID:24853523",plot=T))

# Bosu, Lily, 19039322, 2009 Mar
bosu <- read.csv("/home/colinge/comp-chem-prot/smallCML/pulldowns/bosu-edited-targets-kinases.txt",sep="\t",check.names=F,stringsAsFactors=F)
old.tki <- rbind(old.tki,drugModel(bosu,"Bosutinib","K562",'19039322',2009,"CID:5328940,CHEMBL288441,CHEBI:39112",plot=T))

write.table(old.tki,file="dpi/dpi-old-cemm-papers.txt",quote=F,sep="\t",row.names=F)

# =============================================

# Georg's Nature Chem Biol (mouse data)

kin8 <- read.csv("nchembio.1085-S4.txt",sep="\t",check.names=F,stringsAsFactors=F)
boxplot(10^kin8[,13:20])
plot(density(10^c(kin8[[13]],kin8[[14]],kin8[[17]],kin8[[18]])))
thres <- quantile(10^c(kin8[[13]],kin8[[14]],kin8[[17]],kin8[[18]]),prob=0.95)
kin8.targets <- unlist(lapply(strsplit(kin8[["Ac"]],split="[,-]"),function(x){x[1]}))
kin8.names <- sp.names[kin8.targets]
# Stop because of mouse

# strong targets
node.type <- read.csv("figures/node-type.txt",sep="\t",stringsAsFactors=F)
node.type[[2]][node.type[[2]]=='tar'] <- 'kin'
strong <- read.csv("dpi/strong-unique.txt",sep="\t",header=F,stringsAsFactors=F)[[1]]
node.type[[2]][node.type[[1]]%in%strong] <- 'tar'
write.table(node.type,quote=F,sep="\t",row.names=F,file="figures/node-type.txt")
