setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)
library(multtest)

source("/home/colinge/net-r/expandGraph.R")

# Full PPI & experimental KSI network
full <- read.csv("full-network-20140211.txt",sep="\t",check.names=FALSE,stringsAsFactors=F)
l.no.self <- as.vector(full[["source"]]) != as.vector(full[["target"]])
from.ppi <- as.vector(full[l.no.self & full[["type"]]=="PPI","source"])
to.ppi <- as.vector(full[l.no.self & full[["type"]]=="PPI","target"])
from.ksi <- as.vector(full[l.no.self & full[["type"]]=="KSI","source"])
to.ksi <- as.vector(full[l.no.self & full[["type"]]=="KSI","target"])

# Kinases ---------------------------
kinases <- read.csv("UniProt_human_protein_kinases_formated-atypical-1class.txt",sep="\t",stringsAsFactors=F)
kin.class <- kinases[["Type"]]
names(kin.class) <- kinases[["UniProt_ID"]]
classes <- unique(kinases[["Type"]])

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
kfrom.ksi <- as.vector(kin.part[kl.no.self & kin.part[["type"]]=="KSI","source"])
kto.ksi <- as.vector(kin.part[kl.no.self & kin.part[["type"]]=="KSI","target"])
sinks <- setdiff(kto.ksi,c(kfrom.ppi,kto.ppi,kfrom.ksi))
ft <- cbind(c(kfrom.ppi,kto.ppi,kfrom.ksi,sinks,"sink"),c(kto.ppi,kfrom.ppi,kto.ksi,rep("sink",length(sinks)+1)))
kin.part.g <- ftM2graphNEL(ft, edgemode="directed")
kin.part.cc <- connectedComp(kin.part.g)
kin.part.n.nodes <- length(kin.part.cc[[1]])
kin.part.G <- subGraph(kin.part.cc[[1]],kin.part.g)
n.kin.part.G <- nodes(kin.part.G)
kin.part.am <- graph2SparseM(kin.part.G,useweights=F)
kin.part.P <- adjmat2transprob(kin.part.am)
save(kin.part.P,file="kin-part-exp-transmat.data")
load("kin-part-exp-transmat.data")
deg <- degree(kin.part.G)

# ID mapping ===============================
bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
sp2name.1 <- sp2name[sp2name[["Accession"]] %in% n.kin.part.G,]
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


# Selected kinases from SRPK1 (spliceosome), HIPK1/-2 (angiogenesis, cancer), WEE1 (cancer, novel inhibitor), BMX (cancer, novel probe compound)
dpi <- read.csv("dpi/allstrongdpi-curated-20140305.txt",sep="\t",stringsAsFactors=FALSE)
dpi[["source_name"]] <- dpi[["source_IDs"]]
targets <- names.sp[c("SRPK1","HIPK1","HIPK2","WEE1","BMX")]
lx <- limitProb(P=kin.part.P,nodes=n.kin.part.G,seeds=targets)
p.sink <- lx["sink"]
lx["sink"] <- 0
lx <- lx/(1-p.sink)
names(lx) <- n.kin.part.G
sel <- getNeighborhood(g=kin.part.G,lx=lx,n.add=30,target.mapped=targets,disease.mapped=intersect(kinases[["UniProt_ID"]],n.kin.part.G),symbols=sp.names)
top <- nodes(sel$neighb)
kikin <- full[full[["source"]]%in%top & full[["target"]]%in%top & l.no.self & full[["type"]]!="computeKSI",c("source","target","source_name","target_name","type")]
kikin <- rbind(kikin,unique(dpi[dpi[["target"]]%in%top & dpi[["strong_target"]]=="yes",c("source","target","source_name","target_name","type")]))
write.table(kikin,quote=F,row.names=F,sep="\t",file="test-nets/prop5-deges.txt")





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

ksi <- read.table("ksi-20140211.txt",sep="\t",stringsAsFactors=F,header=T)

prot <- names.sp["APC"]
close <- adj(alldb.G,prot)[[1]]
close.kin <- intersect(close,kinases[["UniProt_ID"]])
modulators <- ksi[ksi[["target"]]%in%close & ksi[["type"]]=="KSI",]
fav.kin <- table(modulators[["source"]])
top <- names(fav.kin[fav.kin>=4])
drugs <- unique(dpi[dpi[[2]]%in%top,c(1,2,7)])
names(drugs) <- c("source","target","type")

mini <- modulators[,c("source","target","type")]
mini <- rbind(mini,data.frame(source=rep(prot,length(close)),target=close,type=rep("PPI",length(close))))
acs <- union(mini[["source"]],mini[["target"]])
mini <- rbind(mini,drugs)
acs <- union(acs,drugs[[2]])
nodes <- data.frame(ac=c(acs,drugs[[1]]),name=c(sp.names[acs],drugs[[1]]),
                    node.type=c(ifelse(acs%in%kinases[["UniProt_ID"]],"kinase","prot"),rep("drug",dim(drugs)[1])),stringsAsFactors=F)
write.table(mini,file="test-nets/apc-edges.txt",quote=F,row.names=F,sep="\t")
write.table(nodes,file="test-nets/apc-nodes.txt",quote=F,row.names=F,sep="\t")
write.table(fav.kin,file="test-nets/apc-kin.txt",quote=F,row.names=F,sep="\t")



