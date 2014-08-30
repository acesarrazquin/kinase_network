setwd("/home/colinge/kinase-network")

library(SparseM)
library(RBGL)
library(Rgraphviz)
library(plotrix)

source("/home/colinge/net-r/expandGraph.R")

# Reference PPI network -------------------------------------
alldb <- read.csv("ppi-20140211.txt",sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
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

# Human kinases --------------------------------------------
kinases <- read.csv("UniProt_human_protein_kinases_formated-atypical-1class.txt",sep="\t",stringsAsFactors=FALSE)
rownames(kinases) <- kinases[["UniProt_ID"]]

class.size <- table(kinases[["Type"]])
kin.class <- kinases[["Type"]]
names(kin.class) <- kinases[["UniProt_ID"]]
classes <- unique(kinases[["Type"]])
n.cl <- length(classes)

fam.0 <- read.csv("uniprot_kin_kinomefam-20140221.txt",sep="\t",stringsAsFactors=FALSE)
fam <- fam.0[!is.na(fam.0[["fam"]]),]
kin.fam <- paste(fam[["uni_type"]],fam[["fam"]],sep=":")[]
names(kin.fam) <- fam[["uniprot_id"]][!is.na(fam[["fam"]])]
fam.size <- table(kin.fam)

# Reference KSI network ----------------------------------
ksi <- read.table("ksi-20140211.txt",sep="\t",stringsAsFactors=FALSE,header=T)
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

is.kinase <- ksi[["source"]]%in%kinases[["UniProt_ID"]] & ksi[["target"]]%in%kinases[["UniProt_ID"]]
from <- as.vector(ksi[["source"]][k.no.self & is.kinase])
to <- as.vector(ksi[["target"]][k.no.self & is.kinase])
ft <- cbind(from,to)
kin.ksi.g <- ftM2graphNEL(ft, edgemode="directed")
kin.ksi.cc <- connectedComp(kin.ksi.g)
kin.ksi.G <- subGraph(kin.ksi.cc[[1]],kin.ksi.g)
n.kin.ksi.g <- nodes(kin.ksi.g)

# ID mapping ===============================
bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
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

# Generation of subsets ========================================================

ppi.deg <- degree(alldb.g)
ksi.deg <- degree(ksi.g)
kin.ksi.deg <- degree(kin.ksi.g)

kin5 <- intersect(names(ksi.deg$outDegree[ksi.deg$outDegree>=5]),kinases[["UniProt_ID"]])
kin10 <- intersect(names(ksi.deg$outDegree[ksi.deg$outDegree>=10]),kinases[["UniProt_ID"]])
kin.set <- kin10

di <- ksi.deg$inDegree
do <- ksi.deg$outDegree
kdi <- kin.ksi.deg$inDegree
kdo <- kin.ksi.deg$outDegree
kinlist <- union(kinases[["UniProt_ID"]],ksi[["source"]])

strong <- unique(read.csv("dpi/allstrongdpi-curated-20140305.txt",sep="\t",stringsAsFactors=FALSE)[[2]])

min.n <- 0

# controlled
nonkin.ctrled <- setdiff(names(di[di>=quantile(di,prob=0.95)]),kinlist)
kr <- kdi/(kdi+kdo)
kin.ctrled <- intersect(names(kr[kr>=0.8 & (kdi+kdo)>=min.n]),kin.set)
length(kin.ctrled)

# controllers
r <- do/(di+do)
nonkin.ctrler <- intersect(names(r[r>=0.8 & (di+do)>=min.n]),kin.set)
kr <- kdo/(kdi+kdo)
kin.ctrler <- intersect(names(kr[kr>=0.8 & (kdi+kdo)>=min.n]),kin.set)
length(kin.ctrler)
length(nonkin.ctrler)
length(intersect(kin.ctrler,nonkin.ctrler))

# integrators
nonkin.int <- intersect(names(r[r<0.66 & r>0.33 & (di+do)>=min.n]),kin.set)
kin.int <- intersect(names(kr[kr<0.66 & kr>0.33 & (kdi+kdo)>=min.n]),kin.set)
length(nonkin.int)
length(kin.int)
length(intersect(kin.int,nonkin.int))
length(intersect(kin.int,nonkin.ctrler))
length(intersect(nonkin.int,nonkin.ctrler))

length(intersect(nonkin.int,kin.ctrled))
length(intersect(nonkin.ctrler,kin.ctrled))

length(intersect(strong,nonkin.ctrler))
length(intersect(strong,kin.ctrler))
length(intersect(strong,nonkin.int))
length(intersect(strong,kin.int))
length(intersect(strong,kin.ctrled))

write.table(kin.ctrler,file="prot-lists-2/kinase-controllers.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.ctrler,file="prot-lists-2/controllers.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(kin.ctrled,file="prot-lists-2/controlled-kinases.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.ctrled,file="prot-lists-2/controlled.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(kin.int,file="prot-lists-2/kinase-integrators.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.int,file="prot-lists-2/integrators.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

symbols(x=di[kin.set],y=do[kin.set],circles=rep(0.2,length(kin.set)),inches=FALSE,bg="black",fg=NULL)
symbols(x=di[kin.ctrled],y=do[kin.ctrled],circles=rep(0.2,length(kin.ctrled)),inches=FALSE,bg="blue",fg=NULL,add=T)
symbols(x=di[kin.ctrler],y=do[kin.ctrler],circles=rep(0.2,length(kin.ctrler)),inches=FALSE,bg="red",fg=NULL,add=T)
symbols(x=di[kin.int],y=do[kin.int],circles=rep(0.2,length(kin.int)),inches=FALSE,bg="green",fg=NULL,add=T)

sp.names[nonkin.ctrler]
A <- sp.names[setdiff(nonkin.ctrler,c(kin.ctrler,nonkin.int,kin.int,kin.ctrled))]
B <- sp.names[setdiff(kin.ctrler,c(nonkin.int,kin.int,kin.ctrled))]
C <- sp.names[setdiff(intersect(nonkin.ctrler,kin.int),nonkin.int)]
D <- sp.names[intersect(intersect(nonkin.ctrler,kin.int),nonkin.int)]
E <- sp.names[intersect(nonkin.int,kin.int)]
Ff <- sp.names[setdiff(intersect(nonkin.ctrler,nonkin.int),c(kin.int,kin.ctrled))]
G <- sp.names[setdiff(intersect(kin.ctrled,nonkin.int),nonkin.ctrler)]
H <- sp.names[intersect(intersect(kin.ctrled,nonkin.int),nonkin.ctrler)]
I <- sp.names[setdiff(intersect(kin.ctrled,nonkin.ctrler),nonkin.int)]
J <- sp.names[setdiff(kin.ctrled,c(nonkin.ctrler,nonkin.int))]
K <- sp.names[setdiff(nonkin.int,c(nonkin.ctrler,kin.ctrled,kin.int))]
L <- sp.names[setdiff(kin.int,c(nonkin.ctrler,nonkin.int))]


# PCA -----
kin.profile <- data.frame(ksi.in=di[kin.set],ksi.out=do[kin.set],ksi.N=(di+do)[kin.set],ppi=ppi.deg[kin.set],ratio=log10((do/(di+do))[kin.set]*2))
#kin.profile <- data.frame(ksi.in=log10(1+di[kin.set]),ksi.out=log10(1+do[kin.set]),ksi.N=log10(1+(di+do)[kin.set]),ppi=log10(1+ppi.deg[kin.set]),ratio=log10((do/(di+do))[kin.set]*2))
kin.profile[["ppi"]][is.na(kin.profile[["ppi"]])] <- median(kin.profile[["ppi"]],na.rm=T)
pairs(kin.profile[,1:4])
pca <- prcomp(~ ksi.in+ksi.out+ratio,kin.profile,scale=T)
pca$rotation
plot(pca)
summary(pca)

symbols(x=pca$x[,1],y=pca$x[,2],circles=rep(0.05,dim(pca$x)[1]),inches=FALSE,bg="black",fg=NULL)
symbols(x=pca$x[kin.ctrled,1],y=pca$x[kin.ctrled,2],circles=rep(0.05,length(kin.ctrled)),inches=FALSE,bg="blue",fg=NULL,add=T)
symbols(x=pca$x[nonkin.ctrler,1],y=pca$x[nonkin.ctrler,2],circles=rep(0.05,length(nonkin.ctrler)),inches=FALSE,bg="red",fg=NULL,add=T)
symbols(x=pca$x[kin.int,1],y=pca$x[kin.int,2],circles=rep(0.05,length(kin.int)),inches=FALSE,bg="green",fg=NULL,add=T)
i <- intersect(nonkin.ctrler,kin.int)
if (length(i)>0)
  symbols(x=pca$x[i,1],y=pca$x[i,2],circles=rep(0.05,length(i)),inches=FALSE,bg="yellow",fg="green",add=T)
i <- intersect(kin.ctrled,kin.int)
if (length(i)>0)
  symbols(x=pca$x[kin.int,1],y=pca$x[kin.int,2],circles=rep(0.05,length(kin.int)),inches=FALSE,bg="cyan",fg=NULL,add=T)

for (i in 1:dim(pca$x)[1]){
  if (pca$x[i,1]^2 + pca$x[i,2]^2 >= 3)
    text(pca$x[i,1],pca$x[i,2],sp.names[rownames(pca$x)[i]])
}

x0 <- 2.5
y0 <- -0.5
fact <- 1.7
x1 <- x0+fact*pca$rotation[1,1]
y1 <- y0+fact*pca$rotation[1,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"KSI in",pos=3)
x1 <- x0+fact*pca$rotation[2,1]
y1 <- y0+fact*pca$rotation[2,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"KSI out",pos=3)
x1 <- x0+fact*pca$rotation[3,1]
y1 <- y0+fact*pca$rotation[3,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"ratio",pos=3)
# x1 <- x0+fact*pca$rotation[4,1]
# y1 <- y0+fact*pca$rotation[4,2]
# arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
# text(x1,y1,"PPI",pos=3)

# GOBP analysis ------------------------------

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

sp2g <- nr.sprot2GO.bp#[intersect(n.kin.part.G,names(nr.sprot2GO.bp))]
go.bps <- setdiff(unique(unlist(sp2g)),c("all","GO:0008150"))

goOneSet <- function(prots){
  go <- NULL
  for (bp in go.bps){
    match <- unlist(lapply(sp2g,function(x){sum(x %in% bp)}))
    inpw <- names(match[match>0])
    pval <- phyper(q=length(intersect(prots,inpw))-1,m=length(inpw),n=20488-length(inpw),k=length(prots),lower.tail=FALSE)
    go <- rbind(go,data.frame(term=bp,p.val=pval,num=length(inpw),num.prots=length(prots),stringsAsFactors=FALSE))
  }
  go
}

cc.go <- row.labels.gen[go.bps]
cc.go <- cbind(cc.go,goOneSet(kin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goOneSet(kin.int)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.int)[,2:4])
cc.go <- cbind(cc.go,goOneSet(kin.ctrled)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.ctrled)[,2:4])

write.table(cc.go,file="go-bp-hm.txt",sep="\t",quote=FALSE,row.names=FALSE)

utile <- cc.go[,c(1,2+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[minp<1e-3,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- matr[[1]]
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")

pdf("figures/fig6-heatmap-BP.pdf",width=10,height=6,pointsize=8)
heatmap(-log(num.matr),scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()


# =========================================================
# Experimental KSIs only
# =========================================================

# Reference KSI network ----------------------------------
eksi <- ksi[ksi[["type"]]=="KSI",]
k.no.self <- as.vector(eksi[["source"]]) != as.vector(eksi[["target"]])
from <- as.vector(eksi[["source"]][k.no.self])
to <- as.vector(eksi[["target"]][k.no.self])
ft <- cbind(from,to)
eksi.g <- ftM2graphNEL(ft, edgemode="directed")
n.eksi.g <- nodes(eksi.g)
eksi.cc <- connectedComp(eksi.g)
eksi.n.nodes <- length(eksi.cc[[1]])
eksi.G <- subGraph(eksi.cc[[1]],eksi.g)
n.eksi.G <- nodes(eksi.G)

is.kinase <- eksi[["source"]]%in%kinases[["UniProt_ID"]] & eksi[["target"]]%in%kinases[["UniProt_ID"]]
from <- as.vector(eksi[["source"]][k.no.self & is.kinase])
to <- as.vector(eksi[["target"]][k.no.self & is.kinase])
ft <- cbind(from,to)
kin.eksi.g <- ftM2graphNEL(ft, edgemode="directed")
kin.eksi.cc <- connectedComp(kin.eksi.g)
kin.eksi.G <- subGraph(kin.eksi.cc[[1]],kin.eksi.g)
n.kin.eksi.g <- nodes(kin.eksi.g)

eksi.deg <- degree(eksi.g)
kin.eksi.deg <- degree(kin.eksi.g)

kin5 <- intersect(names(eksi.deg$outDegree[eksi.deg$outDegree>=5]),kinases[["UniProt_ID"]])
kin10 <- intersect(names(eksi.deg$outDegree[eksi.deg$outDegree>=10]),kinases[["UniProt_ID"]])
kin.set <- kin10

di <- eksi.deg$inDegree
do <- eksi.deg$outDegree
kdi <- kin.eksi.deg$inDegree
kdo <- kin.eksi.deg$outDegree
kinlist <- union(kinases[["UniProt_ID"]],eksi[["source"]])


min.n <- 0

# controlled
nonkin.ctrled <- setdiff(names(di[di>=quantile(di,prob=0.95)]),kinlist)
kr <- kdi/(kdi+kdo)
kin.ctrled <- intersect(names(kr[kr>=0.8 & (kdi+kdo)>=min.n]),kin.set)
length(kin.ctrled)

# controllers
r <- do/(di+do)
nonkin.ctrler <- intersect(names(r[r>=0.8 & (di+do)>=min.n]),kin.set)
kr <- kdo/(kdi+kdo)
kin.ctrler <- intersect(names(kr[kr>=0.8 & (kdi+kdo)>=min.n]),kin.set)
length(kin.ctrler)
length(nonkin.ctrler)
length(intersect(kin.ctrler,nonkin.ctrler))

# integrators
nonkin.int <- intersect(names(r[r<0.66 & r>0.33 & (di+do)>=min.n]),kin.set)
kin.int <- intersect(names(kr[kr<0.66 & kr>0.33 & (kdi+kdo)>=min.n]),kin.set)
length(nonkin.int)
length(kin.int)
length(intersect(kin.int,nonkin.int))
length(intersect(kin.int,nonkin.ctrler))
length(intersect(nonkin.int,nonkin.ctrler))

length(intersect(nonkin.int,kin.ctrled))
length(intersect(nonkin.ctrler,kin.ctrled))

length(intersect(strong,nonkin.ctrler))
length(intersect(strong,kin.ctrler))
length(intersect(strong,nonkin.int))
length(intersect(strong,kin.int))
length(intersect(strong,kin.ctrled))

write.table(kin.ctrler,file="prot-lists-3/kinase-controllers.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.ctrler,file="prot-lists-3/controllers.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(kin.ctrled,file="prot-lists-3/controlled-kinases.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.ctrled,file="prot-lists-3/controlled.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(kin.int,file="prot-lists-3/kinase-integrators.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(nonkin.int,file="prot-lists-3/integrators.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

symbols(x=di[kin.set],y=do[kin.set],circles=rep(0.2,length(kin.set)),inches=FALSE,bg="black",fg=NULL)
symbols(x=di[kin.ctrled],y=do[kin.ctrled],circles=rep(0.2,length(kin.ctrled)),inches=FALSE,bg="blue",fg=NULL,add=T)
symbols(x=di[kin.ctrler],y=do[kin.ctrler],circles=rep(0.2,length(kin.ctrler)),inches=FALSE,bg="red",fg=NULL,add=T)
symbols(x=di[kin.int],y=do[kin.int],circles=rep(0.2,length(kin.int)),inches=FALSE,bg="green",fg=NULL,add=T)

getPie <- function(prots,fname,width=2.5,height=2.5,pointsize=8,with.labels=FALSE){
  c <- table(kinases[kinases[["UniProt_ID"]]%in%prots,"Type"])
  c <- c[order(names(c))]
  cols <- c(rgb(0,1,1,1),rgb(1,0.6,0,1),rgb(0,0,1,1),rgb(1,0.6,0.6,1),rgb(0,1,0,1),
            rgb(0.6,0.8,0,1),rgb(0.8,0.8,0.8,1),rgb(1,1,0,1),rgb(1,0,1,1),rgb(1,0,0,1),rgb(0.4,0.4,1,1))
  names(cols) <- c("AGC","Atypical","CAMK","CK1","CMGC","NEK","Other","STE","TKL","Tyr","RGC")
  if (with.labels)
    labels <- names(c)
  else
    labels <- NA
  pdf(file=fname,width=width,height=height,useDingbats=FALSE,pointsize=pointsize)
  pie(c,labels=labels,col=cols[names(c)])
  dev.off()
}

getPie(nonkin.ctrler,"figures/fig-6-pie-nonkin-ctrler.pdf",with.labels=TRUE)
getPie(kin.ctrler,"figures/fig-6-pie-kin-ctrler.pdf")
getPie(kin.int,"figures/fig-6-pie-kin-int.pdf")
getPie(kin.ctrled,"figures/fig-6-pie-kin-ctrled.pdf")



sp.names[nonkin.ctrler]
A <- sp.names[setdiff(nonkin.ctrler,c(kin.ctrler,nonkin.int,kin.int,kin.ctrled))]
B <- sp.names[setdiff(kin.ctrler,c(nonkin.int,kin.int,kin.ctrled))]
C <- sp.names[setdiff(intersect(nonkin.ctrler,kin.int),nonkin.int)]
D <- sp.names[intersect(intersect(nonkin.ctrler,kin.int),nonkin.int)]
E <- sp.names[intersect(nonkin.int,kin.int)]
Ff <- sp.names[setdiff(intersect(nonkin.ctrler,nonkin.int),c(kin.int,kin.ctrled))]
G <- sp.names[setdiff(intersect(kin.ctrled,nonkin.int),nonkin.ctrler)]
H <- sp.names[intersect(intersect(kin.ctrled,nonkin.int),nonkin.ctrler)]
I <- sp.names[setdiff(intersect(kin.ctrled,nonkin.ctrler),nonkin.int)]
J <- sp.names[setdiff(kin.ctrled,c(nonkin.ctrler,nonkin.int))]
K <- sp.names[setdiff(nonkin.int,c(nonkin.ctrler,kin.ctrled,kin.int))]
L <- sp.names[setdiff(kin.int,c(nonkin.ctrler,nonkin.int))]

cat(intersect(A[order(A)],sp.names[strong]),sep=", ")
cat(setdiff(A[order(A)],sp.names[strong]),sep=", ")
cat(intersect(B[order(B)],sp.names[strong]),sep=", ")
cat(setdiff(B[order(B)],sp.names[strong]),sep=", ")
cat(intersect(C[order(C)],sp.names[strong]),sep=", ")
cat(setdiff(C[order(C)],sp.names[strong]),sep=", ")
cat(intersect(L[order(L)],sp.names[strong]),sep=", ")
cat(setdiff(L[order(L)],sp.names[strong]),sep=", ")
E
Ff
D
cat(intersect(K[order(K)],sp.names[strong]),sep=", ")
cat(setdiff(K[order(K)],sp.names[strong]),sep=", ")
cat(intersect(G[order(G)],sp.names[strong]),sep=", ")
cat(setdiff(G[order(G)],sp.names[strong]),sep=", ")
G
H
cat(intersect(I[order(I)],sp.names[strong]),sep=", ")
cat(setdiff(I[order(I)],sp.names[strong]),sep=", ")
cat(intersect(J[order(J)],sp.names[strong]),sep=", ")
cat(setdiff(J[order(J)],sp.names[strong]),sep=", ")
cat(sp.names[nonkin.ctrled[order(sp.names[nonkin.ctrled])]],sep=", ")


cc.go <- row.labels.gen[go.bps]
cc.go <- cbind(cc.go,goOneSet(kin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goOneSet(kin.int)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.int)[,2:4])
cc.go <- cbind(cc.go,goOneSet(kin.ctrled)[,2:4])
cc.go <- cbind(cc.go,goOneSet(nonkin.ctrled)[,2:4])

write.table(cc.go,file="go-expKSI-bpslim-hm.txt",sep="\t",quote=FALSE,row.names=FALSE)

utile <- cc.go[,c(1,2+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[minp<1e-3,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- matr[[1]]
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")

pdf("figures/fig6-heatmap-expKSI-BPslim.pdf",width=10,height=6,pointsize=8)
heatmap(-log(num.matr),scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()

# Full GO -------

sp2go.p <- read.csv("../databases/sp2go-bp-human.txt",check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
u <- unique(sp2go.p[,4:5])
go.descr.p <- u[[2]]
names(go.descr.p) <- u[[1]]
sp2go.f <- read.csv("../databases/sp2go-mf-human.txt",check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
sp2go.f <- sp2go.f[,c(1,2,3,5,4)] # last 2 columns in wrong order
u <- unique(sp2go.f[,4:5])
go.descr.f <- u[[2]]
names(go.descr.f) <- u[[1]]
sp2go.c <- read.csv("../databases/sp2go-cc-human.txt",check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
sp2go.c <- sp2go.c[,c(1,2,3,5,4)] # last 2 columns in wrong order
u <- unique(sp2go.c[,4:5]) # wrong order again
go.descr.c <- u[[2]]
names(go.descr.c) <- u[[1]]

sp2go <- sp2go.c[,c(1,4)]
go.descr <- go.descr.c


goTwoSet <- function(prots){
  go <- NULL
  for (term in names(go.descr)){
    inpw <- sp2go[sp2go[[2]]==term,1]
    pval <- phyper(q=length(intersect(prots,inpw))-1,m=length(inpw),n=20488-length(inpw),k=length(prots),lower.tail=FALSE)
    go <- rbind(go,data.frame(term=term,p.val=pval,num=length(inpw),num.inter=length(intersect(prots,inpw)),stringsAsFactors=FALSE))
  }
  go
}

cc.go <- names(go.descr)
cc.go <- cbind(cc.go,goTwoSet(kin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goTwoSet(nonkin.ctrler)[,2:4])
cc.go <- cbind(cc.go,goTwoSet(kin.int)[,2:4])
cc.go <- cbind(cc.go,goTwoSet(nonkin.int)[,2:4])
cc.go <- cbind(cc.go,goTwoSet(kin.ctrled)[,2:4])
cc.go <- cbind(cc.go,goTwoSet(nonkin.ctrled)[,2:4])

write.table(cc.go,file="go-expKSI-cc-hm.txt",sep="\t",quote=FALSE,row.names=FALSE)

utile <- cc.go[,c(1,2+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[minp<1e-4,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- paste(matr[[1]],go.descr[as.character(matr[[1]])])
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")
top <- quantile(-log10(num.matr),prob=0.98)
num.matr[-log10(num.matr)>top] <- 10^(-top)
pdf("figures/fig6-heatmap-expKSI-CC.pdf",width=10,height=6,pointsize=8)
heatmap(-log(num.matr),scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()

utile <- cc.go[,c(1,4+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[maxp>=5,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- paste(matr[[1]],go.descr[as.character(matr[[1]])])
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")
pdf("figures/fig6-heatmap-expKSI-CC-num.pdf",width=10,height=6,pointsize=8)
heatmap(num.matr,scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()



sp2go <- sp2go.p[,c(1,4)]
go.descr <- go.descr.p
bp.go <- names(go.descr)
bp.go <- cbind(bp.go,goTwoSet(kin.ctrler)[,2:4])
bp.go <- cbind(bp.go,goTwoSet(nonkin.ctrler)[,2:4])
bp.go <- cbind(bp.go,goTwoSet(kin.int)[,2:4])
bp.go <- cbind(bp.go,goTwoSet(nonkin.int)[,2:4])
bp.go <- cbind(bp.go,goTwoSet(kin.ctrled)[,2:4])
bp.go <- cbind(bp.go,goTwoSet(nonkin.ctrled)[,2:4])
write.table(bp.go,file="go-expKSI-bp-hm.txt",sep="\t",quote=FALSE,row.names=FALSE)

utile <- bp.go[,c(1,2+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[minp<1e-5,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- paste(matr[[1]],go.descr[as.character(matr[[1]])])
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")
top <- quantile(-log10(num.matr),prob=0.95)
num.matr[-log10(num.matr)>top] <- 10^(-top)

pdf("figures/fig6-heatmap-expKSI-BP.pdf",width=10,height=6,pointsize=8)
heatmap(-log(num.matr),scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()

sp2go <- sp2go.f[,c(1,4)]
go.descr <- go.descr.f
mf.go <- names(go.descr)
mf.go <- cbind(mf.go,goTwoSet(kin.ctrler)[,2:4])
mf.go <- cbind(mf.go,goTwoSet(nonkin.ctrler)[,2:4])
mf.go <- cbind(mf.go,goTwoSet(kin.int)[,2:4])
mf.go <- cbind(mf.go,goTwoSet(nonkin.int)[,2:4])
mf.go <- cbind(mf.go,goTwoSet(kin.ctrled)[,2:4])
mf.go <- cbind(mf.go,goTwoSet(nonkin.ctrled)[,2:4])
write.table(mf.go,file="go-expKSI-mf-hm.txt",sep="\t",quote=FALSE,row.names=FALSE)

utile <- mf.go[,c(1,2+3*(0:5))]
minp <- apply(utile[,-1],1,min)
maxp <- apply(utile[,-1],1,max)
matr <- utile[minp<1e-4,]
num.matr <- as.matrix(matr[,-1])
rownames(num.matr) <- paste(matr[[1]],go.descr[as.character(matr[[1]])])
colnames(num.matr) <- c("kin ctrler","ctrler","kin integr","integr","kin ctrled","ctrled")
top <- quantile(-log10(num.matr),prob=0.95)
num.matr[-log10(num.matr)>top] <- 10^(-top)

pdf("figures/fig6-heatmap-expKSI-MF.pdf",width=10,height=6,pointsize=8)
heatmap(-log(num.matr),scale="none",Colv=NA,col=rev(heat.colors(20)))
dev.off()

# Oncogenes/tumor suppressors
onco <- read.csv("../essentialome/gene-families/oncogenes.txt",header=F,stringsAsFactors=F)[[1]]
suppr <- read.csv("../essentialome/gene-families/tumor-suppressors.txt",header=F,stringsAsFactors=F)[[1]]
growth <- read.csv("../essentialome/gene-families/cytokines-and-growth-factors.txt",header=F,stringsAsFactors=F)[[1]]
census <- read.csv("diseases//cancer_gene_census.txt",stringsAsFactors=F,check.names=F,sep="\t")[[1]]
sp.kegg <- read.csv("/home/colinge/databases/uniprot_to_kegg_human-2012-08-23.csv",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c('character','character'))
kegg2names <- read.csv("/home/colinge/databases/kegg_names-2012-08-23.txt",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c("character","character"))
kegg.names <- kegg2names[[2]]
names(kegg.names) <- kegg2names[[1]]
k.cancer <- NULL
for (term in unique(sp.kegg[[2]]))
  if (length(intersect(c('cancer','leukemia','carcinoma','Glioma','Melanoma'),unlist(strsplit(kegg.names[term],split=" ")))) >= 1){
    cat(term,kegg.names[term],"\n")
    mem <- sp.kegg[sp.kegg[[2]]==term,1]
    #inpw <- intersect(kin.names[["AC"]],mem)    
    k.cancer <- c(k.cancer,mem)
  }
k.cancer <- unique(k.cancer)
intersect(onco,sp.names[nonkin.ctrled])
intersect(suppr,sp.names[nonkin.ctrled])
intersect(growth,sp.names[nonkin.ctrled])
intersect(census,sp.names[nonkin.ctrled])
intersect(sp.names[k.cancer],sp.names[nonkin.ctrled])
all.c <- unique(c(census,onco,suppr,sp.names[k.cancer]))
all.c <- unique(c(census,onco,suppr))
intersect(all.c,sp.names[nonkin.ctrled])

phyper(q=length(intersect(all.c,sp.names[nonkin.ctrled]))-1,m=length(all.c),n=20488-length(all.c),k=length(nonkin.ctrled),lower.tail=FALSE)


# Counts of drug targets ==============================

targ <- read.csv("dpi/nb_strong_dpi_kinase.txt",sep="\t",stringsAsFactors=F)
n.strong <- 1261
n.targets <- dim(targ)[1]
n.kinases <- sim(kinases)[1]

cont.tab <- matrix(c(length(nonkin.ctrler),dim(kinases)[1],sum(targ[targ[["id"]]%in%nonkin.ctrler,"dpi"]),
                   sum(targ[!targ[["id"]]%in%nonkin.ctrler,"dpi"])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])

cont.tab <- matrix(c(length(kin.ctrler),dim(kinases)[1],sum(targ[targ[["id"]]%in%kin.ctrler,"dpi"]),
                     sum(targ[!targ[["id"]]%in%kin.ctrler,"dpi"])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])

Ap <- setdiff(nonkin.ctrler,kin.ctrler)
cont.tab <- matrix(c(length(Ap),dim(kinases)[1],sum(targ[targ[["id"]]%in%Ap,"dpi"]),
                     sum(targ[!targ[["id"]]%in%Ap,"dpi"])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])

cont.tab <- matrix(c(length(kin.ctrled),dim(kinases)[1],sum(targ[targ[["id"]]%in%kin.ctrled,"dpi"]),
                     sum(targ[!targ[["id"]]%in%kin.ctrled,"dpi"])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])

cont.tab <- matrix(c(length(kin.int),dim(kinases)[1],sum(targ[targ[["id"]]%in%kin.int,"dpi"]),
                     sum(targ[!targ[["id"]]%in%kin.int,"dpi"])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])

nn <- unlist(lapply(strsplit(names(fam.size),":"),function(x){x[2]}))
fam.size2 <- fam.size
names(fam.size2) <- nn
for (f in unique(targ[["fam"]])){
  cont.tab <- matrix(c(fam.size2[f],dim(kinases)[1],sum(targ[targ[["fam"]]==f,"dpi"],na.rm=T),
                       sum(targ[targ[["fam"]]!=f,"dpi"],na.rm=T)),ncol=2)
  if (sum(is.na(cont.tab))==0){
    res <- chisq.test(cont.tab)
    if (res$p.value < 0.05/129){
      cat(f,res$p.value*129,"\n")
    }
  }
}



