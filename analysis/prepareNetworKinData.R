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
sp.names["Q9Y3C1"] <- "NOP16"
id.sp <- sp2name.3[["Accession"]]
names(id.sp) <- sp2name.3[["Entry name"]]

# Reads NetworKin KSIs =======================

nki0 <- read.csv("formated_networKIN_result.tsv",check.names=F,sep="\t",stringsAsFactors=F)[,c("Kinase","Accession",
                                                                                               "Kinase_Family","Motif_Score",
                                                                                               "Context_Score","Ranking_Score",
                                                                                               "Position_Numeric")]
nki.kin0 <- strsplit(nki0[["Kinase"]],split="\\|")
nki.kin.name <- unlist(lapply(nki.kin0,function(x){x[1]}))
nki.kin <- names.sp[nki.kin.name]
nki0[,"Kinase"] <- nki.kin
good <- !is.na(nki.kin)
nki <- cbind(nki0[good,],data.frame(Kinase.name=nki.kin.name[good],stringsAsFactors=F))
names(nki)[2] <- "Substrate"

# Reads reference KSIs =======================

interactdb <- read.csv("ksi-interactdb-20140211.txt",check.names=F,sep="\t",stringsAsFactors=F)[,c("source","target")]
ksr0 <- read.csv("comKSR.txt",sep="\t",stringsAsFactors=F)
ksr.kin <- names.sp[ksr0[["Kinase"]]]
ksr.sub <- names.sp[ksr0[["Substrate"]]]
ksr1 <- ksr0[!is.na(ksr.kin) & !is.na(ksr.sub),]
ksr <- data.frame(Kinase=names.sp[ksr1[["Kinase"]]],Substrate=names.sp[ksr1[["Substrate"]]],stringsAsFactors=F)
elm <- read.csv("ksi-phosphoELM_vertebrate_2011-11.dump",check.names=F,sep="\t",stringsAsFactors=F)[,1:2]
pplus <- read.csv("ksi-phosphositeplus.txt",check.names=F,sep="\t",stringsAsFactors=F)[,1:2]

# Writes Newman et al. KSIs
ksr.w <- ksr
names(ksr.w) <- c("source","target")
len <- dim(ksr.w)[1]
write.table(cbind(ksr.w,data.frame(source_name=sp.names[ksr.w[["source"]]],target_name=sp.names[ksr.w[["target"]]],
                                   PMIDs=rep("23549483",n=len),dates=rep(2013,n=len),sources=rep("23549483",n=len),
                                   type=rep("KSI",n=len),phospho_positions=rep("NA",n=len))),
            file="ksi-newman-et-al.txt",quote=F,sep="\t",row.names=F)

# Analysis of the overlap
safe.keys <- unique(c(paste(interactdb[[1]],interactdb[[2]],sep="||"),paste(ksr[[1]],ksr[[2]],sep="||"),
                      paste(pplus[[1]],pplus[[2]],sep="||"),paste(elm[[1]],elm[[2]],sep="||")))
nki.keys <- paste(nki[[1]],nki[[2]],sep="||")
length(intersect(nki.keys,safe.keys))
in.safe <- nki.keys %in% safe.keys

pairs(nki[,c("Context_Score","Ranking_Score","Motif_Score")])
pca <- prcomp(nki[,c("Context_Score","Ranking_Score","Motif_Score")],scale=T)
pca$rotation
pdf("figures/fig-1-pca-percvar.pdf",width=2,height=2,pointsize=8)
plot(pca)
dev.off()
summary(pca)
pdf("figures/fig-1-pca-scatter.pdf",width=4,height=4,pointsize=8)
png("figures/fig-1-pca-scatter.png",width=4,height=4,pointsize=8,res=600,units="in")
symbols(x=pca$x[,1],y=pca$x[,2],circles=rep(0.03,dim(pca$x)[1]),inches=F,bg="grey",fg=NULL)
symbols(x=pca$x[in.safe,1],y=pca$x[in.safe,2],circles=rep(0.04,sum(in.safe)),inches=F,bg="red",fg=NULL,add=T)
x0 <- 3
y0 <- -5
fact <- 1.7
x1 <- x0+fact*pca$rotation[1,1]
y1 <- y0+fact*pca$rotation[1,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"context",pos=3)
x1 <- x0+fact*pca$rotation[2,1]
y1 <- y0+fact*pca$rotation[2,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"ranking",pos=2)
x1 <- x0+fact*pca$rotation[3,1]
y1 <- y0+fact*pca$rotation[3,2]
arrows(x0=x0,y0=y0,x1=x1,y1=y1,length=0.15,col="black")
text(x1,y1,"motif",pos=2)
dev.off()

score <- nki[["Context_Score"]]
d.safe <- density(score[in.safe],n=128,from=min(score),to=max(score))
d.all <- density(score,n=128,from=min(score),to=max(score))
plot(density(score[!in.safe]))
lines(density(score[in.safe]),col="red")
lines(x=d.safe$x,y=d.all$y-d.safe$y,col="blue")

score <- nki[["Context_Score"]]
min(score[in.safe])
thres.context <- quantile(score[in.safe],prob=0.15)
thres.context
pdf("figures/fig-1-context-score.pdf",pointsize=8,width=4,height=2)
plot(density(score[!in.safe]),ylim=c(0,25))
lines(density(score[in.safe]),col="red")
abline(v=thres.context,col="blue")
dev.off()

score <- nki[["Ranking_Score"]]
min(score[in.safe])
thres.rank <- quantile(score[in.safe],prob=0.15)
thres.rank
pdf("figures/fig-1-ranking-score.pdf",pointsize=8,width=4,height=2)
plot(density(score[!in.safe]))
lines(density(score[in.safe]),col="red")
abline(v=thres.rank,col="blue")
dev.off()

score <- nki[["Motif_Score"]]
min(score[in.safe])
thres.motif <- quantile(score[in.safe],prob=0.15)
thres.motif
pdf("figures/fig-1-motif-score.pdf",pointsize=8,width=4,height=2)
plot(density(score[!in.safe]))
lines(density(score[in.safe]),col="red")
abline(v=thres.motif,col="blue")
dev.off()

biplot(pca,choices=1:2)
plot(pca$x[,1:2],type="p",pch=".",col="grey")
lines(pca$x[in.safe,1:2],type="p",pch=20,col="red")

# Writing into the ksi- file
nki.trusted <- nki[(nki[["Context_Score"]]>thres.context) & (nki[["Motif_Score"]]>thres.motif) &
                     !is.na(sp.names[nki[["Kinase"]]]) & !is.na(sp.names[nki[["Substrate"]]]),] # non-converted ACs usually indicate delete UniProt entries
len <- dim(nki.trusted)[1]
date <- 2007
write.table(data.frame(source=nki.trusted[["Kinase"]],target=nki.trusted[["Substrate"]],source_name=sp.names[nki.trusted[["Kinase"]]],
                       target_name=sp.names[nki.trusted[["Substrate"]]],PMIDs=rep("17570479",n=len),
                       dates=rep(date,n=len),sources=rep("networkin",n=len),
                       type=rep("computeKSI"),phospho_positions=nki.trusted[["Position_Numeric"]]),
            file="ksi-networkin.txt",quote=F,sep="\t",row.names=F)

# LDA (not appropriate really)
library(MASS)
non.trusted <- sample((1:dim(nki)[1])[!in.safe],sum(in.safe))
trusted <- (1:dim(nki)[1])[in.safe]
indexes <- c(trusted,non.trusted)
mdat <- cbind(nki[indexes,c("Context_Score","Ranking_Score","Motif_Score")],data.frame(trusted=c(rep(TRUE,length(trusted)),rep(FALSE,length(non.trusted))),stringsAsFactors=F))
fit <- lda(trusted ~ Context_Score+Ranking_Score+Motif_Score,mdat,na.action="na.omit",CV=TRUE)
plot(fit)
ct <- table(mdat$trusted, fit$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

library(klaR)
pdf("partimat.pdf")
partimat(mdat[,1:3],grouping=as.factor(mdat[[4]]),method="qda",plot.matrix=T,gs=".") 
dev.off()



