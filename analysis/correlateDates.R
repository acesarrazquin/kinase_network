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

# Loads kinase dates -----------------------------
#kin.dates0 <- read.csv("UniProt_kin_pubmed.txt",check.names=F,sep="\t",stringsAsFactors=F) # mapped searching pubmed
kin.dates0 <- read.csv("pubmed2kinases.txt",check.names=F,sep="\t",stringsAsFactors=F) # NCBI table
#kin.exclude <- read.csv("kinases-dates-to-exclude.txt",header=F,check.names=F,sep="\t",stringsAsFactors=F)[[1]]
kin.mapped <- names.sp[kin.dates0[["Kinase"]]]
#kin.dates <- kin.dates0[!(kin.dates0[["Kinase"]]%in%kin.exclude) & !is.na(kin.mapped),-dim(kin.dates0)[2]] # also removes 2014
kin.dates <- kin.dates0[!is.na(kin.mapped),-dim(kin.dates0)[2]] # also removes 2014
kin.names <- data.frame(Kinase=kin.dates[["Kinase"]],AC=names.sp[kin.dates[["Kinase"]]],stringsAsFactors=F)
rownames(kin.names) <- kin.names[["AC"]]
kin.1sty <- apply(kin.dates[,-1],1,function(x){sum(cumsum(as.numeric(x))==0)})
kin.tot <- apply(kin.dates[,-1],1,sum)

kin.pub <- as.matrix(kin.dates[,-1])
rownames(kin.pub) <- kin.names[["Kinase"]]
colnames(kin.pub) <- names(kin.dates[-1])
colors <- rev(heat.colors(50))
o1 <- order(kin.1sty)
pdf("figures/kin-pub-history.pdf")
heatmap(log(1+kin.pub[o1,]),scale="none",Rowv=NA,Colv=NA,col=colors)
dev.off()
ot <- rev(order(kin.tot))
heatmap(log(1+kin.pub[ot,]),scale="none",Rowv=NA,Colv=NA,col=colors)

# Organizes kinase interaction dates
non.self <- alldb[l.no.self,]
kin.int <- non.self[non.self[["source"]]%in%kin.names[["AC"]] | non.self[["target"]]%in%kin.names[["AC"]],]
kin.int.pub <- matrix(rep(0,dim(kin.pub)[1]*dim(kin.pub)[2]),ncol=dim(kin.pub)[2])
rownames(kin.int.pub) <- rownames(kin.pub)
colnames(kin.int.pub) <- as.numeric(colnames(kin.pub))
for (i in 1:dim(kin.int)[1]){
  dates <- kin.int[i,"dates"]
  if (!is.na(dates)){
    k <- intersect(kin.int[i,c("source","target")],kin.names[["AC"]])
    for (year in trunc(as.numeric(unlist(strsplit(dates,",")))))
      if (!is.na(year) && year >= 1950 && year <= 2013)
        kin.int.pub[kin.names[k,"Kinase"],year-1949] <- kin.int.pub[kin.names[k,"Kinase"],year-1949]+1
  }
}

kin.int.1sty <- apply(kin.int.pub,1,function(x){sum(cumsum(as.numeric(x))==0)})
kin.int.tot <- apply(kin.int.pub,1,sum)

#kin.bait <- non.self[(non.self[["source"]]%in%kin.names[["AC"]] & non.self[["source_is_bait"]]=="yes") | (non.self[["target"]]%in%kin.names[["AC"]] & non.self[["target_is_bait"]]=="yes"),]
kin.bait <- non.self[(non.self[["source"]]%in%kin.names[["AC"]] & non.self[["source_is_bait"]]=="yes" & non.self[["target_is_bait"]]=="no") | (non.self[["target"]]%in%kin.names[["AC"]] & non.self[["target_is_bait"]]=="yes" & non.self[["source_is_bait"]]=="no"),]
kin.bait.pub <- matrix(rep(0,dim(kin.pub)[1]*dim(kin.pub)[2]),ncol=dim(kin.pub)[2])
rownames(kin.bait.pub) <- rownames(kin.pub)
colnames(kin.bait.pub) <- as.numeric(colnames(kin.pub))
for (i in 1:dim(kin.bait)[1]){
  dates <- kin.bait[i,"dates"]
  if (!is.na(dates)){
    k <- intersect(kin.bait[i,c("source","target")],kin.names[["AC"]])
    for (year in trunc(as.numeric(unlist(strsplit(dates,",")))))
      if (!is.na(year) && year >= 1950 && year <= 2013)
        kin.bait.pub[kin.names[k,"Kinase"],year-1949] <- kin.bait.pub[kin.names[k,"Kinase"],year-1949]+1
  }
}

kin.bait.1sty <- apply(kin.bait.pub,1,function(x){sum(cumsum(as.numeric(x))==0)})
kin.bait.tot <- apply(kin.bait.pub,1,sum)

low.kinases <- intersect(kin.names[["AC"]],names(common)[common<=3])
kin.low <- non.self[(non.self[["source"]]%in%low.kinases) | (non.self[["target"]]%in%low.kinases),]
kin.low.pub <- matrix(rep(0,dim(kin.pub)[1]*dim(kin.pub)[2]),ncol=dim(kin.pub)[2])
rownames(kin.low.pub) <- rownames(kin.pub)
colnames(kin.low.pub) <- as.numeric(colnames(kin.pub))
for (i in 1:dim(kin.low)[1]){
  dates <- kin.low[i,"dates"]
  if (!is.na(dates)){
    k <- intersect(kin.low[i,c("source","target")],kin.names[["AC"]])
    for (year in trunc(as.numeric(unlist(strsplit(dates,",")))))
      if (!is.na(year) && year >= 1950 && year <= 2013)
        kin.low.pub[kin.names[k,"Kinase"],year-1949] <- kin.low.pub[kin.names[k,"Kinase"],year-1949]+1
  }
}

kin.low.1sty <- apply(kin.low.pub,1,function(x){sum(cumsum(as.numeric(x))==0)})
kin.low.tot <- apply(kin.low.pub,1,sum)


heatmap(log(1+kin.int.pub[o1,]),scale="none",Rowv=NA,Colv=NA,col=colors)
o1i <- order(kin.int.1sty)
heatmap(log(1+kin.int.pub[o1i,]),scale="none",Rowv=NA,Colv=NA,col=colors)
heatmap(log(1+kin.int.pub[ot,]),scale="none",Rowv=NA,Colv=NA,col=colors)
oti <- rev(order(kin.int.tot))
heatmap(log(1+kin.int.pub[oti,]),scale="none",Rowv=NA,Colv=NA,col=colors)

pdf("figures/fig-2-knowledge-correlation.pdf",width=3,height=3,pointsize=8,useDingbats=FALSE)
plot(x=jitter(1+kin.tot),y=jitter(1+kin.int.tot),pch=20,type="n",log="xy")
rect(10^(par("usr")[1]),10^(par("usr")[3]),10^(par("usr")[2]),10^(par("usr")[4]),col="grey85")
grid(col="white",lwd=2,lty=1)
symbols(x=jitter(1+kin.tot),y=jitter(1+kin.int.tot),circles=rep(0.03,length(kin.int.tot)),inches=F,bg="black",fg=NULL,xlim=c(0,3.5),add=T)
reg <- lm(int.publi~publi,data.frame(publi=log10(1+kin.tot),int.publi=log10(1+kin.int.tot)))
abline(a=reg$coefficients[1],b=reg$coefficients[2],col="red",lwd=2)
symbols(x=jitter(1+kin.tot[kin.names[["AC"]]%in%low.kinases]),y=jitter(1+kin.int.tot[kin.names[["AC"]]%in%low.kinases]),
        circles=rep(0.03,length(kin.int.tot[kin.names[["AC"]]%in%low.kinases])),inches=F,bg="yellow",fg=NULL,xlim=c(0,3.5),add=T)
reg2 <- lm(int.publi~publi,data.frame(publi=log10(1+kin.tot[kin.names[["AC"]]%in%low.kinases]),int.publi=log10(1+kin.int.tot[kin.names[["AC"]]%in%low.kinases])))
abline(a=reg2$coefficients[1],b=reg2$coefficients[2],col="green",lwd=2)
legend("bottomright",lwd=2,col=c("red","green"),legend=c(paste("adjusted R² =",format(summary(reg)$adj.r.squared, digits=2)),paste("adjusted R² =",format(summary(reg2)$adj.r.squared, digits=2),"(low abundant)")))
dev.off()


no.inter <- kin.names[kin.tot>=50 & kin.int.tot<1,"Kinase"]
one.inter <- kin.names[kin.tot>=50 & kin.int.tot<=1,"Kinase"]
no.pm <- kin.names[kin.tot<1 & kin.int.tot>=1,"Kinase"]
no.inter
one.inter
no.pm

# Simple based on pairs
vect.pub <- unlist(apply(kin.pub,1,list))
vect.int.pub <- unlist(apply(kin.int.pub,1,list))
no0 <- (vect.pub != 0) | (vect.int.pub != 0)
plot(x=log(1+vect.pub[no0]),y=log(1+vect.int.pub[no0]),pch=20)
reg <- lm(int.publi~publi,data.frame(publi=log(1+vect.pub[no0]),int.publi=log(1+vect.int.pub[no0])))
summary(reg)
abline(a=reg$coefficients[1],b=reg$coefficients[2],col="red",lwd=2)

t <- NULL
for (i in 1:dim(kin.pub)[1])
  t <- c(t,cumsum(kin.pub[i,]))
c.kin.pub <- matrix(t,nrow=dim(kin.pub)[1],byrow=T)
vect.pub <- t
t <- NULL
for (i in 1:dim(kin.int.pub)[1])
  t <- c(t,cumsum(kin.int.pub[i,]))
c.kin.int.pub <- matrix(t,nrow=dim(kin.int.pub)[1],byrow=T)
vect.int.pub <- t
no0 <- (vect.pub != 0) & (vect.int.pub != 0)
pdf("figures/time-correlation-acc-1.pdf",,width=3,height=3,pointsize=8)
#plot(x=log(1+vect.pub[no0]),y=log(1+vect.int.pub[no0]),pch=".")
symbols(x=log10(1+vect.pub[no0]),y=log10(1+vect.int.pub[no0]),circles=rep(0.015,sum(no0)),inches=F,bg="black",fg=NULL)
reg <- lm(int.publi~publi,data.frame(publi=log(1+vect.pub[no0]),int.publi=log(1+vect.int.pub[no0])))
#abline(a=reg$coefficients[1],b=reg$coefficients[2],col="red",lwd=2)
summary(reg)
legend("bottomright",bty="n",legend=paste("R²=",format(summary(reg)$adj.r.squared, digits=2)))
dev.off()

cr <- sample(rainbow(50),size=50,replace=F)
postscript("figures/time-correlation-acc-2.pdf",pointsize=10)
plot(x=1,y=1,type="n",xlim=c(0,4),ylim=c(0,3))
for (i in 1:dim(kin.pub)[1]){
  s <- min(kin.1sty[i],kin.int.1sty[i])
  if (s < dim(kin.pub)[2]){
    x0 <- log10(1+c.kin.pub[i,(s+1):dim(kin.pub)[2]])
    x <- c(rep(0,64-length(x0)),x0)#(x0-x0[1])#/(x0[length(x0)]-x0[1])
    y0 <- log10(1+c.kin.int.pub[i,(s+1):dim(kin.pub)[2]])
    y <- c(rep(0,64-length(y0)),y0)#y0-y0[1]
    lines(x=runif(64-s,min=-0.05,max=0.05)+x,
          y=runif(64-s,min=-0.05,max=0.05)+y,type="l",
          col=cr[1+i%%50])
  }
}
abline(a=reg$coefficients[1],b=reg$coefficients[2],col="black",lwd=4)
dev.off()

mat <- matrix(rep(0,150*150),nrow=150)
for (k in 1:sum(no0)){
  i <- 1+trunc(10*log(1+vect.int.pub[no0][k]))
  j <- 1+trunc(10*log(1+vect.pub[no0][k]))
  mat[i,j] <- mat[i,j]+1
}
image(log(1+mat),col=colors)

# Temporal correlation and classes of profiles
t <- NULL
for (i in 1:dim(kin.pub)[1])
  #  t <- c(t,cumsum(kin.pub[i,]),cumsum(kin.int.pub[i,]))
  t <- c(t,cumsum(kin.pub[i,26:64]),cumsum(kin.int.pub[i,26:64]))
c.both <- matrix(t,nrow=dim(kin.int.pub)[1],byrow=T)
rownames(c.both) <- kin.names[["Kinase"]]
#colnames(c.both) <- c(colnames(kin.pub),colnames(kin.pub))
colnames(c.both) <- c(colnames(kin.pub)[26:64],colnames(kin.pub)[26:64])
hmlog <- heatmap(c.both,col=colors,Colv=NA)
heatmap(log10(1+c.both),col=colors,Colv=NA,scale="none",Rowv=hmlog$rowInd)
pdf("figures/fig2-pub-history.pdf",pointsize=8,width=6,height=8)
png("figures/fig2-pub-history.png",pointsize=8,width=6,height=8,res=600,units="in")
hm <- heatmap(log10(1+c.both),col=colors,Colv=NA,scale="none",keep.dendro=T)
dev.off()

write.table(data.frame(kinase=kin.names[rev(hm$rowInd),"Kinase"],pub=c.both[rev(hm$rowInd),39],int=c.both[rev(hm$rowInd),78]),file="figures/fig2-kinase-list.txt",quote=F,sep="\t")
max(c.both)

fact <- max(c.kin.pub)/max(c.kin.int.pub)
c.both2 <- c.both
c.both2[,65:128] <- fact*c.both2[,65:128]
hmlog2 <- heatmap(c.both2,col=colors,Colv=NA)
heatmap(log(1+c.both2),col=colors,Colv=NA,scale="none")
heatmap(log(1+c.both2),col=colors,Colv=NA,scale="none",Rowv=hmlog2$rowInd)

dc <- dist(log10(1+c.both))
hi <- hclust(dc)
plot(hi)

n.clust <- 3
grow <- 1+c.both
km <- kmeans(log10(grow),centers=n.clust)
#cr <- c("red","blue","green","cyan","orange")#rainbow(n.clust)
cr <- c("red","green","yellow")
cr <- c(rgb(1,0,0,0.4),rgb(0,1,0,0.4),rgb(1,1,0,0.4))
pdf("figures/time-correlation-acc-2.pdf",width=3,height=3,pointsize=8,useDingbats=FALSE)
plot(x=grow[,1:39],y=grow[,40:78],type="n",log="xy",xlog=T,ylog=T)
rect(10^(par("usr")[1]),10^(par("usr")[3]),10^(par("usr")[2]),10^(par("usr")[4]),col="grey85")
grid(col="white",lwd=2,lty=1)
for (i in 1:dim(grow)[1]){
  lines(x=jitter(grow[i,1:39]),
        y=jitter(grow[i,40:78]),type="l",
        col=cr[km$cluster[i]])
}
cr2 <- c("red","green","yellow")
for (i in c(2,3,1)){
  lines(x=10^km$centers[i,1:39],
        y=10^km$centers[i,40:78],
        col="black",lwd=9)
  lines(x=10^km$centers[i,1:39],
        y=10^km$centers[i,40:78],
        col=cr2[i],lwd=5)
}
dev.off()

hm <- heatmap(log10(1+c.both),col=colors,Colv=NA,scale="none",RowSideColors=cr[km$cluster[kin.names[,"Kinase"]]])



# Publication content =================================

# Cancer
census <- read.csv("diseases//cancer_gene_census.txt",stringsAsFactors=F,check.names=F,sep="\t")
in.cen <- kin.names[["Kinase"]] %in% census[[1]]
length(unique(census[[1]]))
sum(in.cen)
cont.tab <- matrix(c(sum(kin.tot[in.cen]),sum(kin.tot[!in.cen]),sum(kin.int.tot[in.cen]),sum(kin.int.tot[!in.cen])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2]) # less cancer in PPI

cont.tab <- matrix(c(sum(kin.tot[in.cen]),sum(kin.tot[!in.cen]),sum(kin.int.tot[in.cen]),sum(kin.int.tot[!in.cen]),
                     sum(kin.bait.tot[in.cen]),sum(kin.bait.tot[!in.cen]),sum(kin.low.tot[in.cen]),sum(kin.low.tot[!in.cen])),ncol=4)
cont.tab
chisq.test(cont.tab[,1:2])
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2]) # less cancer in PPI
chisq.test(cont.tab[,c(1,3)])
cont.tab[1,1]*cont.tab[2,3]/(cont.tab[2,1]*cont.tab[1,3]) # less cancer in PPI
chisq.test(cont.tab[,c(1,4)])
cont.tab[1,1]*cont.tab[2,4]/(cont.tab[2,1]*cont.tab[1,4]) # less cancer in PPI

# KEGG
sp.kegg <- read.csv("/home/colinge/databases/uniprot_to_kegg_human-2012-08-23.csv",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c('character','character'))
kegg2names <- read.csv("/home/colinge/databases/kegg_names-2012-08-23.txt",check.names=F,sep="\t",stringsAsFactors=F,colClasses=c("character","character"))
kegg.names <- kegg2names[[2]]
names(kegg.names) <- kegg2names[[1]]
#sp.kegg <- sp.kegg[!is.na(sp.kegg[[2]]),]
#for (i in 1:dim(sp.kegg)[1]){
#  l <- nchar(sp.kegg[i,2])
#  sp.kegg[i,2] <- paste(paste(rep("0",5-l),collapse=""),sp.kegg[i,2],sep="")
#}

keggSets <- function(map,descr){  
  kegg <- NULL
  for (term in unique(map[[2]])){
    mem <- map[map[[2]]==term,1]
    inpw <- kin.names[["AC"]] %in% mem
    cont.tab <- matrix(c(sum(kin.tot[inpw]),sum(kin.tot[!inpw]),sum(kin.int.tot[inpw]),sum(kin.int.tot[!inpw])),ncol=2)
    pval <- chisq.test(cont.tab)$p.value
    if (!is.nan(pval)){
      effect <- cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2])
      de <- descr[as.character(term)]
      if (is.na(de))
        de <- term
      kegg <- rbind(kegg,data.frame(term=term,descr=de,p.val=pval,num.pw=length(mem),in.pw=sum(inpw),effect=effect,stringsAsFactors=F))
    }
  }
  kegg
}

kegg <- keggSets(sp.kegg,kegg.names)
kegg[kegg[["p.val"]]<1.0e-4,]

write.table(kegg[kegg[["effect"]]>1.5 & kegg[["p.val"]]<1.0e-4 & kegg[["in.pw"]]>=3,],file="figures/fig2-kegg-up-list.txt",quote=F,sep="\t")  # Cancer again
write.table(kegg[kegg[["effect"]]<1/1.5 & kegg[["p.val"]]<1.0e-4 & kegg[["in.pw"]]>=3,],file="figures/fig2-kegg-down-list.txt",quote=F,sep="\t")  # Almost no cancer

kegg.tab <- matrix(c(0,8,18,66),ncol=2)
chisq.test(kegg.tab)

k.cancer <- NULL
for (term in unique(sp.kegg[[2]]))
  if (length(intersect(c('cancer','leukemia','carcinoma','Glioma','Melanoma'),unlist(strsplit(kegg.names[term],split=" ")))) >= 1){
    cat(term,kegg.names[term],"\n")
    mem <- sp.kegg[sp.kegg[[2]]==term,1]
    #inpw <- intersect(kin.names[["AC"]],mem)    
    k.cancer <- c(k.cancer,mem)
  }
k.cancer <- unique(k.cancer)
length(k.cancer)
in.cen <- kin.names[["AC"]] %in% k.cancer
sum(in.cen)
cont.tab <- matrix(c(sum(kin.tot[in.cen]),sum(kin.tot[!in.cen]),sum(kin.int.tot[in.cen]),sum(kin.int.tot[!in.cen])),ncol=2)
cont.tab
chisq.test(cont.tab)
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2]) # less cancer in PPI

cont.tab <- matrix(c(sum(kin.tot[in.cen]),sum(kin.tot[!in.cen]),sum(kin.int.tot[in.cen]),sum(kin.int.tot[!in.cen]),
                     sum(kin.bait.tot[in.cen]),sum(kin.bait.tot[!in.cen]),sum(kin.low.tot[in.cen]),sum(kin.low.tot[!in.cen])),ncol=4)
cont.tab
chisq.test(cont.tab[,1:2])
cont.tab[1,1]*cont.tab[2,2]/(cont.tab[2,1]*cont.tab[1,2]) # less cancer in PPI
chisq.test(cont.tab[,c(1,3)])
cont.tab[1,1]*cont.tab[2,3]/(cont.tab[2,1]*cont.tab[1,3]) # less cancer in PPI
chisq.test(cont.tab[,c(1,4)])
cont.tab[1,1]*cont.tab[2,4]/(cont.tab[2,1]*cont.tab[1,4]) # less cancer in PPI



