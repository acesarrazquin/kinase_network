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
deg <- degree(alldb.G)
plot(density(sqrt(deg)))
lines(density(sqrt(deg[kin.in.net])),col="green")
write.table(data.frame(kinases,common=kinases[["UniProt_ID"]]%in%high.prot,stringsAsFactors=F),quote=F,row.names=F,sep="\t",file="figures/kinase-annot.txt")
rownames(kinases) <- kinases[["UniProt_ID"]]

multhist(list(all=common,kinases=common[kinases[["UniProt_ID"]]]),freq=F,col=c("black","blue"),breaks=11,ylab="Probability density",xlab="Num CL")

low.prot <- intersect(n.alldb.G,names(common)[common<=3])
high.prot <- intersect(n.alldb.G,names(common)[common==11])
low.kin <- intersect(low.prot,kin.in.net)
high.kin <- intersect(high.prot,kin.in.net)
pdf("figures/fig-3-kin-degree.pdf",width=3,height=2,pointsize=8,useDingbats=F)
plot(density(sqrt(deg)),main="",xlab=expression(sqrt("degree")),type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(sqrt(deg)))
lines(density(sqrt(deg[high.prot])),col="green")
lines(density(sqrt(deg[kin.in.net])),col="red")
#lines(density(sqrt(deg[low.prot])),col="orange")
#lines(density(sqrt(deg[low.kin])),col="cyan")
lines(density(sqrt(deg[high.kin])),col="orange")
legend(x="topright",lwd=1,legend=c("all proteins","common proteins","kinases","common kinases"),col=c("black","green","red","orange"),bg="white")
dev.off()

# Kinases connect to themselves =======================================================

getNeighb <- function(g,ACs,u=TRUE){
  mapped <- ACs[ACs %in% nodes(g)]
  if (u)
    unique(unlist(adj(g,mapped)))
  else
    unlist(adj(g,mapped))
  
} # getNeighb

closestFreq <- function(freq,deg){
  freq.deg <- as.numeric(names(freq))
  if (deg <= freq.deg[1])
    freq[1]
  else
    if (deg >= freq.deg[length(freq)])
      freq[length(freq)]
  else{
    i <- 2
    while (deg > freq.deg[i])
      i <- i+1
    freq[i-1]+(freq[i]-freq[i-1])*(deg-freq.deg[i-1])/(freq.deg[i]-freq.deg[i-1])
  }
} # closestFreq

# through PPIs --------------------

kin.neighb <- getNeighb(alldb.G,kin.in.net,u=F)
kin.neighb.in <- kin.neighb %in% kin.in.net
kin.neighb.out <- !kin.neighb.in
len.in <- sum(kin.neighb.in)/length(kin.in.net) # 4.013
len.out <- sum(kin.neighb.out)/length(kin.in.net) # 27.84
u.len.in <- length(intersect(kin.neighb,kin.in.net))
u.len.out <- length(setdiff(kin.neighb,kin.in.net))

dd.e <- table(deg[kin.in.net])
dd.w <- table(deg[setdiff(n.alldb.G,kin.in.net)])
dens.whole <- deg[setdiff(n.alldb.G,kin.in.net)]
for (i in 1:length(dens.whole))
  dens.whole[i] <- closestFreq(dd.e,dens.whole[i])/dd.w[as.character(dens.whole[i])]

l.in <- l.out <- u.l.in <- u.l.out <- NULL
for (r in 1:1000){
  rkin.ac <- sample(setdiff(n.alldb.G,kin.in.net),length(kin.in.net),prob=dens.whole)
  rkin.neighb <- getNeighb(alldb.G,rkin.ac,u=F)
  rin <- rkin.neighb %in% rkin.ac
  l.in <- c(l.in,sum(rin)/length(rkin.ac))
  l.out <- c(l.out,sum(!rin)/length(rkin.ac))
  u.l.in <- c(u.l.in,length(intersect(rkin.neighb,rkin.ac)))
  u.l.out <- c(u.l.out,length(setdiff(rkin.neighb,rkin.ac)))
}
write.table(data.frame(l.in,l.out,u.l.in,u.l.out),file="random-in-out.txt",sep="\t",quote=F,row.names=F,col.name=F)
len.inout <- read.csv(file="random-in-out.txt",header=F,sep="\t")
l.in <- len.inout[[1]]
l.out <- len.inout[[2]]
u.l.in <- len.inout[[3]]
u.l.out <- len.inout[[4]]

#pdf("figures/fig-3-self-kin-ppi-indiv.pdf",width=1.5,height=2,pointsize=8,useDingbats=F)
len.rat <- len.in/len.out
l.rat <- l.in/l.out
plot(density(l.rat),xlim=c(min(c(l.rat,len.rat)),max(c(l.rat,len.rat))),main="individual kinases",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(l.rat))
abline(v=len.rat,col="red")
sum(l.rat>len.rat)/length(l.rat) # 0
#dev.off()

pdf("figures/fig-3-self-kin-ppi-global.pdf",width=1.5,height=2,pointsize=8,useDingbats=F)
len.rat <- u.len.in/u.len.out
l.rat <- u.l.in/u.l.out
plot(density(l.rat),xlim=c(min(c(l.rat,len.rat)),max(c(l.rat,len.rat))),main="global coverage",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(l.rat))
abline(v=len.rat,col="red")
sum(l.rat>len.rat)/length(l.rat) # 0
dev.off()

# Export kinase-kinase PPI network -----------
kikin <- alldb[alldb[["source"]]%in%kin.in.net & alldb[["target"]]%in%kin.in.net,]
write.table(kikin,quote=F,row.names=F,sep="\t",file="figures/kinase-kinase-ppi.txt")

# Export kinase-kinase PPI & KSI network
full <- read.csv("full-network-20140211.txt",sep="\t",check.names=FALSE,stringsAsFactors=F)
f.no.self <- as.vector(full[["source"]]) != as.vector(full[["target"]])
kikin2 <- full[full[["source"]]%in%kinases[["UniProt_ID"]] & full[["target"]]%in%kinases[["UniProt_ID"]] & f.no.self,]
write.table(data.frame(kikin2[,c("source","target","type")],graph.weight=rep(1,dim(kikin2)[1]),stringsAsFactors=F),quote=F,row.names=F,sep="\t",file="figures/kinase-kinase-all.txt")
# keys <- paste(kikin2[["source"]],kikin2[["target"]])
# gr.weight <- 0.05
# target.type <- kinases[kikin2[["target"]],"Type"]
# for (r in 1:3000){
#   s <- sample(kikin2[["source"]],1)
#   sc <- kinases[s,"Type"]
#   t <- sample(setdiff(kikin2[target.type==sc,"target"],s),1)
#   key <- paste(s,t)
#   if (!key%in%keys){
#     cat(s,"\t",t,"\tlayout\t",gr.weight,"\n",file="figures/kinase-kinase-all.txt",append=T,sep="")
#     keys <- c(keys,key)
#   }
# }



# Through KSIs ------------------------------------
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

kin.w.subs <- intersect(kinases[["UniProt_ID"]],ksi[["source"]])
kin.s.neighb <- getNeighb(ksi.G,kin.w.subs,u=F)
kin.s.neighb.in <- kin.s.neighb %in% kinases[["UniProt_ID"]]
kin.s.neighb.out <- !kin.neighb.in
len.s.in <- sum(kin.s.neighb.in)/length(kin.w.subs) # 4.013
len.s.out <- sum(kin.s.neighb.out)/length(kin.w.subs) # 27.84
su.len.in <- length(intersect(kin.s.neighb,kinases[["UniProt_ID"]]))
su.len.in
1-phyper(q=su.len.in,m=dim(kinases)[1],n=21400-dim(kinases)[1],k=length(unique(kin.s.neighb)))
su.len.out <- length(setdiff(kin.s.neighb,kinases[["UniProt_ID"]]))

# Level of control through substrates ==============================================
pdf("figures/fig-3-subs-degree.pdf",width=3,height=2,pointsize=8,useDingbats=F)
plot(density(sqrt(deg)),xlab=expression(sqrt("degree")),type="n",main="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(sqrt(deg)))
lines(density(sqrt(deg[intersect(n.alldb.G,ksi[["target"]])])),col="blue")
lines(density(sqrt(deg[intersect(n.alldb.G,ksi[ksi[["type"]]=="KSI","target"])])),col="cyan")
#lines(density(sqrt(deg[intersect(n.alldb.G,ksi[ksi[["source"]]%in%high.prot,"target"])])),col="red")
lines(density(sqrt(deg[intersect(n.alldb.G,ksi[ksi[["source"]]%in%low.prot,"target"])])),col="violet")
legend(x="topright",lwd=1,legend=c("all proteins","kinase substrates","kinase subs. exp.","common substrates"),col=c("black","blue","cyan","violet"),bg="white")
dev.off()

n.kin.s <- length(intersect(kinases[["UniProt_ID"]],ksi[["target"]]))
n.kin.s.exp <- length(intersect(kinases[["UniProt_ID"]],ksi[ksi[["type"]]=="KSI","target"]))
n.kin.high.s <- length(intersect(kinases[["UniProt_ID"]],ksi[ksi[["source"]]%in%high.kin,"target"]))
n.kin.low.s <- length(intersect(kinases[["UniProt_ID"]],ksi[ksi[["source"]]%in%low.kin,"target"]))

subs.in.net <- intersect(n.alldb.G,ksi[["target"]])
dd.e <- table(deg[subs.in.net])
dd.w <- table(deg[setdiff(n.alldb.G,subs.in.net)])
dens.whole <- deg[setdiff(n.alldb.G,subs.in.net)]
for (i in 1:length(dens.whole))
  dens.whole[i] <- closestFreq(dd.e,dens.whole[i])/dd.w[as.character(dens.whole[i])]

nks <- NULL
for (r in 1:1000){
  rsubs.ac <- sample(setdiff(n.alldb.G,subs.in.net),length(subs.in.net),prob=dens.whole)
  nks <- c(nks,length(intersect(rsubs.ac,kin.in.net))/length(subs.in.net))
}
write.table(nks,file="random-nks.txt",sep="\t",quote=F,row.names=F,col.name=F)
nks <- read.csv(file="random-nks.txt",header=F,sep="\t")[[1]]

pdf("figures/fig-3-subs-kin-ksi-global.pdf",width=1.5,height=2,pointsize=8,useDingbats=F)
plot(density(nks),xlim=c(min(c(nks,n.kin.s/length(subs.in.net))),max(c(nks,n.kin.s/length(subs.in.net)))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(nks))
abline(v=n.kin.s/length(subs.in.net),col="red")
abline(v=n.kin.s.exp/length(subs.in.net),col="red",lty=2)
dev.off()

num.subs <- length(unique(ksi[[2]]))
num.ksi <- dim(ksi)[1]
num.kin <- length(unique(ksi[[1]]))
kin.cover <- sum(ksi[[2]]%in%kin.in.net)/num.ksi
kin.cover.exp <- sum(ksi[ksi[["type"]]=="KSI","target"]%in%kin.in.net)/sum(ksi[["type"]]=="KSI")
nks <- NULL
for (r in 1:1000){
  rsubs.ac <- sample(setdiff(n.alldb.G,subs.in.net),num.subs,prob=dens.whole)
  assigned <- sample(rsubs.ac,num.ksi,replace=T)
  num.cover <- sum(assigned%in%kin.in.net)
  nks <- c(nks,num.cover/num.ksi)
}
write.table(nks,file="random-nks-indiv.txt",sep="\t",quote=F,row.names=F,col.name=F)
nks <- read.csv(file="random-nks-indiv.txt",header=F,sep="\t")[[1]]

pdf("figures/fig-3-subs-kin-ksi-indiv.pdf",width=1.5,height=2,pointsize=8,useDingbats=F)
plot(density(nks),xlim=c(min(c(nks,kin.cover)),max(c(nks,kin.cover.exp))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(nks))
abline(v=kin.cover,col="red")
abline(v=kin.cover.exp,col="red",lty=2)
dev.off()



# Interconnectivity between kinase families ===========================

# By PPIs --------------------
fam.size <- table(kinases[["Type"]])
kin.class <- kinases[["Type"]]
names(kin.class) <- kinases[["UniProt_ID"]]
classes <- unique(kinases[["Type"]])
n.cl <- length(classes)
thres <- 1

class.net <- NULL
for (i in 1:length(classes)){
  n.drawn <- sum(kin.class[kikin[["source"]]]==classes[i] | kin.class[kikin[["target"]]]==classes[i])
  for (j in i:length(classes)){
    n.links <- sum((kin.class[kikin[["source"]]]==classes[i] & kin.class[kikin[["target"]]]==classes[j]) |
                     (kin.class[kikin[["source"]]]==classes[j] & kin.class[kikin[["target"]]]==classes[i]))
    pval <- 1-phyper(q=n.links,m=fam.size[classes[j]],n=dim(kinases)[1]-fam.size[classes[j]],k=n.drawn)
    if ((pval <= thres) && (n.links > 1))
      class.net <- rbind(class.net,data.frame(from=classes[i],to=classes[j],pval=pval,lpval=-log10(pval),num=n.links,stringsAsFactors=F))
  }
}
class.net[is.infinite(class.net[["lpval"]]),"lpval"] <- 2*max(class.net[!is.infinite(class.net[["lpval"]]),"lpval"])
write.table(class.net,quote=F,row.names=F,sep="\t",file="figures/class-class-net.txt")
write.table(data.frame(class=classes,class.size=fam.size[classes],stringsAsFactors=F),file="figures/class-nodes.txt",quote=F,row.names=F,sep="\t")

# By partners (1-edge PPIs)
neighb <- adj(alldb.G,kin.in.net)
partners <- list(AGC=c(),CAMK=c(),CK1=c(),CMGC=c(),NEK=c(),RGC=c(),STE=c(),TKL=c(),Tyr=c(),Other=c(),Atypical=c())
for (k in names(neighb)){
  cl <- kinases[k,"Type"]
  partners[[cl]] <- union(partners[[cl]],neighb[[k]])
}
n.partners <- lapply(partners,length)
n.balls <- length(n.alldb.G)
n.balls <- length(unique(unlist(partners)))
p.class.net <- NULL
for (i in 1:(length(classes)-1)){
  n.drawn <- n.partners[[classes[i]]]
  for (j in (i+1):length(classes)){
    n.white <- n.partners[[classes[j]]]
    s.overlap <- length(intersect(partners[[classes[i]]],partners[[classes[j]]]))
    pval <- 1-phyper(q=s.overlap,m=n.white,n=n.balls-n.white,k=n.drawn)
    if ((pval <= thres) && (s.overlap >= 5))
      p.class.net <- rbind(p.class.net,data.frame(from=classes[i],to=classes[j],pval=pval,lpval=-log10(pval),num=s.overlap,stringsAsFactors=F))
  }
}
p.class.net[is.infinite(p.class.net[["lpval"]]),"lpval"] <- 2*max(p.class.net[!is.infinite(p.class.net[["lpval"]]),"lpval"])
write.table(p.class.net,quote=F,row.names=F,sep="\t",file="figures/class-class-1step-net.txt")

# By KSIs --------------------
thres <- 1

k.class.net <- NULL
for (cli in classes){
  drawn <- intersect(ksi[kin.class[ksi[["source"]]]==cli,"target"],kinases[["UniProt_ID"]])
  n.drawn <- length(drawn)
  for (clj in classes){
    n.links <- sum(kin.class[drawn]%in%clj,na.rm=T)
    pval <- 1-phyper(q=n.links,m=fam.size[clj],n=dim(kinases)[1]-fam.size[clj],k=n.drawn)
    if ((pval <= thres) && (n.links > 1))
      k.class.net <- rbind(k.class.net,data.frame(from=cli,to=clj,pval=pval,lpval=-log10(pval),num=n.links,stringsAsFactors=F))
  }
}
k.class.net[is.infinite(k.class.net[["lpval"]]),"lpval"] <- 2*max(k.class.net[!is.infinite(k.class.net[["lpval"]]),"lpval"])
write.table(k.class.net,quote=F,row.names=F,sep="\t",file="figures/class-class-ksi-net.txt")

# By substrates
n.balls <- length(unique(ksi[["target"]]))
s.class.net <- NULL
for (cli in classes){
  drawn <- unique(ksi[kin.class[ksi[["source"]]]==cli,"target"])
  n.drawn <- length(drawn)
  for (clj in classes){
    if (cli != clj){
      whites <- unique(ksi[kin.class[ksi[["source"]]]==clj,"target"])
      n.white <- length(whites)
      s.overlap <- length(intersect(whites,drawn))
      pval <- 1-phyper(q=s.overlap,m=n.white,n=n.balls-n.white,k=n.drawn)
      if ((pval <= thres) && (s.overlap >= 5))
        s.class.net <- rbind(s.class.net,data.frame(from=cli,to=clj,pval=pval,lpval=-log10(pval),num=s.overlap,stringsAsFactors=F))
    }
  }
}
s.class.net[is.infinite(s.class.net[["lpval"]]),"lpval"] <- 2*max(s.class.net[!is.infinite(s.class.net[["lpval"]]),"lpval"])
write.table(s.class.net,quote=F,row.names=F,sep="\t",file="figures/class-class-1step-ksi-net.txt")

