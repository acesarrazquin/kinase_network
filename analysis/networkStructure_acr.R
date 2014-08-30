
library(SparseM)
library(RBGL)
library(Rgraphviz)
library(plotrix)

# PPI network ----------------------
alldb <- read.csv("ppi-20140211.txt",
                  sep="\t", check.names=FALSE, stringsAsFactors=F)
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

# Mann's 11 cell line core proteomes ------------------
core <- read.csv("inputfiles/core-prots.csv", sep="\t", stringsAsFactors=F)
present <- NULL
names <- NULL
for (j in seq(11, 41, by=3)){ #IBAQ values
  pres <- (core[[j]] > 0) | (core[[j+1]] > 0) | (core[[j+2]] > 0) # triplicates
  present <- cbind(present,pres)
  names <- c(names,strsplit(names(core)[j], split="[._]", perl=T)[[1]][2])
}
colnames(present) <- names
present[is.na(present)] <- FALSE
num.cl <- apply(present, 1, sum) #counts in how many tissues the protein is present

#sp.acs <- read.csv("/home/colinge/haplogene/human-sp-acs.txt",sep="\t",stringsAsFactors=F)[[1]]
up.acs <- read.csv("inputfiles/up2name_human-20140812.tab", 
                   sep="\t", stringsAsFactors=F)[[1]]
acs <- strsplit(core[["Uniprot"]]," ")
acs.up <- lapply(acs, function(x){intersect(up.acs, gsub("-.+","",x))[1]}) #only intersects first elemnt of list
num.names <- unlist(acs.up)

good <- !is.na(num.names)
common <- num.cl[good]
names(common) <- num.names[good]

# Human kinases --------------------------------------------
kinases <- read.csv("inputfiles/UniProt_human_protein_kinases_formated-atypical-1class.txt",
                    sep="\t", stringsAsFactors=F)
kin.in.net <- intersect(kinases[["UniProt_ID"]], n.alldb.G)

deg <- degree(alldb.G)
# plot(density(sqrt(deg)))
# lines(density(sqrt(deg[kin.in.net])),col="green")
# write.table(data.frame(kinases,common=kinases[["UniProt_ID"]]%in%high.prot,stringsAsFactors=F),quote=F,row.names=F,sep="\t",file="figures/kinase-annot.txt")
# rownames(kinases) <- kinases[["UniProt_ID"]]
# 
# multhist(list(all=common,kinases=common[kinases[["UniProt_ID"]]]),freq=F,col=c("black","blue"),breaks=11,ylab="Probability density",xlab="Num CL")

low.prot <- intersect(n.alldb.G, names(common)[common<=3])
high.prot <- intersect(n.alldb.G, names(common)[common==11])
low.kin <- intersect(low.prot, kin.in.net)
high.kin <- intersect(high.prot, kin.in.net)

#pdf("figures/fig-3-kin-degree.pdf",width=3,height=2,pointsize=8,useDingbats=F)

deg.df <- data.frame(degrees=c(deg, deg[high.prot], deg[kin.in.net], deg[high.kin]),
                     types=c(rep(1, length(deg)), 
                              rep(2, length(deg[high.prot])), 
                              rep(3, length(deg[kin.in.net])), 
                              rep(4, length(deg[high.kin]))))

p <- ggplot(deg.df) + geom_line(aes(x=sqrt(degrees), color=factor(types)), 
                                   stat='density',
                                size=2.5)
p <- p + xlim(c(-0.5, 20)) 
p <- p + xlab(expression(sqrt('PPI degree'))) + labs(color="")
p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.position='none',
               
               legend.text=element_text(size=20, color="grey10", vjust=1))
p <- p + scale_color_manual(breaks=c(1,2,3,4),
                              values=c("#6E6E6E", "#0DBF0D", "#FF000D", "#FF9900"),
                              labels=c("all proteins", "common proteins",
                                       "kinases", "common kinases"))

p
pdf("analysis/kin-degree.pdf", width=12, height=8)
print(p)
dev.off()

# substrates
deg.kin.sub  <- deg[intersect(n.alldb.G, ksi[["target"]])]
deg.kin.sub.exp <- deg[intersect(n.alldb.G,ksi[ksi[["type"]]=="KSI","target"])]
#deg.kin.sub.comm <- deg[intersect(n.alldb.G,ksi[ksi[["source"]]%in%high.prot,"target"])]
deg.kin.sub.comm <- deg[intersect(n.alldb.G,ksi[ksi[["source"]]%in%low.prot,"target"])]

deg.sub.df <- data.frame(degrees=c(deg, deg.kin.sub, deg.kin.sub.exp, deg.kin.sub.comm),
                         types=c(rep(1, length(deg)), 
                                 rep(2, length(deg.kin.sub)), 
                                 rep(3, length(deg.kin.sub.exp)), 
                                 rep(4, length(deg.kin.sub.comm))))

p <- ggplot(deg.sub.df) + geom_line(aes(x=sqrt(degrees), color=factor(types)), 
                                stat='density',
                                size=2.5)
p <- p + xlim(c(-0.5, 20)) 
p <- p + xlab(expression(sqrt('PPI degree'))) + labs(color="")
p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.key=element_blank(),
               
               legend.text=element_text(size=20, color="grey10", vjust=1))
p <- p + scale_color_manual(breaks=c(1,2,3,4),
                            values=c("#6E6E6E", "#3F18FF", "#6DE0FF", "#B30FE8"),
                            labels=c("all proteins", "kinase substrates",
                                     "kinase subs. exp.", "common substrates"))

p

pdf("analysis/kin-subs-degree.pdf", width=12, height=8)
print(p)
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

kin.neighb <- getNeighb(alldb.G, kin.in.net, u=F)
kin.neighb.in <- kin.neighb %in% kin.in.net
kin.neighb.out <- !kin.neighb.in

len.in <- sum(kin.neighb.in)/length(kin.in.net) # 4.013
len.out <- sum(kin.neighb.out)/length(kin.in.net) # 27.84

u.len.in <- length(intersect(kin.neighb, kin.in.net))
u.len.out <- length(setdiff(kin.neighb, kin.in.net))

dd.e <- table(deg[kin.in.net])
dd.w <- table(deg[setdiff(n.alldb.G, kin.in.net)])
dens.whole <- deg[setdiff(n.alldb.G, kin.in.net)] #degree of nonkinases
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

# individual
len.rat <- len.in/len.out
l.rat <- l.in/l.out
l.rat.df <- data.frame(val=l.rat)

p <- ggplot(l.rat.df) + geom_line(aes(val), stat="density", size=2, color="#6E6E6E")
p <- p + geom_vline(xintercept=len.rat, color="#FF000D", size=2.5)
p <- p + xlab(expression('average individual PPI ratio'~sum(K[i])/sum(N[i]))) + labs(color="")
p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.key=element_blank(),
               
               legend.text=element_text(size=20, color="grey10", vjust=1))

pdf("analysis/self-kin-ppi-indiv.pdf", width=12, height=6)
print(p)
dev.off()

## global

len.rat <- u.len.in/u.len.out
l.rat <- u.l.in/u.l.out
l.rat.df <- data.frame(val=l.rat)
l.rat.df.fake <- data.frame(val=len.rat)

p <- ggplot(l.rat.df) + geom_line(aes(val), stat="density", size=2, color="#6E6E6E")
#p <- ggplot(l.rat.df.fake) + geom_point(aes(val,val)) 
p <- p + geom_vline(xintercept=len.rat, color="#FF000D", size=2.5)
p <- p + xlab('ratio kinases/non-kinases PPIs (K/N)') + labs(color="") +
         ylab("density")
p <- p + xlim(c(0.0475, len.rat+0.025*len.rat)) +ylim(c(0,210))
p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.key=element_blank(),
               
               legend.text=element_text(size=20, color="grey10", vjust=1))

pdf("analysis/self-kin-ppi-glob.pdf", width=12, height=6)
print(p)
dev.off()

plot(density(l.rat),xlim=c(min(c(l.rat,len.rat)),max(c(l.rat,len.rat))),main="global coverage",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(l.rat))
abline(v=len.rat,col="red")
sum(l.rat>len.rat)/length(l.rat) # 0
dev.off()

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

1-phyper(q=su.len.in,m=dim(kinases)[1],n=21400-dim(kinases)[1],k=length(unique(kin.s.neighb)))
su.len.out <- length(setdiff(kin.s.neighb,kinases[["UniProt_ID"]]))

# Level of control through substrates ==============================================

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
# write.table(nks,file="random-nks.txt",sep="\t",quote=F,row.names=F,col.name=F)
# nks <- read.csv(file="random-nks.txt",header=F,sep="\t")[[1]]

nks.df <- data.frame(val=nks)

intercept <- n.kin.s/length(subs.in.net)
intercept.exp <- n.kin.s.exp/length(subs.in.net)

nks.df.fake <- data.frame(val=intercept)
p <- ggplot(nks.df) + geom_line(aes(val), stat="density", size=2, color="#6E6E6E")
#p <- ggplot(nks.df.fake) + geom_point(aes(val,val)) 
p <- p + geom_vline(xintercept=intercept, color="#FF000D", size=2.5)
p <- p + geom_vline(xintercept=intercept.exp, color="#FF000D", size=2.5,
                                                  linetype=6)
p <- p + xlab('ratio kinases/non-kinases KSIs (K/N)') + labs(color="") +
  ylab("density")
p <- p +ylim(c(0,310))
p <- p + scale_x_continuous(breaks=c(0.02, 0.04, 0.06, 0.08, 0.10),
                            labels=c('0.02', "0.04", '0.06', "0.08", ""),
                            limits=c(0.0125, intercept+0.025*intercept))
p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.key=element_blank(),
               
               legend.text=element_text(size=20, color="grey10", vjust=1))

pdf("analysis/self-kin-ksi-glob.pdf", width=12, height=6)
print(p)
dev.off()

#individual
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
# write.table(nks,file="random-nks-indiv.txt",sep="\t",quote=F,row.names=F,col.name=F)
# nks <- read.csv(file="random-nks-indiv.txt",header=F,sep="\t")[[1]]

nks.df <- data.frame(val=nks)

intercept <- kin.cover
intercept.exp <- kin.cover.exp

p <- ggplot(nks.df) + geom_line(aes(val), stat="density", size=2, color="#6E6E6E")
p <- p + geom_vline(xintercept=intercept, color="#FF000D", size=2.5)
p <- p + geom_vline(xintercept=intercept.exp, color="#FF000D", size=2.5, linetype=6)
p <- p + xlab(expression('average individual KSI ratio'~sum(K[i])/sum(N[i]))) + labs(color="")

p <- p + scale_x_continuous(breaks=c(0.02, 0.04, 0.06, 0.08, 0.10, 0.12),
                            labels=c('0.02', "", '0.06', "", '0.10', ""),
                            limits=c(0.01, intercept.exp+0.025*intercept.exp))

p <- p + theme(panel.grid.minor=element_blank(),
               axis.title=element_text(size=30, color="grey10"),
               axis.text=element_text(size=26, color="grey10"),
               axis.title.x=element_text(vjust=-1),
               axis.title.y=element_text(vjust=1.5),
               axis.line=element_line(color="grey10"),
               axis.ticks=element_line(color="grey10"),
               legend.key=element_blank(),
               
               legend.text=element_text(size=20, color="grey10", vjust=1))

pdf("analysis/self-kin-ksi-indiv.pdf", width=12, height=6)
print(p)
dev.off()



nks.df.fake <- data.frame(val=intercept)
plot(density(nks),xlim=c(min(c(nks,kin.cover)),max(c(nks,kin.cover.exp))),main="",type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey85")
grid(col="white",lwd=1,lty=1)
lines(density(nks))
abline(v=kin.cover,col="red")
abline(v=kin.cover.exp,col="red",lty=2)
dev.off()


