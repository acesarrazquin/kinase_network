#########################################################################################
# Script to filter Networkin KSI predictions comparing them with experimental ones.
#
# Based on "networkin_score": 
#
#########################################################################################

# open networkin predictions
nkin <- read.table("networkin_formated.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# open known KSI (for the moment I don't consider KSI from pathways (reactome, PID))
pelm <- read.table("phosphoelm/phosphoelm_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
psplus <- read.table("phosphositeplus/phosphositeplus_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pnets <- read.table("phosphonetworks/phosphonetworks_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
intdb <- read.table("interactdb/interactdb_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# check overlap
all.ksi <- data.frame(kinase=c(intdb[,1], pelm[,1], psplus[,1], pnets[,1]), 
                      substrate = c(intdb[,2], pelm[,2], psplus[,2], pnets[,2]))

uniq.ksi.keys <- unique(paste(all.ksi$kinase, all.ksi$substrate, sep="|"))
nkin.ksi.keys <- paste(nkin[, 1], nkin[, 2], sep="|")
nkin$key <- nkin.ksi.keys

nkin.safe <- nkin.ksi.keys %in% uniq.ksi.keys
nkin.true <- nkin[nkin.safe,]

# compare scores

#1. Try PCA with the three scores -> NOT GOOD IDEA

# pca <- prcomp(nkin[,c("networkin_score","netphorest_score","string_score")],scale=T)
# summary(pca)
# 
# symbols(x=pca$x[,1],y=pca$x[,2],circles=rep(0.03,dim(pca$x)[1]),inches=F,bg="grey",fg=NULL)
# symbols(x=pca$x[nkin.safe,1],y=pca$x[nkin.safe,2],circles=rep(0.04,sum(nkin.safe)),inches=F,bg="red",fg=NULL,add=T)
# # the overlap is a bit everywhere...

#2. Try selecting only the networking_score (in the literature they say that >1 is already relatively good)
nkin.true$ksi_type <- "experimental"
nkin$ksi_type  <- "computational"
exp_comp_ksi <- rbind(nkin.true, nkin)

p <- ggplot(exp_comp_ksi)
p <- p + geom_boxplot(aes(y=as.numeric(networkin_score), x=factor(ksi_type)))
p <- p + coord_cartesian(ylim=c(0,15))
p

# chose as threshold 50% of experimental (very conservative but otherwise too many!)
thres <- quantile(nkin.true$networkin_score, 0.5)
nkin.filt <- nkin[nkin$networkin_score >= thres, ]

# merge equal interactions with different phosphoposition, and output file
uniq.keys <- unique(nkin.filt$key)
netkin.ksi.final <- data.frame()

for (key in uniq.keys){
  print(key)
  key.df <- data.frame()
  ints <- nkin.filt[nkin.filt$key == key, ]

  phospho_positions <- paste(sort(unique(ints$phospho_positions)), 
                             collapse=",")
  uniq.df <- unique(ints[, c('source', 'target', 'source_name', 'target_name', 'PMIDs',
                             'dates', 'sources', 'type')])
  uniq.df$phospho_positions <- phospho_positions
  
  netkin.ksi.final <- rbind(netkin.ksi_final, uniq.df)
}
write.table(netkin.ksi.final, "networkin_filtered_ksi.txt", row.names=FALSE, quote=FALSE,
            sep="\t", col.names=TRUE)
netkin.ksi.final
  