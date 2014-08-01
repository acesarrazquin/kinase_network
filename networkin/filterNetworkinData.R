

# open networkin predictions
nkin <- read.table("networkin_formated.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# open known KSI (for the moment I don't consider KSI from pathways (reactome, PID))
pelm <- read.table("../phosphoelm/phosphoelm_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
psplus <- read.table("../phosphositeplus/phosphositeplus_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pnets <- read.table("../phosphonetworks/phosphonetworks_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
intdb <- read.table("../interactdb/interactdb_ksi.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# check overlap
all.ksi <- data.frame(kinase=c(intdb[,1], pelm[,1], psplus[,1], pnets[,1]), 
                      substrate = c(intdb[,2], pelm[,2], psplus[,2], pnets[,2]))
uniq.ksi <- unique(paste(all.ksi$kinase, all.ksi$substrate, sep="|"))
nkin.ksi <- paste(nkin[, 1], nkin[, 2], sep="|")