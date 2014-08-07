require(rcdk)

clusterDrugsSmile <- function(drugs_file)
  drugs  <- read.csv(drugs_file, header = TRUE, sep = '\t')
  smiles  <- as.character(drugs$CanonicalSmiles)
  names(smiles)  <- as.character(drugs$CompoundName)

  sp <- get.smiles.parser()
  mols <- parse.smiles(smiles)
  fplist <- lapply(mols, get.fingerprint, type = 'standard')

  fpdist <- fp.sim.matrix(fplist, method = 'tanimoto')
  fpdist <- as.dist(1-fpdist)

  
  clus <- hclust(fpdist)

  plot(clus, labels = names(fplist), main = 'davis2011_compounds')
