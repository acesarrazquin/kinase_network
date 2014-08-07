name2uniprotid <- function(namesFile){
  outname <- paste(strsplit(namesFile, '\\.')[[1]][1], '-ids.txt', sep = '')
  # read file
  genes <- as.character(read.csv(namesFile, header=FALSE, sep='\t')[,1])
                        
  # connecto to database
  require(RPostgreSQL)
  drv <- dbDriver("PostgreSQL")
  con  <- dbConnect(drv, dbname="cemmdev", host="dev1", user="acesarrazquin", password="#4acesarrazquin")
  
  gene_ids = data.frame(uniprot.id=as.character())
  for (gene in genes){
    myQuery <- 'SELECT a.accession FROM biodb.dbentry a, biodb.gene_name b, biodb.dbentry_gene_name c WHERE c.gene_name_id=b.gene_name_id AND c.dbentry_id=a.dbentry_id AND (a.database_id=3 OR a.database_id=4) AND a.species_id=6 AND b.name=$1'
    uniprot_ids <- dbGetQuery(con, myQuery, gene) # will have isoforms separated by a -
    uniprot_id <- unique(sapply(uniprot_ids, function(x) strsplit(x, '-')[[1]][1])) # will always get something?
    
    if (nrow(uniprot_ids) == 0){
      cat('UniProt ID not found!')
      uniprot_id = ''
    }
    print(gene)
    print(uniprot_id)
    gene_ids = rbind(gene_ids, data.frame(uniprot.id = uniprot_id))
    
    write.table(gene_ids, outname, quote=FALSE, row.names=FALSE, col.names=FALSE)
    
  }
  
  dbDisconnect(con)
  #dbUnloadDriver(drv)
  return(gene_ids)
}