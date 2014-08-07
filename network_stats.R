##########################################################################
# Script to calculate statistics of the interaction network 
#
# E.g. how many proteins from each source, etc, overlaps...
#
#
##########################################################################

netCountsAbs <- function(infile, type){
  options(stringsAsFactors=TRUE)
  fullnet <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        check.names=FALSE)
  
  dpi.mapping <- data.frame(name=c('drugbank', 'davis2011', 'anastassiadis2011', 'bantscheff2007',
                                   'rix2013', 'remsing2009', 'rix2007', 'rix2010'),
                            pid=c('drugbank', '22037378', '22037377', '17721511',
                                  '24130846', '19039322', '17720881', '19890374'))
  
  
  int.type <- fullnet[fullnet$type == type, ]
  
  if (type=="KSI"){
    sources <- c("networkin", "phospho.elm","phosphositeplus","phosphonetworks","intact",
                 "pid(NCI-Nature_Curated)","pid(Reactome)","reactome","PID","dip","innatedb","mint",
                 "hprd","corum","Behrends2010" )
    
  }else if (type=="PPI"){
    sources <- c("intact","mint","innatedb","dip","matrixdb","corum","hprd",
                 "PID","pid(BioCarta)","pid(NCI-Nature_Curated)","pid(Reactome)","reactome","Behrends2010")
    
    
  }else if (type=="DPI"){
    options(stringsAsFactors=FALSE)
    sources <- c('drugbank', '22037378', '22037377', '17721511',
                     '24130846', '19039322', '17720881', '19890374')
    #sources <- unique(strsplit(paste(int.type$sources, collapse=","), split=",")[[1]])
  }

  df.type <- data.frame()
  for (src in sources){

    overlaps <- vector()
    
    # for the parenthesis of pid...
    src_ <- gsub("\\(", "\\\\(", src)
    src_ <- gsub("\\)", "\\\\)", src_)
    
    int.src <- int.type[grep(src_, int.type$sources),]

    total <- nrow(int.src)
    names(total) <- 'total'

    
    int.self <- nrow(int.src[int.src$sources==src, ])
    names(int.self) <- src

    
    print(paste(src, ": ", total, " (unique: ", int.self, ")", sep=""))
    
    overlaps <- c(overlaps, int.self)
    
    index <- which(sources==src)
    sources2 <- sources[-index]
    for (src2 in sources2){
      
      src2_ <- gsub("\\(", "\\\\(", src2)
      src2_ <- gsub("\\)", "\\\\)", src2_)
      
      int.pair <- nrow(int.src[grep(src2_, int.src$sources), ])
      names(int.pair) <- src2
      
      overlaps <- c(overlaps, int.pair)
      
    }
    overlaps <- c(overlaps, total)
    
    df.type <- rbind(df.type, data.frame(t(overlaps)))
    
  }
  
  
  if (type == "DPI"){  
     colnames(df.type) <- c(dpi.mapping$name, 'total')
     row.names(df.type)  <- dpi.mapping$name
  }else{
    row.names(df.type) <- sources
  }
  return(df.type)
}

netCountsRel <- function(infile, type){
  options(stringsAsFactors=TRUE)
  fullnet <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        check.names=FALSE)
  
  dpi.mapping <- data.frame(name=c('drugbank', 'davis2011', 'anastassiadis2011', 'bantscheff2007',
                                   'rix2013', 'remsing2009', 'rix2007', 'rix2010'),
                            pid=c('drugbank', '22037378', '22037377', '17721511',
                                  '24130846', '19039322', '17720881', '19890374'))
  
  

  int.type <- fullnet[fullnet$type == type, ]
  total.type <- nrow(int.type)
  if (type=="KSI"){
    sources <- c("networkin", "phospho.elm","phosphositeplus","phosphonetworks","intact",
                 "pid(NCI-Nature_Curated)","pid(Reactome)","reactome","PID","dip","innatedb","mint",
                 "hprd","corum","Behrends2010" )
    
  }else if (type=="PPI"){
    sources <- c("intact","mint","innatedb","dip","matrixdb","corum","hprd",
                 "PID","pid(BioCarta)","pid(NCI-Nature_Curated)","pid(Reactome)","reactome","Behrends2010")
    
  }else if (type=="DPI"){
    options(stringsAsFactors=FALSE)
    sources <- c('22037377', '22037378', '24130846', '19039322', '17720881', 
                 '17721511', '19890374', 'drugbank')
    #sources <- unique(strsplit(paste(int.type$sources, collapse=","), split=",")[[1]])
  }
  
  df.type <- data.frame(db1=c(), db2=c(), overlap=c())
  for (src in sources){
       
    # for the parenthesis of pid...
    src_ <- gsub("\\(", "\\\\(", src)
    src_ <- gsub("\\)", "\\\\)", src_)
    
    int.src <- int.type[grep(src_, int.type$sources),]
    total <- nrow(int.src)

    int.self <- nrow(int.src[int.src$sources==src, ])
    df.type <- rbind(df.type, data.frame(db1=src, db2=src, overlap=int.self/total))
    
    index <- which(sources==src)
    sources2 <- sources[-index]
    for (src2 in sources2){
      
      src2_ <- gsub("\\(", "\\\\(", src2)
      src2_ <- gsub("\\)", "\\\\)", src2_)
      
      int.pair <- nrow(int.src[grep(src2_, int.src$sources), ])
      df.type <- rbind(df.type, data.frame(db1=src, db2=src2, overlap=int.pair/total))
      
    }
    #df.type <- rbind(df.type, data.frame(db1=src, db2='total', ov=total))
  }
  
  if (type == "DPI"){
    for (pid in dpi.mapping$pid){
      db.name <- as.character(dpi.mapping[dpi.mapping$pid==pid, 'name'])

      df.type[as.character(df.type$db1) == pid, 'db1'] <- db.name
      df.type[as.character(df.type$db2) == pid, 'db2'] <- db.name
    }
  }
  
  #df.type <- rbind(df.type, data.frame(db1='total', db2='total', ov=total.type))

  return(df.type)
}

plotOverlap <- function(infile, type="PPI"){
  dataset <- netCountsRel(infile, type)
  
  
  p <- ggplot(dataset, aes(y=db1, x=db2))
  p <- p + geom_tile(aes(fill = overlap), colour = "white") 
  
  p <- p + scale_fill_gradient(low = "white", high = "steelblue", 
                               limits=c(0,1))
  p <- p + theme(panel.grid=element_blank(), panel.background=element_blank())
  
  p <- p + theme(axis.text.y=element_text(size=10),
                 axis.text.x=element_text(angle=45, hjust=1, size=10))
  p
}