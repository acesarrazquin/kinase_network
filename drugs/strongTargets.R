strongTargets <- function(dpiFile, threshold){

  outliers = 0
  minimum = 0
  # Read file
  dpi <- read.csv(dpiFile, sep = '\t', header = TRUE)
  dpi$strong_target = 'no'
  
  # Get unique drugs
  uniq.drugs <- as.character(unique(dpi[, 1]))

  # Loop over each drug
  for (i in 1:length(uniq.drugs)){
    print(uniq.drugs[i])
    drug.dpi <- dpi[dpi$source == uniq.drugs[i], ]  
    print(threshold)
    
    values <- as.numeric(drug.dpi$strength)
    log.values <- log10(values)
    log.values[is.nan(log.values)] <- -10
    
    # 1. consider strong targets if the values are outliers
    log.out <- boxplot.stats(log.values)$out
    indexes <- which(log.values %in% log.out)
    
    out <- values[indexes]
    out <- out[out <= threshold]
    print(length(out))
    if (length(out) != 0){
      dpi[dpi$source == uniq.drugs[i] & dpi$strength %in% out, ]$strong_target = 'yes'
      outliers = outliers +1
    }
    
    # 2. consider strong targets if not outliers BUT at least 5 times more than min 
    thres <- min(threshold, 5*min(values))
    
    #print(length(dpi[dpi$source == uniq.drugs[i] & dpi$strength <= thres, ]$strong_target))
    len = length(dpi[dpi$source == uniq.drugs[i] & dpi$strength <= thres, ]$strong_target)
    print(len)
    if (len != 0){
      dpi[dpi$source == uniq.drugs[i] & dpi$strength <= thres, ]$strong_target = 'yes'
      minimum = minimum +1
    }  
  }
  print(outliers)
  print(minimum)
  #outname = paste(dpiFile, '-filtered_both', sep = '')
  #write.table(dpi, outname, sep = '\t', na = '', quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(dpi)
}


strongest2Targets <- function(dpiFile){ # for each dataset separately, and for Bantscheff take care with the two datasets (cellculture and lysate...)
  
  # Read file
  dpi <- read.csv(dpiFile, sep = '\t', header = TRUE)
  
  # Choose strongest two targts per compound
  #strongest <- data.frame(source=as.character(), target1=as.character(),target2=as.character())
  strongest <- data.frame()
  unique.sources <- c("22037378", "22037377", "17721511")
  for (sources in unique.sources){
    unique.compounds <- unique(as.character(dpi[dpi$sources == sources, ]$source))
    for (compound in unique.compounds){
      print(compound)
      print(sources)
      targets <- dpi[dpi$source == compound & dpi$sources == sources,]
    
      strength.order <- order(targets$strength, decreasing=FALSE)
      ordered.targets <- targets[strength.order,]
    
      target1 <- ordered.targets[1,]
      
      n=2
      target2 <- ordered.targets[n,]
          
      while (!is.na(target2$target) & target1$target == target2$target){
          n <- n+1
          target2 <- ordered.targets[n,]
      }   
      
      if (as.character(target1$target_form) != ''){
        target1.name <- paste(target1$target_name, target1$target_form, sep='-')
      }else{
        target1.name <- target1$target_name
      }
      
      
      if (!(is.na(target2$target_form)) & as.character(target2$target_form) != ''){
        target2.name <- paste(target2$target_name, target2$target_form, sep='-')
      }else{
        target2.name <- target2$target_name
      }
      
      #strongest <- rbind(strongest, data.frame(source=compound, target1=target1.name, target2=target2.name, sources=sources))
      strongest <- rbind(strongest, target1)
      strongest <- rbind(strongest, target2)
    }
  }
  return(strongest)
}