
ExtractHaploSHAPEIT <- function(results.file, make.ref.alleles = T){
  
  message("Reading and formatting shapeit results file...")
  
  #~~ read in the results file
  
  shapeit.res <- read.table(results.file, skip=6, stringsAsFactors = F)
  
  shapeit.refalleles <- shapeit.res[,4:5]
  
  #~~ transpose data and retain only the haplotype information
  haplotypes.raw <- data.frame(t(shapeit.res[,10:ncol(shapeit.res)]))
  
  #~~ add ID information and tidy up col names
  
  shapeittab <- read.table(results.file, skip = 5, nrows = 1, comment.char = "")
  shapeittab <- t(shapeittab[,10:ncol(shapeittab)])
  shapeittab <- data.frame(shapeittab)

  shapeitmap <- read.table(results.file, skip=6, stringsAsFactors = F)[,1:3]
  
  haplotypes.raw <- cbind(shapeittab, haplotypes.raw)
  names(haplotypes.raw) <- c("ID", as.character(shapeitmap[,3]))
  
  #~~ convert to character
  
  for(i in 1:ncol(haplotypes.raw)) haplotypes.raw[,i] <- as.character(haplotypes.raw[,i])
  
  #~~ The first two characters (i.e. 0|0) is the phased genotype. 
  #   Retain these and substitute in the reference alleles (if specified).
  
  message("Extracting phased genotypes...")
  
  haplotypes.extract <- data.frame(matrix(NA,
                                          nrow = nrow(haplotypes.raw),
                                          ncol = ncol(haplotypes.raw)))
  
  names(haplotypes.extract) <- names(haplotypes.raw)
  haplotypes.extract$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    haplotypes.extract[,i] <- substr(haplotypes.raw[,i], 1, 3)
    
    if(make.ref.alleles == TRUE){
      haplotypes.extract[,i] <- gsub(0, as.character(shapeit.refalleles[i-1,1]), haplotypes.extract[,i])
      haplotypes.extract[,i] <- gsub(1, as.character(shapeit.refalleles[i-1,2]), haplotypes.extract[,i])
    }
  }
  
  
  #~~ extract the genotype probabilities
  
  message("Extracting genotype probabilities...")
  
  
  genotype.probs <- data.frame(matrix(NA,
                                      nrow = nrow(haplotypes.raw),
                                      ncol = ncol(haplotypes.raw)))
  names(genotype.probs) <- names(haplotypes.raw)
  genotype.probs$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    genotype.probs[,i] <- sapply(haplotypes.raw[,i], function(x) strsplit(x, split = ":")[[1]][2])
  }
  
  #~~ create a results table for the haplotypes
  
  message("Identifying haplotypes...")
  
  haplo.results <- data.frame(ID = haplotypes.extract$ID,
                              Haplo1 = NA,
                              Haplo2 = NA)
  
  for(i in 1:nrow(haplotypes.extract)){
    haplo.results$Haplo1[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 1, 1), collapse = ""))
    haplo.results$Haplo2[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 3, 3), collapse = ""))
  }
  
  #~~ melt into one column
  require(reshape2)
  haplo.results.1col <- melt(haplo.results, id.vars="ID")
  names(haplo.results.1col) <- c("ID", "HaploID", "Haplo")
  
  for(i in 1:ncol(haplo.results.1col)) haplo.results.1col[,i] <- as.character(haplo.results.1col[,i])
  
  return(list(haplotypes = haplo.results.1col,
              phased.genotypes = haplotypes.extract,
              genotype.probabilities = genotype.probs))
  message("...done.")
}
