




parse_dat_file <- function(x, familyids, snpids){
  
  #~~ Determine which lines are SNPs, familyIDs, and the phasing output
  
  x.snps <- which(x %in% snpids)
  x.family <- which(x %in% familyids)
  x.nchar <- matrix(sapply(x, nchar))
  x.cmplines <- which(x.nchar == length(snpids))
  
  #~~ Get rid of the initial unID'ed phasing
  
  x <- x[-x.cmplines[which(x.cmplines < min(x.family))]]
  
  #~~ Re-determine which lines are SNPs, familyIDs, and the phasing output
  
  x.snps <- which(x %in% snpids)
  x.family <- which(x %in% familyids)
  
  familyids <- x[x.family]
  snpids <- unique(x[x.snps])
  
  x.nchar <- matrix(sapply(x, nchar))
  x.cmplines <- which(x.nchar == length(snpids))
  
  #~~ Remove SNP and family lines and de-determine phasing output lines
  
  x <- x[-c(x.snps, x.family)]
  x.nchar <- matrix(sapply(x, nchar))
  x.cmplines <- which(x.nchar == length(snpids))
  
  #~~ Find out where there are gaps between phased lines
  
  recombframe <- data.frame(Lines = x.cmplines,
                            Gap = c(7, diff(x.cmplines)))
  
  recombframe$Family <- NA
  recombframe$Family[which(recombframe$Gap > 1)] <- familyids
  recombframe$Family <- na.locf(recombframe$Family)
  
  recombframe$Gap <- NULL
  
  recombframe$data <- x[recombframe$Lines]
  
  rm(x, x.nchar, x.snps, x.family)
  
  recombframe <- recombframe[-which(recombframe$data == paste(rep("X", length(snpids)), collapse = "")),]
  
  y <- data.frame(table(recombframe$Family))
  names(y) <- c("Family", "NoChromatids")
  
  recombframe <- join(recombframe, y)
  
  recombframe$analysisID <- paste0(i, AnalysisSuffix)
  
  recombframe <- subset(recombframe, NoChromatids == 1)
  head(recombframe)
  
  recombframe$Lines <- NULL
  recombframe$NoChromatids <- NULL
  
  recombframe$data <- gsub("X", "-", recombframe$data)
  
  recombframe$ANIMAL = NA
  recombframe$RecombCount = NA
  recombframe$parent = NA
  
  head(recombframe)
  
  #~~ Get the IDs
  
  recombframe$ANIMAL <- unlist(lapply(recombframe$Family, function(foo) strsplit(foo, split = "_")[[1]][2]))
  recombframe$parent <- unlist(lapply(recombframe$Family, function(foo) strsplit(foo, split = "_")[[1]][3]))
  recombframe$RRID <- unlist(lapply(recombframe$Family, function(foo) strsplit(foo, split = "_")[[1]][4]))
  
  
  #~~ get the number of recombinations for each individual
  
  recExtract <- function(CMPstring){
    x <- gsub("-", "", CMPstring)
    x <- strsplit(x, split = "")[[1]]
    if(length(x) > 2){
      x <- length(which(c("X", x) != c(x, "X")))-2
    } else {
      NA
    }
    x
  }
  
  
  
  recombframe$RecombCount <- unlist(lapply(recombframe$data, recExtract))
  
  recombframe$data2 <- gsub("-" , "", recombframe$data)
  recombframe$data2 <- gsub("c" , "", recombframe$data2)
  recombframe$data2 <- gsub(":" , "", recombframe$data2)
  
  recombframe$No.Inf.Loci <- nchar(recombframe$data2)
  
  recombframe <- subset(recombframe, select = -data2)
  
  recombframe$RecombCount <- as.numeric(recombframe$RecombCount)
  
  #~~ Get first and list informative positions
  
  recombframe$First.Inf.Order <- unlist(lapply(recombframe$data, InfLengthFunc))
  recombframe$Last.Inf.Order <-  unlist(lapply(recombframe$data, function(x) InfLengthFunc(x, position = "Last")))
  
  #~~ add the RRID information - the individual in which the recombination took place
  
  suppressMessages(recombframe <- join(recombframe, famped))
  
  #~~ add analysisID
  
  recombframe$UniqueID <- paste(recombframe$analysisID, recombframe$Family, sep = "_")
  
  recombframe$data2 <- NULL
  recombframe
  
}

recExtract <- function(CMPstring){
  x <- gsub("-", "", CMPstring)
  x <- strsplit(x, split = "")[[1]]
  if(length(x) > 2){
    x <- length(which(c("X", x) != c(x, "X")))-2
  } else {
    NA
  }
  x
}

