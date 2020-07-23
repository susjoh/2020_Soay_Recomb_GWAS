#
# Extract linkage maps and crossover information
# Susan Johnston
# April 2017
#

library(ggplot2)
library(reshape2)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)
library(zoo)

source("2.1_Dat_file_function.R")

auto.vec <- 1:26    # 1:29
sex.vec <- 27
chr.vec <- 1:27

#~~ Add CriMAP requirements

AnalysisSuffix <- "b"

firstRun <- TRUE

#~~ Load GenABEL gwaa.data for all genotyped deer


load("data/Plates_1to87_QC3.GenABEL.RData")
abeldata <-soay87
rm(soay87)

#~~ Read in pedigree files

pedigree <- read.table("data/20190711_Soay_Pedigree.txt", header = T, stringsAsFactors = F)
famped <- read.table(paste0("results/1_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

names(pedigree)[1] <- "ANIMAL"
pedigree$MOTHER[is.na(pedigree$MOTHER)] <- 0
pedigree$FATHER[is.na(pedigree$FATHER)] <- 0

#~~ read in SNP positions and PAR SNPs

mapdata <- read.table("data/2_Linkage_Map_Positions_g.txt", header = T, stringsAsFactors = F)
mapdata <- arrange(mapdata, Chr, Order)

mapdata <- subset(mapdata, SNP.Name %in% snpnames(abeldata))

head(mapdata)

pseudoautoSNPs <- readLines("data/1_Pseudoautosomal_SNPs_in_X.txt")

#~~ Crimap cannot handle chromosome 1 on my machine as it is too large. Split it into two:

mapdata$Chr2 <- mapdata$Chr

mapdata$Chr2[mapdata$Chr == 1][1:2000] <- "1p"
mapdata$Chr2[mapdata$Chr == 1][2001:length(mapdata$Chr2[mapdata$Chr == 1])] <- "1q"

table(mapdata$Chr2)

chr.vec <- c("1p", "1q", 2:26)
sex.vec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create CRI-Map Files and run maps/chrompic  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun == TRUE){
  
  system.time({
    
    for(i in chr.vec){
      
      print(paste("Running chromosome", i, "of", length(chr.vec)))
      
      create_crimap_input(gwaa.data = abeldata, 
                          familyPedigree = famped,
                          analysisID = paste0(i, AnalysisSuffix),
                          snplist = subset(mapdata, Chr2 == i)$SNP.Name, 
                          is.X = ifelse(i == sex.vec, TRUE, FALSE),
                          pseudoautoSNPs = pseudoautoSNPs,
                          outdir = paste0("crimap/crimap_", AnalysisSuffix),
                          clear.existing.analysisID = TRUE,
                          use.mnd = TRUE)
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
    }
    
  })
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse crossovers and double-recombinants       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rectab <- NULL
doub.xovers <- NULL

for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  #~~ Create vectors of the family and snp IDs to remove from .dat file
  
  familyids <- unique(famped$Family)
  snpids <- subset(mapdata, Chr2 == i)$SNP.Name
  
  #~~ Read in the .dat file
  
  x <- readLines(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".dat"))
  
  y <- parse_dat_file(x, familyids, snpids)
  
  rectab <- rbind(rectab, y)
  
  doub.xovers <- rbind(doub.xovers, check_double_crossovers(y))
  
  rm(x)
  
}

#~~ Merge with the span positions
chrom.map <- subset(mapdata, SNP.Name %in% snpnames(abeldata) & Chr %in% 1:26)
chrom.map$Order[1] <- 1
for(i in 2:nrow(chrom.map)) chrom.map$Order[i] <- ifelse(chrom.map$Chr2[i] == chrom.map$Chr2[i-1], chrom.map$Order[i-1] + 1, 1)

chrom.map$analysisID <- paste0(chrom.map$Chr2, AnalysisSuffix)
x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)

write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_raw_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_raw_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_raw_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)

# rectab <- read.table("results/2_Per_Chromosome_Recomb_raw_b.txt", header = T, sep = "\t", stringsAsFactors = F)
# doub.xovers <- read.table("results/2_Double_Xovers_raw_b.txt", header = T, sep = "\t", stringsAsFactors = F)
# chrom.map <- read.table("results/2_cM_Map_raw_b.txt", header = T, sep = "\t", stringsAsFactors = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. QC and deal with singleton double xovers        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(doub.xovers)
head(chrom.map)

#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)

#~~ SAVE FIGURE

ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + 
  geom_histogram(binwidth = 1) + 
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Span Distance (cM)",
       y = "Count")
ggsave("figs/2_Double_Crossover_Histogram_QC1.png", width = 6, height = 4.5)


#~~ Remove all singletons phased runs shorter than 10cM

bad.doubs <- rbind(subset(doub.xovers, Singleton == "yes"))
head(bad.doubs)

x <- subset(chrom.map, select = c(analysisID, Order, SNP.Name))
names(x) <- c("analysisID", "StartPos", "SNP.Name")

bad.doubs <- join(bad.doubs, x)

head(bad.doubs)
bad.doubs <- subset(bad.doubs, select = c(SNP.Name, Family))
bad.doubs$RRID <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][4])
bad.doubs$ANIMAL <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][2])
bad.doubs$Family <- NULL

bad.doubs <- melt(bad.doubs, id.vars = "SNP.Name")
bad.doubs$variable <- NULL
names(bad.doubs) <- c("SNP.Name", "ANIMAL")

#~~ Make a master .mnd file and add bad doubles (singletons)

mnd.master <- NULL

for(i in chr.vec){
  mnd.file <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mnd"), header = T)
  mnd.master <- rbind(mnd.master, mnd.file)
  rm(mnd.file)
}

mnd.master <- rbind(mnd.master, bad.doubs)
mnd.master <- unique(mnd.master)

write.table(mnd.master, paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Rerun without Singletons                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- NULL
doub.xovers <- NULL

system.time({
  
  for(i in chr.vec){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    create_crimap_input(gwaa.data = abeldata, 
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        snplist = subset(mapdata, Chr2 == i)$SNP.Name, 
                        is.X = ifelse(i == sex.vec, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    #~~ Create vectors of the family and snp IDs to remove from .dat file
    
    familyids <- unique(famped$Family)
    snpids <- subset(mapdata, Chr2 == i)$SNP.Name
    
    #~~ Read in the .dat file
    
    x <- readLines(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".dat"))
    
    y <- parse_dat_file(x, familyids, snpids)
    
    rectab <- rbind(rectab, y)
    
    doub.xovers <- rbind(doub.xovers, check_double_crossovers(y))
    
    rm(x)
    
    
  }
  
})



#~~ Merge with the span positions

chrom.map <- subset(mapdata, SNP.Name %in% snpnames(abeldata) & Chr %in% 1:26)
chrom.map$Order[1] <- 1
for(i in 2:nrow(chrom.map)) chrom.map$Order[i] <- ifelse(chrom.map$Chr2[i] == chrom.map$Chr2[i-1], chrom.map$Order[i-1] + 1, 1)

chrom.map$analysisID <- paste0(chrom.map$Chr2, AnalysisSuffix)
x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)


write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_QC1_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_QC1_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_QC1_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. QC and deal with singleton double xovers        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(doub.xovers)
head(chrom.map)

#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)

#~~ SAVE FIGURE

ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + 
  geom_histogram(binwidth = 1) + 
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Span Distance (cM)",
       y = "Count")
ggsave("figs/2_Double_Crossover_Histogram_QC1.png", width = 6, height = 4.5)


#~~ Remove all singletons phased runs shorter than 10cM

bad.doubs <- rbind(subset(doub.xovers, Singleton == "yes"),
                   subset(doub.xovers, SpancMDist < 10))
head(bad.doubs)

x <- subset(chrom.map, select = c(analysisID, Order, SNP.Name))
names(x) <- c("analysisID", "StartPos", "SNP.Name")

bad.doubs <- join(bad.doubs, x)

head(bad.doubs)
bad.doubs <- subset(bad.doubs, select = c(SNP.Name, Family))
bad.doubs$RRID <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][4])
bad.doubs$ANIMAL <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][2])
bad.doubs$Family <- NULL

bad.doubs <- melt(bad.doubs, id.vars = "SNP.Name")
bad.doubs$variable <- NULL
names(bad.doubs) <- c("SNP.Name", "ANIMAL")


mnd.master <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), header = T, stringsAsFactors = F)
mnd.master <- unique(rbind(mnd.master, bad.doubs))


write.table(mnd.master, paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Rerun without Singletons and short segments, FINAL  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rectab <- NULL
doub.xovers <- NULL

system.time({

  for(i in chr.vec){

    print(paste("Running chromosome", i, "of", length(chr.vec)))

    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        snplist = subset(mapdata, Chr2 == i)$SNP.Name,
                        is.X = ifelse(i == sex.vec, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE,
                        use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"))

    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))

    #~~ Create vectors of the family and snp IDs to remove from .dat file

    familyids <- unique(famped$Family)
    snpids <- subset(mapdata, Chr2 == i)$SNP.Name

    #~~ Read in the .dat file

    x <- readLines(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".dat"))

    y <- parse_dat_file(x, familyids, snpids)

    rectab <- rbind(rectab, y)

    doub.xovers <- rbind(doub.xovers, check_double_crossovers(y))

    rm(x)


  }

})



#~~ Merge with the span positions

chrom.map <- subset(mapdata, SNP.Name %in% snpnames(abeldata) & Chr %in% 1:26)
chrom.map$Order[1] <- 1
for(i in 2:nrow(chrom.map)) chrom.map$Order[i] <- ifelse(chrom.map$Chr2[i] == chrom.map$Chr2[i-1], chrom.map$Order[i-1] + 1, 1)

chrom.map$analysisID <- paste0(chrom.map$Chr2, AnalysisSuffix)
x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)

write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_QC2_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_QC2_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. Extract and format crossover data                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab      <- read.table(paste0("results/2_Per_Chromosome_Recomb_QC2_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)
doub.xovers <- read.table(paste0("results/2_Double_Xovers_QC2_", AnalysisSuffix, ".txt")         , header = T, stringsAsFactors = F)
chrom.map <- read.table(paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt")               , header = T, stringsAsFactors = F)



#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)


#~~ Create new famped!

famped <- subset(famped, Family %in% rectab$Family)


#~~ remove the remaining short crossovers directly from the rectab object

doub.xovers <- rbind(subset(doub.xovers, Singleton == "yes"),
                   subset(doub.xovers, SpancMDist < 10))

for(i in 1:nrow(doub.xovers)){
  
  row.n <- which(rectab$UniqueID == doub.xovers$UniqueID[i])
  if(length(row.n) > 0){
    x <- rectab[row.n,"data"]
    x <- strsplit(x, split = "")[[1]]
    x[doub.xovers$StartPos[i]:doub.xovers$StopPos[i]] <- "-"
    
    #~~ revise crossover number
    
    x1 <- x[which(x %in% c(0, 1))]
    rectab$RecombCount[row.n] <- length(which(c(-999, x1) != c(x1, -999)))-2
    
    rectab$data[row.n] <- paste(x, collapse = "")
    
    rm(x, x1, row.n)
  }
}

head(rectab)

rectab$No.Inf.Loci <- sapply(rectab$data, function(x) nchar(gsub("-", "", x)))

rectab$Chr <- gsub(AnalysisSuffix, "", rectab$analysisID)

#~~ Deal with IDs where there were no informative loci

head(rectab)
table(is.na(rectab$No.Inf.Loci))


#~~ Deal with chromosome 1

x <- subset(rectab, Chr %in% c("1p", "1q"))

x1 <- data.frame(No.Inf.Loci = tapply(x$No.Inf.Loci, x$Family, sum),
                 First.Inf.Order = tapply(x$First.Inf.Order, x$Family, min),
                 Last.Inf.Order = tapply(x$Last.Inf.Order, x$Family, max))
x1$Family <- row.names(x1)

x2 <- subset(x, select = c( Family, Chr, data))
x2 <- dcast(x2, Family ~ Chr)
head(x2)
x2$data <- paste0(x2$`1p`, x2$`1q`)
x2$Family <- as.character(x2$Family)

x1 <- join(x1, x2)
x1 <- subset(x1, select = -c(`1p`, `1q`))

x <- subset(x, select = -c(No.Inf.Loci, First.Inf.Order, Last.Inf.Order, Chr, UniqueID, data, RecombCount, analysisID)) %>% unique
x <- join(x, x1)
x$analysisID <- paste0(1, AnalysisSuffix)
x$UniqueID <- paste0(x$analysisID, "_", x$Family)

names(rectab)[which(!names(rectab) %in% names(x))]
x$Chr <- gsub(AnalysisSuffix, "", x$analysisID)

x$RecombCount <- unlist(lapply(x$data, recExtract))

rectab <- subset(rectab, !Chr %in% c("1p", "1q"))
rectab <- rbind(rectab, x)

rm(x, x1, x2)

x <- data.frame(table(rectab$Family))
x <- subset(x, Freq < 26)

rectab <- subset(rectab, !Family %in% x$Var1)
famped <- subset(famped, !Family %in% x$Var1)

write.table(rectab, paste0("results/2_Per_Chromosome_Recomb_Final_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(famped, paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)

#~~ Create the basic recombination phenotypes

rectab <- read.table(paste0("results/2_Per_Chromosome_Recomb_Final_", AnalysisSuffix, ".txt"), header = T, sep = "\t")
rectab <- subset(rectab, Chr != 27)

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalChrCount    = tapply(rectab$RecombCount, rectab$Family, length),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum),
                      NonRecombCount   = tapply(rectab$RecombCount, rectab$Family, function(foo) length(which(foo == 0))))
recsumm$Family <- row.names(recsumm)

names(rectab)

x <- subset(rectab, select = c(ANIMAL, parent, Family, RRID)) %>% unique

recsumm <- join(recsumm, x)

recsumm

head(recsumm)

write.table(recsumm, paste0("results/2_Per_Individual_Recomb_Final_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)

#~~ remove the .gen files

setwd(paste0("crimap/crimap_", AnalysisSuffix))
x <- dir()
#x <- x[-grep(".cmp", x)]
x <- x[-grep(".mnd", x)]
x <- x[-grep(".loc", x)]

for(i in x) system(paste0("rm ", i))

