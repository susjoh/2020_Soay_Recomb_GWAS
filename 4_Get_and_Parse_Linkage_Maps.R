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

source("2.1_Dat_file_function.R")

auto.vec <- 1:26    # 1:29
sex.vec <- 27
chr.vec <- 1:27

#~~ Add CriMAP requirements

AnalysisSuffix <- "b"

#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Plates_1to87_QC3.GenABEL.RData")
abeldata <-soay87
rm(soay87)

#~~ Read in pedigree files

pedigree <- read.table("data/20190711_Soay_Pedigree.txt", header = T, stringsAsFactors = F)
famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ read in SNP positions and PAR SNPs

mapdata <- read.table("data/2_Linkage_Map_Positions_g.txt", header = T, stringsAsFactors = F)
mapdata <- arrange(mapdata, Chr, Order)
mapdata <- subset(mapdata, SNP.Name %in% snpnames(abeldata))
head(mapdata)

pseudoautoSNPs <- readLines("data/1_Pseudoautosomal_SNPs_in_X.txt")

#~~ Crimap cannot handle chromosome 1 & 2 on my machine as it is too large. Split it into two:

mapdata$Chr2 <- mapdata$Chr
mapdata$Chr2[mapdata$Chr == 1][1:2000] <- "1p"
mapdata$Chr2[mapdata$Chr == 1][2001:length(mapdata$Chr2[mapdata$Chr == 1])] <- "1q"
mapdata$Chr2[mapdata$Chr == 2][1:2000] <- "2p"
mapdata$Chr2[mapdata$Chr == 2][2001:length(mapdata$Chr2[mapdata$Chr == 2])] <- "2q"


table(mapdata$Chr2)

chr.vec <- c("1p", "1q", "2p", "2q", 3:26)
sex.vec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Rerun without Singletons and short segments, FINAL  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
firstRun <- FALSE

if(firstRun == T){
  system.time({
    
    for(i in chr.vec){
      
      print(paste("Running chromosome", i, "of", length(chr.vec)))
      
      if(i %in% c("1p", "1q", "2p", "2q")){
        if(i %in% c("1p", "2p")){
          x111 <- which(mapdata$Chr2 == i)
          x111 <- c(x111, max(x111)+1:100)
        } else {
          x111 <- which(mapdata$Chr2 == i)
          x111 <- c(min(x111)-100:1, x111)
        }
        
        tempsnplist <- mapdata[x111, "SNP.Name"]
      } else {
        tempsnplist <- subset(mapdata, Chr2 == i)$SNP.Name
        
      }
      
      create_crimap_input(gwaa.data = abeldata,
                          familyPedigree = famped,
                          analysisID = paste0(i, AnalysisSuffix),
                          snplist = tempsnplist,
                          is.X = ifelse(i == sex.vec, TRUE, FALSE),
                          pseudoautoSNPs = pseudoautoSNPs,
                          outdir = paste0("crimap/crimap_", AnalysisSuffix),
                          clear.existing.analysisID = TRUE,
                          use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"))
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
      run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
    }
    
  })
  
  i = 27
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(i, AnalysisSuffix),
                      snplist = subset(mapdata, Chr2 == i)$SNP.Name,
                      is.X = ifelse(i == sex.vec, TRUE, FALSE),
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  parse_mend_err(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                 paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                 familyPedigree = famped, is.X = T, pseudoautoSNPs = pseudoautoSNPs,
                 genabel.phdata = phdata(abeldata))
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(i, AnalysisSuffix),
                      snplist = subset(mapdata, Chr2 == i)$SNP.Name,
                      is.X = ifelse(i == sex.vec, TRUE, FALSE),
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE, use.mnd = T)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
}

maptab <- NULL

for(i in c(chr.vec, 27)){
  
  x <- parse_map(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".map"))
  head(x)
  
  maptab <- rbind(maptab, x)
  
}


mapdata %>% group_by(Chr) %>% summarise(Max.cm = max(cMPosition)) %>% data.frame
maptab %>% group_by(analysisID) %>% summarise(Max.cm = max(cMPosition)) %>% data.frame


#~~ Fix chromosomes 1 & 2

write.table(maptab, "test.txt", row.names = F, sep = "\t", quote = F)





#~~~~ Get the cumulative information

maptab <- add_count(maptab, SNP.Name)
subset(maptab, n == 2)
maptab$Chr <- gsub("b", "", maptab$analysisID)
maptab$Chr <- gsub("p", "", maptab$Chr)
maptab$Chr <- gsub("q", "", maptab$Chr)


mapdata$OldOrder <- NA
mapdata$OldOrder[1] <- 1
for(i in 2:nrow(mapdata)) mapdata$OldOrder[i] <- ifelse(mapdata$Chr[i] != mapdata$Chr[i-1], 1, mapdata$OldOrder[i-1] + 1)

x <- subset(mapdata, select = c(OldOrder, SNP.Name))

head(x)
maptab <- join(maptab, x)

head(maptab)

maptab <- arrange(maptab, Chr, OldOrder)

maptab <- maptab[-which(maptab$n == 2 & is.na(maptab$Female.r)),]
maptab <- add_count(maptab, SNP.Name)

maptab <- maptab[-which(maptab$n == 2 & maptab$Order < 300),]

dput(names(maptab))

maptab <- maptab[,c("Order", "OldOrder", "SNP.Name", "analysisID", "Chr",
                    "cMPosition", "r", "cMdiff", 
                    "cMPosition.Female", "Female.r", "cMdiff.Female",
                    "cMPosition.Male", "Male.r", "cMdiff.Male")]

head(maptab)

maptab$cMPosition.NEW <- NA
maptab$cMPosition.NEW[1] <- 0
for(i in 2:nrow(maptab)) maptab$cMPosition.NEW[i] <- ifelse(maptab$Chr[i] != maptab$Chr[i-1], 0, maptab$cMPosition.NEW[i-1] + maptab$cMdiff[i-1])


maptab[which(maptab$cMPosition != maptab$cMPosition.NEW),] %>% head




maptab$cMPosition.Female.NEW <- NA
maptab$cMPosition.Female.NEW[1] <- 0
for(i in 2:nrow(maptab)) maptab$cMPosition.Female.NEW[i] <- ifelse(maptab$Chr[i] != maptab$Chr[i-1], 0, maptab$cMPosition.Female.NEW[i-1] + maptab$cMdiff.Female[i-1])

maptab$cMPosition.Male.NEW <- NA
maptab$cMPosition.Male.NEW[1] <- 0
for(i in 2:nrow(maptab)) maptab$cMPosition.Male.NEW[i] <- ifelse(maptab$Chr[i] != maptab$Chr[i-1], 0, maptab$cMPosition.Male.NEW[i-1] + maptab$cMdiff.Male[i-1])


maptab %>% group_by(Chr) %>% summarise(Max.cm = max(cMPosition),
                                       Max.cm.NEW = max(cMPosition.NEW),
                                       Max.cm.Female = max(cMPosition.Female),
                                       Max.cm.NEW.Female = max(cMPosition.Female.NEW),
                                       Max.cm.Male = max(cMPosition.Male),
                                       Max.cm.NEW.Male = max(cMPosition.Male.NEW)) %>% data.frame

maptab$cMPosition        <- maptab$cMPosition.NEW
maptab$cMPosition.Female <- maptab$cMPosition.Female.NEW
maptab$cMPosition.Male   <- maptab$cMPosition.Male.NEW

maptab$cMPosition.NEW         <- NULL
maptab$cMPosition.Female.NEW  <- NULL
maptab$cMPosition.Male.NEW    <- NULL

maptab$Order <- maptab$OldOrder
maptab$OldOrder <- NULL
maptab$analysisID <- NULL


ggplot(maptab, aes(Order, cMPosition)) + geom_point() + facet_wrap(~Chr, scales = "free")
ggplot(maptab, aes(Order, cMPosition.Female)) + geom_point() + facet_wrap(~Chr, scales = "free")
ggplot(maptab, aes(Order, cMPosition.Male)) + geom_point() + facet_wrap(~Chr, scales = "free")


head(mapdata)



write.table(maptab, "results/4_20200504_Full_Linkage_Map.txt", row.names = F, sep = "\t", quote = F)



