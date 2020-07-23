#
# Prepare family structures and run Cri-MAP to extract crossover positions
# Susan Johnston
# April 2017
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

auto.vec <- 1:26    # 1:29
sex.vec <- 27
chr.vec <- 1:27

#~~ load old models

recsumm <- read.table("data/2_TotalSummIndivRR_FullCleanPostSim_g.txt", header = T, stringsAsFactors = F)
head(recsumm)

#~~ Load GenABEL gwaa.data for all genotyped sheep

load("data/Plates_1to87_QC3.GenABEL.RData")

#~~ Read in pedigree file

pedigree <- read.table("data/20190711_Soay_Pedigree.txt", header = T, stringsAsFactors = F)
names(pedigree)[1] <- "ANIMAL"
pedigree$MOTHER[is.na(pedigree$MOTHER)] <- 0
pedigree$FATHER[is.na(pedigree$FATHER)] <- 0

pedigree$MOTHER[grep("F", pedigree$MOTHER)] <- 0
pedigree$FATHER[grep("M", pedigree$FATHER)] <- 0
pedigree <- pedigree[-grep("F", pedigree$ANIMAL),]
pedigree <- pedigree[-grep("M", pedigree$ANIMAL),]

abeldata <- soay87
rm(soay87)


#~~ deal with wrong sexes

x <- phdata(abeldata)
x$id[which(x$sex == 0)][which(x$id[which(x$sex == 0)] %in% pedigree$FATHER)]
x$id[which(x$sex == 1)][which(x$id[which(x$sex == 1)] %in% pedigree$MOTHER)]

#phdata(abeldata[x$id[which(x$sex == 0)][which(x$id[which(x$sex == 0)] %in% pedigree$FATHER)],])$sex <- 1

#~~ read in SNP positions and PAR SNPs

mapdata <- read.table("data/2_Linkage_Map_Positions_g.txt", header = T, stringsAsFactors = F)
mapdata <- arrange(mapdata, Chr, Order)

mapdata <- subset(mapdata, SNP.Name %in% snpnames(abeldata))

head(mapdata)

pseudoautoSNPs <- readLines("data/1_Pseudoautosomal_SNPs_in_X.txt")

#~~ Define thresholds for pedigree construction

mend.err.locus.threshold <- 0.01
mend.err.pair.threshold  <- 0.001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Determine Working Pedigrees for Crimap    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove non-genotyped parents from the pedigree

pedigree    <- subset(pedigree, ANIMAL %in% idnames(abeldata))
pedigree$MOTHER[which(!pedigree$MOTHER %in% idnames(abeldata))] <- 0
pedigree$FATHER[which(!pedigree$FATHER %in% idnames(abeldata))] <- 0

head(pedigree)

#~~ Determine all P-O pairs where the parent has both parents and mate are genotyped.

po.pairs <- melt(pedigree, id.vars = "ANIMAL")
po.pairs <- subset(po.pairs, value != 0)

head(po.pairs)
po.pairs$Selected <- NA

famped <- NULL

for(i in 1:nrow(po.pairs)){

  if(i %in% seq(1, nrow(po.pairs), 500)) print(paste("Running line", i, "of", nrow(po.pairs)))
  
  #~~ Extract the focal parent
  
  x1 <- pedigree[which(pedigree$ANIMAL == po.pairs$value[i]),]
  
  #~~ Extract the offspring
  
  x2 <- pedigree[which(pedigree$ANIMAL == po.pairs$ANIMAL[i]),]
  
  #~~ Extract focal parent's parents and mate
  
  x3 <- data.frame(ANIMAL = c(unlist(x1[1,2:3]), ifelse(po.pairs$variable[i] == "MOTHER", x2$FATHER[1], x2$MOTHER[1])), 
                   MOTHER = 0, FATHER = 0)
  
  x4 <- rbind(x3, x1, x2)
  x4$Family <- paste("Offspring", po.pairs$ANIMAL[i], po.pairs$variable[i], po.pairs$value[i], sep = "_")
  
  x4 <- subset(x4, ANIMAL != 0)
  
  if(nrow(x4) == 5){
    po.pairs$Selected[i] <- "yes"
    famped <- rbind(famped, x4)
  } else {
    po.pairs$Selected[i] <- "no"
  }
  
  rm(x1, x2, x3, x4)
  
}

head(famped)
rm(i)

nrow(famped)/5

famped$ANIMAL <- as.character(famped$ANIMAL)

AnalysisSuffix <- "b"

sort(table(chromosome(abeldata)))

#~~ Crimap cannot handle chromosome 1 on my machine as it is too large. Split it into two:

head(mapdata)
mapdata$Chr2 <- mapdata$Chr
length(mapdata$Chr2[mapdata$Chr == 1])

mapdata$Chr2[mapdata$Chr == 1][1:2000] <- "1p"
mapdata$Chr2[mapdata$Chr == 1][2001:length(mapdata$Chr2[mapdata$Chr == 1])] <- "1q"

table(mapdata$Chr2)

chr.vec <- c("1p", "1q", 2:27)
sex.vec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create crimap files & Run 1st instance    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   is.X = ifelse(i == sex.vec, TRUE, FALSE),
                   pseudoautoSNPs = pseudoautoSNPs,
                   genabel.phdata = phdata(abeldata),
                   save.mendfile = TRUE)
    
  }
  
})

rm(po.pairs, x, i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Deal with Mendelian Errors                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

menderr <- NULL

for(i in chr.vec){
  
  temperr <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mndverbose"),
                        header = T, sep = "\t", stringsAsFactors = F) 
  
  if(i == sex.vec) temperr <- subset(temperr, select = -sex)
  
  menderr <- rbind(menderr,
                   cbind(temperr, analysisID = paste0(i, AnalysisSuffix)))
  
  rm(temperr)
  
}

#~~ Do particular pairs perform badly?

menderr <- subset(menderr, SNP.Name %in% snpnames(abeldata))

poor.ids <- subset(menderr, select = c(ANIMAL, MOTHER, Mat.Mismatch, Family))
poor.ids2 <- subset(menderr, select = c(ANIMAL, FATHER, Pat.Mismatch, Family))

names(poor.ids) <- c("ANIMAL", "PARENT", "Mismatch", "Family")
names(poor.ids2) <- c("ANIMAL", "PARENT", "Mismatch", "Family")

poor.ids$PARENT.Type <- "MOTHER"
poor.ids2$PARENT.Type <- "FATHER"

poor.ids <- rbind(poor.ids, poor.ids2)

rm(poor.ids2)

poor.ids <- droplevels(subset(poor.ids, Mismatch == "yes"))
poor.ids$ID.Parent.Family <- paste(poor.ids$ANIMAL, poor.ids$PARENT, poor.ids$Family, sep = "_")

poor.ids <- data.frame(table(poor.ids$ID.Parent))
names(poor.ids)[1] <- "ID.Parent"

ggplot(poor.ids, aes(Freq)) + geom_histogram(binwidth = 1, col = "grey")

table(poor.ids$Freq)

poor.ids <- poor.ids[which(poor.ids$Freq > 50),]

#~~ which families have these bad ids?

poor.ids$Family <- sapply(as.character(poor.ids$ID.Parent), function (x) paste(strsplit(x, split = "_")[[1]][3:6], collapse = "_"))

famped <- subset(famped, !Family %in% poor.ids$Family) #85 families removed

#~~ which SNPs are bad??

head(menderr)
menderr <- subset(menderr, !Family %in% poor.ids$Family)

table(table(menderr$SNP.Name))

write.table(poor.ids, paste0("results/1_ID_pairs_with_high_Mend_errors_", AnalysisSuffix, ".txt"),
            row.names = F, quote = F)

write.table(famped, 
            paste("results/1_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)




