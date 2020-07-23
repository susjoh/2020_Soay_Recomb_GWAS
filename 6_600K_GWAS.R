#
# GWAS for meiosis traits
#
#
#

library(GenABEL)
library(RepeatABEL)
library(magrittr)
library(plyr)
library(ggplot2)
library(asreml)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load Data and Functions                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("r/rGLSadj.R")

#~~ Read in recombination data and prepare for "cis/trans" and "trans only"

recsumm <- read.table("results/3_20200721_Individual_Recombination_Rates.txt", header = T, stringsAsFactors = F)

rectab <- read.table("results/3_20200721_Individual_Recombination_Rates_per_Chr.txt", header = T, stringsAsFactors = F)
names(rectab)[which(names(rectab) == "ANIMAL")] <- "Offspring.ID"

recsumm_trans <- NULL

for(i in 1:26){
  
  x <- subset(rectab, Chr == i)
  x <- subset(x, select = c(RRID, Offspring.ID, RecombCount))
  x$Chromosome <- i
  
  x <- join(recsumm, x)
  head(x)
  x$TotalRecombCountTrans <- x$TotalRecombCount - x$RecombCount
  
  recsumm_trans <- rbind(recsumm_trans, x)
  
}

rm(rectab)

recsumm_trans <- subset(recsumm_trans, select = c(RRID, Offspring.ID, Chromosome, TotalRecombCountTrans))

#~~ Load the genomic data and add the X chromosome information

# load("../../Soay Sheep Genomic Data/20200517_imputation/sheep_geno_imputed_oar31_17052020.GenABEL.RData")
# load("data/Plates_1to87_QC3.GenABEL.RData")
# 
# snp50 <- snpnames(soay87)
# snpHD <- snpnames(soayimp_oar3)
# snp50X <- snpnames(soay87[,chromosome(soay87) == 27])
# 
# #~~ Add SNPs that weren't imputed (includes the entire X chromosome)
# 
# soay87 <- soay87[,which(!snpnames(soay87) %in% snpnames(soayimp_oar3))]
# soayimp <- merge.gwaa.data(soayimp_oar3, soay87)
# rm(soay87, soayimp_oar3, snp50, snpHD, snp50X)
# save(soayimp, file = "soayimp_genotype_data.RData")

load("soayimp_genotype_data.RData")

#~~ Format pedigree information

pedigree <- read.table("data/20190711_Soay_Pedigree.txt", header = T, stringsAsFactors = F)
names(pedigree)[1] <- "RRID"

pedigree$MOTHER[grep("F", pedigree$MOTHER)] <- NA
pedigree$FATHER[grep("M", pedigree$FATHER)] <- NA
pedigree <- pedigree[-grep("F", pedigree$RRID),]
pedigree <- pedigree[-grep("M", pedigree$RRID),]

recsumm <- join(recsumm, pedigree)

recsumm$RRID <- factor(recsumm$RRID)
recsumm$MOTHER <- factor(recsumm$MOTHER)

#~~ Load GRM

# grm.auto <- read.table("../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Plates_1to87_QC3.GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
# ids.auto <- read.table("../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Plates_1to87_QC3.GRM.adj.grm.id")  # CONTAINS ID LIST
# 
# grminv <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$RRID) # vector of IDs from the datasset that you use for the asreml model
# dim(grminv)
# save(grminv, file = "data/GRM.RData")

load("data/GRM.RData")

grminv2 <- grminv
grminv2[lower.tri(grminv2)] = t(grminv2)[lower.tri(grminv2)]


#~~ Get minor allele frequency data

mafinfo <- summary.snp.data(gtdata(soayimp))
mafinfo$SNP.Name <- row.names(mafinfo)
mafinfo <- subset(mafinfo, select = c(Q.2, SNP.Name))
names(mafinfo)[1] <- "MAF"

#~~ Get Map info (relative to Rambouillet)

mapdata <- data.frame(Chr = chromosome(soayimp),
                      Position = map(soayimp),
                      SNP.Name = snpnames(soayimp))
mapdata$Chr <- as.numeric(as.character(mapdata$Chr))

mapdata <- arrange(mapdata, Chr, Position)
mapdata$Chr <- as.character(mapdata$Chr)
mapdata$SNP.Name <- as.character(mapdata$SNP.Name)
mapdata$Diff <- c(0, diff(mapdata$Position))
mapdata$Diff[which(mapdata$Diff < 0)] <- 1000
mapdata$Cumu <- cumsum(mapdata$Diff)

#~~ Get Chromosome Info

chrinfo <- NULL

for(i in na.omit(unique(mapdata$Chr))){
  
  temp1 <- arrange(subset(mapdata, Chr == i), Cumu)
  
  temp2 <- data.frame(Chr = i,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)


mapdata <- subset(mapdata, select = c(SNP.Name, Cumu))

rm(x, i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run the GWASs                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

firstRun <- FALSE

if(firstRun){
  
  # Both Sexes

  recsumm <- subset(recsumm, select = c(RRID, RRID.SEX, RRID.Fhat3, TotalRecombCount)) %>% na.omit
  recsumm <- subset(recsumm, RRID %in% dimnames(grminv2)[1][[1]])
  
  
  recsumm.m <- subset(recsumm, RRID.SEX == "M") %>% droplevels
  recsumm.f <- subset(recsumm, RRID.SEX == "F") %>% droplevels
  
  bothrec <- rGLSadj(formula.FixedEffects = TotalRecombCount ~ RRID.SEX + RRID.Fhat3,
                      genabel.data = soayimp,
                      phenotype.data = recsumm,
                      id = "RRID",
                      GRM = grminv2)
  
  bothrec <- process_rGLSadj_results(bothrec, soayimp)
  
  # Males
  
  
  malerec <- rGLSadj(formula.FixedEffects = TotalRecombCount ~ RRID.Fhat3,
                     genabel.data = soayimp,
                     phenotype.data = recsumm.m,
                     id = "RRID",
                     GRM = grminv2)
  
  malerec <- process_rGLSadj_results(malerec, soayimp)
  
  
  # Females
  
  femalerec <- rGLSadj(formula.FixedEffects = TotalRecombCount ~ RRID.Fhat3,
                       genabel.data = soayimp,
                       phenotype.data = recsumm.f,
                       id = "RRID",
                       GRM = grminv2)
  
  femalerec <- process_rGLSadj_results(femalerec, soayimp)
  
  
  save(bothrec, malerec, femalerec, file = "results/gwas_recrate.RData")
  rm(bothrec, malerec, femalerec)
  
  gc()
  
  
}

load("results/gwas_recrate.RData")

gwas.results <- rbind(cbind(bothrec, Sex = "Both Sexes"),
                      cbind(malerec, Sex = "Males"),
                      cbind(femalerec, Sex = "Females"))

rm(bothrec, malerec, femalerec)

#~~ Get MAF information

idvec <- unique(recsumm$RRID) %>% as.character
idvec <- idvec[which(idvec %in% idnames(soayimp))]

snpinfo <- summary.snp.data(gtdata(soayimp[idvec,]))
head(snpinfo)

snpinfo$SNP.Name <- row.names(snpinfo)
snpinfo <- subset(snpinfo, select = c(SNP.Name, Q.2))
head(snpinfo)

gwas.results <- join(gwas.results, snpinfo)
table(gwas.results$Chromosome)
gwas.results$Chromosome <- as.numeric(as.character(gwas.results$Chromosome))

gwas.results <- join(gwas.results, mapdata)

gwas.results <- gwas.results[,c("Sex", "SNP.Name", "Chromosome", "Position", "A1", "A2", "effB", "se_effB", 
                                "chi2.1df", "P1df", "Pc1df",  "ExpP",  
                                "Q.2", "Cumu")]
gwas.results$Sex <- as.character(gwas.results$Sex)

write.table(gwas.results, file = "results/6_GWAS_Results_600K.txt", row.names = F, sep = "\t", quote = F)

gwas.results <- subset(gwas.results, Pc1df < 1.28e-6)

ggplot(gwas.results, aes(Cumu, -log10(Pc1df), colour = Chromosome %% 2)) +
  geom_point() +
  facet_wrap(~Sex) +
  geom_hline(yintercept = -log10(1.28e-6), linetype = "dashed")


