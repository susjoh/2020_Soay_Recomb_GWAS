#
#
# Animal models for meiosis traits
#
#
#

library(GenABEL)
library(magrittr)
library(plyr)
library(ggplot2)
library(asreml)
library(reshape2)
library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load data and set up working environment     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ LOad recombination info
recsumm <- read.table("results/3_20200721_Individual_Recombination_Rates.txt", header = T, stringsAsFactors = F)

recsumm$RRID <- factor(recsumm$RRID)
recsumm$RRID.SEX <- factor(recsumm$RRID.SEX)

recsumm <- subset(recsumm, !is.na(RRID.SEX) & !is.na(RRID.Fhat3))

recsumm.m <- subset(recsumm, RRID.SEX == "M") %>% droplevels
recsumm.f <- subset(recsumm, RRID.SEX == "F") %>% droplevels

#~~ Load asreml related stuff

load("data/GRM.RData")

source("r/ASReml4.EstEffects.R")
source("r/makeGRM.R")

attr(grminv, which = "INVERSE") <- TRUE


#~~ Get the genotypes and add to recsumm

load("soayimp_genotype_data.RData")

topsnpsvec <- c("oar3_OAR7_21368818", "oar3_OAR6_116402578")

temp <- as.character.gwaa.data(soayimp[,topsnpsvec]) %>% data.frame
temp$RRID <- row.names(temp)
head(temp)

recsumm <- join(recsumm, temp)
head(recsumm)

rm(soayimp, temp, i)
gc()

#~~ Get 50K data to extract SNPs

mapdata <- read.table("../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Previous_Versions/Plates_1to87_QC3.bim")
head(mapdata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Make GRMs                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get 20K SNP lists

x7 <- subset(mapdata, V1 == 7)
which(x7$V4 > 21368818)[1]
which(x7$V4 > 21116299)[1]
snplist7 <- x7$V2[361:380] %>% as.character

x6 <- subset(mapdata, V1 == 6)
nrow(x6)
snplist6 <- x6$V2[1837:1856] %>% as.character

#~~ Make relatedness matrices

# setwd("../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Previous_Versions/")
# 
# writeLines(snplist6, "snplist_6.txt")
# system("gcta64.exe --bfile Plates_1to87_QC3 --autosome --autosome-num 26 --extract snplist_6.txt --make-grm-gz --out snplist_6_GRM")
# system("gcta64.exe --grm-gz snplist_6_GRM --grm-adj 0 --make-grm-gz --out snplist_6_GRM.adj")
# 
# writeLines(snplist7, "snplist_7.txt")
# system("gcta64.exe --bfile Plates_1to87_QC3 --autosome --autosome-num 26 --extract snplist_7.txt --make-grm-gz --out snplist_7_GRM")
# system("gcta64.exe --grm-gz snplist_7_GRM --grm-adj 0 --make-grm-gz --out snplist_7_GRM.adj")
# 
# setwd("../../../Recombination Projects/20200721_Soay_Recomb_GWAS/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Run Models to get regional Va                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.data <- recsumm[,c("TotalRecombCount", "RRID", "RRID.SEX", "RRID.Fhat3")] %>% na.omit
x.data$RRID2 <- x.data$RRID
x.data$RRID3 <- x.data$RRID


res.regh2 <- NULL

for(i in topsnpsvec){
  

  #~~ fit the GRMs
  
  if(i == "oar3_OAR6_116402578"){
    
    grm.region <- read.table(paste0("animal_models/snplist_6_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
    ids.region <- read.table(paste0("animal_models/snplist_6_GRM.adj.grm.id"))  # CONTAINS ID LIST
    
  }
  
  if(i == "oar3_OAR7_21368818"){
    
    grm.region <- read.table(paste0("animal_models/snplist_7_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
    ids.region <- read.table(paste0("animal_models/snplist_7_GRM.adj.grm.id"))  # CONTAINS ID LIST
    
  }
  
  # crop grm to only IDs I need
  
  grmregion <- makeGRM(grm.region, ids.region, id.vector = recsumm$RRID)
  attr(grmregion, which = "INVERSE") <- TRUE
  
  rm(grm.region, ids.region)
  
  x.both.reg <- NULL
  x.male.reg <- NULL
  x.female.reg <- NULL
  
  try(x.both.reg <- asreml(fixed = TotalRecombCount ~ RRID.SEX ,
                           random = ~ vm(RRID, grminv) + vm(RRID2, grmregion) ,
                           data = x.data,
                           residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))
  
  try(res.regh2 <- rbind(res.regh2,
                         cbind(asreml4pin(x.both.reg), Sex = "Both Sexes", SNP = i)))
  
  try(x.male.reg <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                           random = ~ vm(RRID, grminv) + vm(RRID2, grmregion) ,
                           data = droplevels(subset(x.data, RRID.SEX == "M")),
                           residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))
  try(res.regh2 <- rbind(res.regh2,
                         cbind(asreml4pin(x.male.reg), Sex = "Males", SNP = i)))
  
  try(x.female.reg <- asreml(fixed = TotalRecombCount ~  1,
                             random = ~ vm(RRID, grminv) + vm(RRID2, grmregion) ,
                             data = droplevels(subset(x.data, RRID.SEX == "F")),
                             residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))
  
  try(res.regh2 <- rbind(res.regh2,
                         cbind(asreml4pin(x.female.reg), Sex = "Females", SNP = i)))
  
  beepr::beep()
  rm(x.both.reg, x.male.reg, x.female.reg, grmregion)
  gc()
  
  
}

res.regh2

write.table(res.regh2, "results/8_Regional_Heritability_Single_Window.txt", row.names = F, sep = "\t", quote = F)

grm.region6 <- read.table(paste0("animal_models/snplist_6_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.region6 <- read.table(paste0("animal_models/snplist_6_GRM.adj.grm.id"))  # CONTAINS ID LIST

grm.region7 <- read.table(paste0("animal_models/snplist_7_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.region7 <- read.table(paste0("animal_models/snplist_7_GRM.adj.grm.id"))  # CONTAINS ID LIST

grm.region6 <- makeGRM(grm.region6, ids.region6, id.vector = recsumm$RRID)
attr(grm.region6, which = "INVERSE") <- TRUE

grm.region7 <- makeGRM(grm.region7, ids.region7, id.vector = recsumm$RRID)
attr(grm.region7, which = "INVERSE") <- TRUE




try(x.both.reg <- asreml(fixed = TotalRecombCount ~ RRID.SEX ,
                         random = ~ vm(RRID, grminv) + vm(RRID2, grm.region6) + vm(RRID3, grm.region7) ,
                         data = x.data,
                         residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))

try(res.regh2 <- rbind(res.regh2,
                       cbind(asreml4pin(x.both.reg), Sex = "Both Sexes", SNP = i)))

try(x.male.reg <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                         random = ~ vm(RRID, grminv) + vm(RRID2, grm.region6) + vm(RRID3, grm.region7) ,
                         data = droplevels(subset(x.data, RRID.SEX == "M")),
                         residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))
try(res.regh2 <- rbind(res.regh2,
                       cbind(asreml4pin(x.male.reg), Sex = "Males", SNP = i)))

try(x.female.reg <- asreml(fixed = TotalRecombCount ~  1,
                           random = ~ vm(RRID, grminv) + vm(RRID2, grm.region6) + vm(RRID3, grm.region7) ,
                           data = droplevels(subset(x.data, RRID.SEX == "F")),
                           residual = ~idv(units), workspace = 500e+6, pworkspace = 500e+6))

try(res.regh2 <- rbind(res.regh2,
                       cbind(asreml4pin(x.female.reg), Sex = "Females", SNP = i)))


x.both.reg$vparameters

asreml4pin(x.both.reg)
vpredict(x.both.reg, h1 ~ V1/(V1+V2+V3))
vpredict(x.both.reg, h1 ~ V2/(V1+V2+V3))

asreml4pin(x.male.reg)
vpredict(x.male.reg, h1 ~ V1/(V1+V2+V3))
vpredict(x.male.reg, h1 ~ V2/(V1+V2+V3))

asreml4pin(x.female.reg)
vpredict(x.female.reg, h1 ~ V1/(V1+V2+V3))
vpredict(x.female.reg, h1 ~ V2/(V1+V2+V3))

beepr::beep()



write.table(res.regh2, "results/8_Regional_Heritability_FULL.txt", row.names = F, sep = "\t", quote = F)



