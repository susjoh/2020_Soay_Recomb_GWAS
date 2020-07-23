#
# How much Va is attributed to associated regions?
# SEJ, AMS
# Jan 2019
#
# 


library(GenABEL)
library(magrittr)
library(plyr)
library(ggplot2)
library(asreml)

source("r/ASReml.EstEffects.R")
source("r/rGLSadj.R")
source("r/pin.R")
source("r/makeGRM.R")

load("eddie_gwas/Recomb_GWAS.RData")

load("data/GRM.RData")

focalsnps <- read.table("results/4_Top_SNPs_for_Paper.txt", header = T, sep = "\t", stringsAsFactors = F)
focalsnpshold <- focalsnps
focalsnps <- subset(focalsnps, Model != "Both Sexes")
focalsnps <- subset(focalsnps, select = c(SNP.Name, Chromosome, Position)) %>% unique
# focalsnps <- rbind(cbind(Model = "Females", focalsnps),
#                    cbind(Model = "Males", focalsnps))
focalsnps$Trait <- "TotalRecombCount"

maptab <- data.frame(SNP.Name = snpnames(soayimp),
                     Chromosome = chromosome(soayimp),
                     Position = map(soayimp), stringsAsFactors = F)

maptab <- arrange(maptab, Chromosome, Position)

runModels <- TRUE

recsumm$RRID2 <- recsumm$RRID

#~~ How much variance does the region explain? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regvar <- NULL

for(i in 1:nrow(focalsnps)){
  
  snplist.name <- paste0("animal_models\\snplist__", focalsnps$SNP.Name[i])
  
  if(runModels == TRUE){

    models <- data.frame(Model = c("1+RRID.Fhat3,random=~giv(RRID)+ide(RRID)+giv(RRID2)"),
                         Sex = c("Females", "Males"),
                         Response = focalsnps$Trait[i], stringsAsFactors = F)

    focsnp <- focalsnps$SNP.Name[i]
    
    x2 <- subset(maptab, Chromosome == focalsnps$Chromosome[i])
    
    if(!focsnp %in% maptab$SNP.Name){
      x1 <- strsplit(focsnp, split = "_")[[1]]
      focsnp <- x2$SNP.Name[which(x2$Position > as.numeric(x1[3]))[1]]
    } 
    
    snplist <- x2$SNP.Name[(which(x2$SNP.Name == focsnp) - 25):(which(x2$SNP.Name == focsnp) + 25)]
    
    #~~ Deal with SNPs at ends of chromosomes
    
    if(which(x2$SNP.Name == focsnp) > nrow(x2) - 10) snplist <- x2$SNP.Name[(nrow(x2)-20):nrow(x2)]
    if(which(x2$SNP.Name == focsnp) < 11) snplist <- x2$SNP.Name[1:21]
    
    print(snplist)
    
    writeLines(snplist, con = paste0(snplist.name, ".txt"))
    
    system(paste0("gcta64.exe --bfile data\\sheep_geno_imputed_ram_27092019 --autosome --autosome-num 26 --extract ", snplist.name, ".txt --make-grm-gz --out ", snplist.name, "_GRM" ))
    system(paste0("gcta64.exe --grm-gz ", snplist.name, "_GRM --grm-adj 0 --make-grm-gz --out ", snplist.name, "_GRM.adj"))
    
    grm.region <- read.table(paste0(snplist.name, "_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
    ids.region <- read.table(paste0(snplist.name, "_GRM.adj.grm.id"))  # CONTAINS ID LIST
    
    # crop grm to only IDs I need
    
    grmregion <- makeGRM(grm.region, ids.region, id.vector = recsumm$RRID)
    
    #~~ Make data for model
    
    x.vec <- as.character(models$Model[1])
    x.vec <- gsub(",random=~", "+", x.vec)
    x.vec <- strsplit(x.vec, split = "\\+")[[1]]
    x.vec[grep("giv", x.vec)] <- "RRID"
    x.vec[grep("ide", x.vec)] <- "RRID"
    x.vec <- c(x.vec, as.character(focalsnps$Trait[i]), "RRID2")
    x.vec <- unique(x.vec[-which(x.vec == 1)])
    
    x.data <- subset(recsumm, RRID.SEX == "F")
    
    x.data <- droplevels(na.omit(x.data[,x.vec]))
    
    x.data <- droplevels(subset(x.data, RRID %in% dimnames(grminv)[[1]]))
    x.data <- droplevels(subset(x.data, RRID %in% dimnames(grmregion)[[1]]))
    
    
    try({
      eval(parse(text = paste0("fit1 <- asreml(fixed=", focalsnps$Trait[i], "~",
                               models$Model[1],"
                             , data=x.data, ginverse=list(RRID=grminv, RRID2 = grmregion),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)")))
      
      save(fit1, file = paste0(snplist.name, "_Females.RData"))

    })
    
    beepr::beep()
    
    x <- cbind(Trait = focalsnps$Trait[i],
               Sex = "Females",
               Chromsome = focalsnps$Chromosome[i],
               SNP.Name = focalsnps$SNP.Name[i],
               print(ASReml.EstEffects(fit1)))
    regvar <- rbind(regvar, x)
    
    
    x.data <- subset(recsumm, RRID.SEX == "M")
    
    x.data <- droplevels(na.omit(x.data[,x.vec]))
    
    x.data <- droplevels(subset(x.data, RRID %in% dimnames(grminv)[[1]]))
    x.data <- droplevels(subset(x.data, RRID %in% dimnames(grmregion)[[1]]))
    
    
    try({
      eval(parse(text = paste0("fit1 <- asreml(fixed=", focalsnps$Trait[i], "~",
                               models$Model[1],"
                             , data=x.data, ginverse=list(RRID=grminv, RRID2 = grmregion),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)")))
      
      save(fit1, file = paste0(snplist.name, "_Males.RData"))
      
    })
    
    rm(x.data, x.vec, grm.region, grmregion, ids.region, snplist)
    
  } else {
    
    load(file = paste0(snplist.name, ".RData"))
    
  }
  
  x <- cbind(Trait = focalsnps$Trait[i],
             Sex = "Males",
             Chromsome = focalsnps$Chromosome[i],
             SNP.Name = focalsnps$SNP.Name[i],
             print(ASReml.EstEffects(fit1)))
  
  
  regvar <- rbind(regvar, x)
  
  rm(snplist.name, fit1, x)
}





regvar$Model <- paste0(regvar$Trait, "__", regvar$Sex, "__", regvar$SNP.Name)


recodetab <- data.frame(variable = c("giv(RRID).giv",
                                     "giv(RRID2).giv",
                                     "ide(RRID)!id",
                                     "R!variance"),
                        new.variable = c("Additive Genetic",
                                         "Regional Genetic",
                                         "Permanent Environment",
                                         "Residual"))

regvar <- join(regvar, recodetab)

regvar$new.variable2 <- factor(regvar$new.variable,
                               levels = rev(c("Additive Genetic", "Regional Genetic", "Permanent Environment", "Residual")))

ggplot(regvar, aes(SNP.Name, Effect, fill = new.variable2)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(fill = "Variable") +
  facet_wrap(~Sex) +
  theme(axis.text.x  = element_text (size = 12, angle = 270, vjust = 0.2),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14))


x <- reshape2::dcast(subset(regvar, select = c(Model, new.variable, Effect)), Model ~ new.variable)
x$Vp <- rowSums(x[,2:ncol(x)], na.rm = T)
x$GenVar <- x$`Regional Genetic`/(x$`Regional Genetic` + x$`Additive Genetic`)
x$GenVar <- round(x$GenVar, digits = 3)
x$`Regional Genetic` <- round(x$`Regional Genetic`, digits = 4)

subset(x, select = c(Model,  `Regional Genetic`, GenVar))

# Females RNF212B explains 16% of genetic variance
# Females RNF212 explains 32% of genetic variance
# Females ZWINT explains 6% of genetic variance

# Males RNF212B explains 25% of genetic variance


write.table(regvar, "results/5_Variance_Attributed_by_associated_Regions.txt", row.names = F, sep = "\t", quote = F)

