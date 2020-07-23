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

load("data/GRM.RData")

source("r/ASReml4.EstEffects.R")
source("r/makeGRM.R")


attr(grminv, which = "INVERSE") <- TRUE

pedigree <- read.table("data/20190711_Soay_Pedigree.txt", header = T, stringsAsFactors = F)

for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])

ainv <- ainverse(pedigree)

recsumm <- read.table("results/3_20200721_Individual_Recombination_Rates.txt", header = T, stringsAsFactors = F)

recsumm$RRID <- factor(recsumm$RRID)
recsumm$RRID.SEX <- factor(recsumm$RRID.SEX)

recsumm <- subset(recsumm, !is.na(RRID.SEX) & !is.na(RRID.Fhat3))

recsumm.m <- subset(recsumm, RRID.SEX == "M") %>% droplevels
recsumm.f <- subset(recsumm, RRID.SEX == "F") %>% droplevels


# Top GWAS Hits

gwas.results <- read.table("results/6_GWAS_Results_600K.txt", header = T, stringsAsFactors = F, sep = "\t")
gwas.results.to.keep <- subset(gwas.results, SNP.Name %in% c("oar3_OAR7_21368818", "oar3_OAR6_116402578"))

gwas.results$Sex[which(gwas.results$Sex == "Both Sexes")] <- "A. Both Sexes"
gwas.results$Sex[which(gwas.results$Sex == "Females")] <- "B. Females"
gwas.results$Sex[which(gwas.results$Sex == "Males")] <- "C. Males"
gwas.results <- subset(gwas.results, Pc1df < 1.28e-6)

load("soayimp_genotype_data.RData")

#~~ Get the top SNPS

tophits <- arrange(gwas.results, Chromosome, Position)
tophits <- subset(tophits, Chromosome != 0)

tophits$Group <- NA
tophits$Group[1] <- 1

for(i in 2:nrow(tophits)) tophits$Group[i] <- ifelse(tophits$Chromosome[i] == tophits$Chromosome[i-1] &
                                                       tophits$Position[i] - tophits$Position[i-1] < 2e6 &
                                                       tophits$Position[i] - tophits$Position[i-1] >= 0, 
                                                     tophits$Group[i-1], tophits$Group[i-1] + 1)

topregions <- tophits %>% group_by(Group, Sex) %>% 
  summarise(Chromosome = mean(Chromosome),
            MinPos = min(Position),
            MaxPos = max(Position),
            Pc1df = min(Pc1df)) %>%
  na.omit %>%
  filter(MinPos != 0)

topsnps <- NULL

for(i in 1:nrow(topregions)){
  
  topsnps <- rbind(topsnps,
                   subset(tophits, Sex == topregions$Sex[i] & Chromosome == topregions$Chromosome[i] & Pc1df == topregions$Pc1df[i]))
  
}


topsnpsvec <- c(unique(topsnps$SNP.Name))

#~~ Get the genotypes and add to recsumm

temp <- as.character.gwaa.data(soayimp[,topsnpsvec]) %>% data.frame
temp$RRID <- row.names(temp)
head(temp)

recsumm <- join(recsumm, temp)
head(recsumm)

rm(soayimp, temp, i)
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run Models to get genotype effect sizes      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

res.fixef <- NULL
res.wald <- NULL

for(i in topsnpsvec){
  
  x.data <- recsumm[,c("TotalRecombCount", "RRID", "RRID.SEX", "RRID.Fhat3", i)] %>% na.omit
  names(x.data)[5] <- "SNP"
  x.data$RRID2 <- x.data$RRID
  
  x.both <- asreml(fixed = TotalRecombCount ~ RRID.SEX + SNP + RRID.Fhat3,
                   random = ~ vm(RRID, ainv) + ide(RRID),
                   data = x.data,
                   residual = ~idv(units))
  
  x.male <- asreml(fixed = TotalRecombCount ~ SNP + RRID.Fhat3,
                   random = ~ vm(RRID, ainv) + ide(RRID),
                   data = droplevels(subset(x.data, RRID.SEX == "M")),
                   residual = ~idv(units))
  
  
  x.female <- asreml(fixed = TotalRecombCount ~ SNP + RRID.Fhat3,
                     random = ~ vm(RRID, ainv) + ide(RRID),
                     data = droplevels(subset(x.data, RRID.SEX == "F")),
                     residual = ~idv(units))
  
  x.both.1 <- summary.asreml(x.both, coef = T)$coef.fixed %>% data.frame
  x.male.1 <- summary.asreml(x.male, coef = T)$coef.fixed %>% data.frame
  x.female.1 <- summary.asreml(x.female, coef = T)$coef.fixed %>% data.frame
  
  x.both.1$Effect <- row.names(x.both.1)
  x.male.1$Effect <- row.names(x.male.1)
  x.female.1$Effect <- row.names(x.female.1)
  
  x.both.1$Sex <- "Both Sexes"
  x.male.1$Sex <- "Male"
  x.female.1$Sex <- "Female"
  
  x.both.1$SNP <- i
  x.male.1$SNP <- i
  x.female.1$SNP <- i
  
  res.fixef <- rbind(res.fixef, x.both.1, x.male.1, x.female.1)
  
  x.both.2 <- wald.asreml(x.both) %>% data.frame
  x.male.2 <- wald.asreml(x.male) %>% data.frame
  x.female.2 <- wald.asreml(x.female) %>% data.frame
  
  x.both.2$Effect <- row.names(x.both.2)
  x.male.2$Effect <- row.names(x.male.2)
  x.female.2$Effect <- row.names(x.female.2)
  
  x.both.2$Sex <- "Both Sexes"
  x.male.2$Sex <- "Male"
  x.female.2$Sex <- "Female"
  
  x.both.2$SNP <- i
  x.male.2$SNP <- i
  x.female.2$SNP <- i
  
  res.wald <- rbind(res.wald, x.both.2, x.male.2, x.female.2)
  
  rm(x.both.2, x.male.2, x.female.2, x.both.1, x.male.1, x.female.1, x.both, x.male, x.female, x.data)
  gc()
  
}

res.fixef.hold <- res.fixef

res.fixef <- subset(res.fixef, Sex != "Both Sexes")


#~~ Plot the effect sizes from the models.

x.snp <- res.fixef[grep("SNP", res.fixef$Effect),]
x.snp$Effect <- gsub("SNP_", "", x.snp$Effect)

x.int <- res.fixef[grep("Intercept", res.fixef$Effect),]
x.int

x.plot <- NULL

for(i in 1:nrow(x.int)){
  
  x1 <- subset(x.snp, Sex == x.int$Sex[i] & SNP == x.int$SNP[i])
  x1$solution <- x1$solution + x.int$solution[i]
  x.plot <- rbind(x.plot, x1)
  rm(x1)
  
}

rm(x.snp, x.int)

x.plot <- subset(x.plot, Sex != "Both Sexes")
x.plot

ggplot(x.plot, aes(Effect, solution)) +
  geom_point() +
  geom_errorbar(aes(ymin = solution - std.error, ymax = solution + std.error), width = 0.1) +
  facet_grid(SNP ~ Sex)

x.plot <- subset(x.plot, SNP %in% c("oar3_OAR7_21368818", "oar3_OAR6_116402578"))

x.plot$Gene <- ifelse(x.plot$SNP == "oar3_OAR6_116402578", "RNF212 (oar3_OAR6_116402578)", "RNF212B (oar3_OAR7_21368818)")


head(recsumm)

newx <- subset(recsumm, select = c(RRID, RRID.SEX, TotalRecombCount, oar3_OAR6_116402578, oar3_OAR7_21368818))
newx <- melt(newx, id.vars = c("RRID", "RRID.SEX", "TotalRecombCount"))
head(newx)
newx$RRID.SEX <- ifelse(newx$RRID.SEX == "F", "Female", "Male")
newx$Gene <- ifelse(newx$variable == "oar3_OAR6_116402578", "RNF212 (oar3_OAR6_116402578)", "RNF212B (oar3_OAR7_21368818)")
newx <- na.omit(newx)
newx2 <- subset(newx, select = -TotalRecombCount) %>% unique


samplesizes <- left_join(
  newx %>% group_by(RRID.SEX, variable, value) %>% summarize(count = n()),
  
  newx2 %>% group_by(RRID.SEX, variable, value) %>% summarize(idcount = n())
)

samplesizes$idcount <- paste0("(", samplesizes$idcount, ")")
names(samplesizes)[1:3] <- c("Sex", "Locus", "Effect")


samplesizes <- melt(samplesizes, id.vars = c("Sex", "Locus", "Effect"))


samplesizes$Gene <- ifelse(samplesizes$Locus == "oar3_OAR6_116402578", "RNF212 (oar3_OAR6_116402578)", "RNF212B (oar3_OAR7_21368818)")
samplesizes$y <- ifelse(samplesizes$variable == "count", 42.3, 41)

ggplot() +
  geom_point(data = x.plot, aes(Effect, solution), size = 2) +
  geom_errorbar(data = x.plot, aes(x = Effect, y = solution, ymin = solution - std.error, ymax = solution + std.error), width = 0.1) +
  geom_rect(data = x.plot, mapping=aes(ymin = 39.5, ymax = 43, xmin = 0, xmax = 4), fill = "white", col = "black") +
  geom_text(data = samplesizes, aes(Effect, y, label = value)) +
  facet_grid(Gene ~ Sex) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, colour = "black"),
        axis.text.y  = element_text (size = 12, colour = "black"),
        strip.text = element_text (size = 12),
        axis.title.y = element_text (size = 12, vjust = 2),
        axis.title.x = element_text (size = 12),
        legend.position = "none") +
  labs(x = "SNP Genotype", y = "Autosomal Crossover Count") +
  scale_y_continuous(breaks = seq(20, 50, 2.5)) +
  coord_cartesian(ylim = c(22, 43), expand = F)
ggsave("figs/8_SNP_Effect_Sizes_Fixef.png", height = 7, width = 7)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Run Models to get regional Va                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

res.regh2 <- NULL

for(i in c("oar3_OAR7_21368818", "oar3_OAR6_116402578")){
  
  x.data <- recsumm[,c("TotalRecombCount", "RRID", "RRID.SEX", "RRID.Fhat3", i)] %>% na.omit
  names(x.data)[5] <- "SNP"
  x.data$RRID2 <- x.data$RRID
  
  #~~ fit the GRMs
  
  grm.region <- read.table(paste0("animal_models/snplist__", i, "_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
  ids.region <- read.table(paste0("animal_models/snplist__", i, "_GRM.adj.grm.id"))  # CONTAINS ID LIST

  
  # crop grm to only IDs I need
  
  grmregion <- makeGRM(grm.region, ids.region, id.vector = recsumm$RRID)
  attr(grmregion, which = "INVERSE") <- TRUE
  
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
  
  
  rm(x.both.reg, x.male.reg, x.female.reg, grmregion)
  gc()
  
  
}

res.regh2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Make Table for Paper                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.snp <- res.fixef.hold[grep("SNP", res.fixef.hold$Effect),]
x.snp$Effect <- gsub("SNP_", "", x.snp$Effect)

x.snp <- subset(x.snp, SNP %in% c("oar3_OAR7_21368818", "oar3_OAR6_116402578"))


gwas.results.to.keep <- arrange(gwas.results.to.keep, Chromosome, Position, Sex)
gwas.results.to.keep <- subset(gwas.results.to.keep, select = -c(A1, A2, effB, se_effB, chi2.1df, P1df, ExpP, Cumu))
names(gwas.results.to.keep) <- c("Sex", "SNP", "Chromosome", "Position", "P", "MAF")

x <- subset(x.snp, select = c(Sex, SNP, Effect, solution, std.error))
x <- subset(x, Effect != "A/A")

x1 <- dcast(x, Sex + SNP ~ Effect, value.var = "solution")
x2 <- dcast(x, Sex + SNP ~ Effect, value.var = "std.error")
names(x2)[3:4] <- c("G/A_se", "G/G_se")

x <- join(x1, x2)
rm(x1, x2)
x$Sex[which(x$Sex == "Male")] <- "Males"
x$Sex[which(x$Sex == "Female")] <- "Females"

gwas.results.to.keep <- join(gwas.results.to.keep, x)
write.table(gwas.results.to.keep, "results/8_Top_Hits.txt", row.names = F, sep = "\t", quote = F)
