#
# Compile Meiosis Data from Johnston et al Genetics (2016)
# and updated dataset from 2019 runs.
# SEJ May 2019
#

library(plyr)
library(ggplot2)
library(tidyr)
library(GenABEL)
library(dplyr)

auto.vec <- 1:26
sex.vec <- 27

AnalysisSuffix <- "b"

load("data/Plates_1to87_QC3.GenABEL.RData")

abeldata <- soay87

soayibc <- read.table("data/Plates_1to87_QC3.ibc", header = T, stringsAsFactors = F)

basedata <- read.table("data/20191002_SoayBaseData.txt", header = T, stringsAsFactors = F, sep = "\t")

recsumm <- read.table("results/2_Per_Individual_Recomb_Final_b.txt", header = T, stringsAsFactors = F)


#~~ add ibc info

soayibc <- subset(soayibc, select = c("IID", "Fhat3"))
names(soayibc)[1] <- "ID"

basedata <- join(basedata, soayibc)
head(basedata)

#~~ Compile recsumm

recsumm <- subset(recsumm, select = c(RRID, ANIMAL, TotalInfLoci, TotalRecombCount, NonRecombCount))
names(recsumm)[2] <- "Offspring.ID"

basedata <- subset(basedata, select = c(ID, BirthYear, Sex, Fhat3))
basedata$Sex <- ifelse(basedata$Sex == 1, "F", ifelse(basedata$Sex == 2, "M", NA))

names(basedata) <- c("RRID", "RRID.BYEAR", "RRID.SEX", "RRID.Fhat3")
recsumm <- join(recsumm, basedata)

names(basedata) <- c("Offspring.ID", "Offspring.BYEAR", "Offspring.SEX", "Offspring.Fhat3")
recsumm <- join(recsumm, basedata)

recsumm$RRID.Age <- recsumm$Offspring.BYEAR - recsumm$RRID.BYEAR

ggplot(recsumm, aes(TotalInfLoci, TotalRecombCount)) + geom_point()

recsumm <- subset(recsumm, TotalInfLoci > 6500)
recsumm$RRID.Age.2 <- recsumm$RRID.Age
recsumm$RRID.Age.2[which(recsumm$RRID.SEX == "M")] <- recsumm$RRID.Age.2[which(recsumm$RRID.SEX == "M")] + 0.2

recsumm$test <- paste(recsumm$RRID, "_", recsumm$Offspring.ID)
table(table(recsumm$test))

recsumm$RRID.SEX[which(recsumm$RRID == 11807)] <- "F"

table(recsumm$RRID.SEX, useNA = "always")



ggplot(recsumm, aes(RRID.Age.2, TotalRecombCount, col = RRID.SEX)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth() +
  labs(x = "Age", y = "Autosomal Crossover Count", col = "Sex") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.text = element_text (size = 14),
        legend.title = element_text (size = 14),
        legend.position = "top")

ggsave("figs/2_Age_Association.png", width = 6, height = 5)

recsumm$test <- NULL

#~~~~~ Compile rectab

rectab <- read.table("results/2_Per_Chromosome_Recomb_Final_b.txt", header = T, stringsAsFactors = F)

table(rectab$RecombCount)

#~~ Compile a table of the crossover positions.

maptab <- read.table("data/2_Linkage_Map_Positions_g.txt", header = T, stringsAsFactors = F)
maptab <- subset(maptab, Chr != sex.vec)

maptemp <- data.frame(SNP.Name = snpnames(abeldata),
                      Position = map(abeldata))
maptab <- join(maptab, maptemp)
maptab <- subset(maptab, !is.na(Position))

sum(nchar(subset(rectab, Family == "Offspring_10000_FATHER_7492")$data)) == nrow(maptab)

recpos <- list()
head(rectab)
rectab$Offspring.ID <- rectab$ANIMAL


for(i in 1:nrow(rectab)){
  
  if(i %in% seq(1, nrow(rectab), 1000)) print(paste("Running row", i, "of", nrow(rectab)))
  
  if(rectab$RecombCount[i] > 0){
    x <- rectab$data[i]
    x1 <- data.frame(Base = strsplit(x, split = "")[[1]],
                     Position = subset(maptab, Chr == rectab$Chr[i])$Position)
    x1 <- subset(x1, Base != "-")
    x1 <- cbind(x1[-nrow(x1),], x1[-1,])
    names(x1) <- c("Base1", "Pos1", "Base2", "Pos2")
    x1 <- subset(x1, Base1 != Base2)
    x1$RRID <- rectab$RRID[i]
    x1$Offspring.ID <- rectab$Offspring.ID[i]
    x1$Chr <- rectab$Chr[i]
    x1$Order <- 1:nrow(x1)
    
    recpos[[i]] <- x1
    
    rm(x, x1)
  }
  
}

recpos <- bind_rows(recpos)

head(recpos)

recpos$Distance <- recpos$Pos2 - recpos$Pos1

maxchr <- data.frame(MaxChr = tapply(maptab$Position, maptab$Chr, max))
maxchr$Chr <- 1:nrow(maxchr)

recpos <- join(recpos, maxchr)

recpos$EstPosCO <- recpos$Pos1 + (0.5*recpos$Distance)

#~~ Get statistics 

# Distance to first crossover from centromere

x <- subset(recpos, Order == 1 & Chr %in% 4:26)
head(x)
x$DistanceFirstCO <- x$EstPosCO
x <- subset(x, select = c(RRID, Offspring.ID, Chr, DistanceFirstCO))

rectab <- join(rectab, x)

#~~ Mean distance between double crossovers

recpos$Index <- paste(recpos$RRID, recpos$Offspring.ID, recpos$Chr, sep = "_")
head(recpos)
recpos$Difference <- c(diff(recpos$EstPosCO), 0)
recpos <- recpos[which(recpos$Index[1:(nrow(recpos)-1)] == recpos$Index[2:(nrow(recpos))]),]
x <- data.frame(DistanceDoubleCO = tapply(recpos$Difference, recpos$Index, mean))
x$Index <- row.names(x)

x <- separate(x, Index, c("RRID", "Offspring.ID", "Chr"), "_")
head(x)

rectab <- join(rectab, x)

head(rectab)
rectab$Index <- paste0(rectab$RRID, "_", rectab$Offspring.ID)

x <- data.frame(MeanDistanceFirstCO = tapply(rectab$DistanceFirstCO, rectab$Index, mean, na.rm = T),
                MeanDistanceDoubleCO = tapply(rectab$DistanceDoubleCO, rectab$Index, mean, na.rm = T))
x$Index <- row.names(x)

head(x)
x <- separate(x, Index, c("RRID", "Offspring.ID"), "_")

recsumm <- join(recsumm, x)
head(recsumm)

write.table(recsumm, "results/3_20200721_Individual_Recombination_Rates.txt", row.names = F, sep = "\t", quote = F)

write.table(rectab, "results/3_20200721_Individual_Recombination_Rates_per_Chr.txt", row.names = F, sep = "\t", quote = F)
