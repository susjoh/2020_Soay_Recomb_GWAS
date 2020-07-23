#
# Figures for Paper
# SEJ July 2020
#
#

library(GenABEL)
library(magrittr)
library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)

source("r/multiplot.R")

#~~ Autosomal Crossover Counts

recsumm <- read.table("results/3_20200721_Individual_Recombination_Rates.txt", header = T, stringsAsFactors = F)

recsumm$RRID <- factor(recsumm$RRID)
recsumm$RRID.SEX <- factor(recsumm$RRID.SEX)

recsumm <- subset(recsumm, !is.na(RRID.SEX) & !is.na(RRID.Fhat3))

recsumm.m <- subset(recsumm, RRID.SEX == "M") %>% droplevels
recsumm.f <- subset(recsumm, RRID.SEX == "F") %>% droplevels

#~~ Linkage Map

newcm <- read.table("results/4_20200504_Full_Linkage_Map_Full_Details.txt", header = T, stringsAsFactors = F)
head(newcm)

oldcm <- read.table("data/2_Linkage_Map_Positions_g.txt", header = T, stringsAsFactors = F)
head(oldcm)
oldcm <- subset(oldcm, select = c(SNP.Name, cMPosition, cMPosition.Female, cMPosition.Male))
names(oldcm)[2:4] <- paste0(names(oldcm)[2:4], ".2016")

newcm <- join(newcm, oldcm)
rm(oldcm)

head(newcm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Basic ACC Statistics                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm$Sex2 <- ifelse(recsumm$RRID.SEX == "F", "Female", "Male")

p1 <- ggplot(recsumm, aes(Sex2, TotalRecombCount, fill = Sex2)) +
  geom_boxplot(notch = 1, width = 0.5) +
  theme_bw() +
  labs(x = "SEX", y = "Autosomal Crossover Count") +
  ggtitle("A") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12, colour = "black"),
        axis.text.y  = element_text (size = 12, colour = "black"),
        strip.text.x = element_text (size = 14, hjust = 0),
        axis.title.y = element_text (size = 12, vjust = 2),
        axis.title.x = element_text (size = 12),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")

p2 <- ggplot(recsumm, aes(RRID.Fhat3, TotalRecombCount, col = Sex2)) +
  geom_point(alpha = 0.2) + 
  stat_smooth() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Genomic Inbreeding Coefficient (Fhat3)", y = "Autosomal Crossover Count", col = "Sex") +
  ggtitle("B") +
  theme(axis.text.x  = element_text (size = 12, colour = "black"),
        axis.text.y  = element_text (size = 12, colour = "black"),
        strip.text.x = element_text (size = 14, hjust = 0),
        axis.title.y = element_text (size = 12, vjust = 2),
        axis.title.x = element_text (size = 12),
        strip.background = element_rect(fill = "white"), 
        legend.position = "none")

p3 <- ggplot(recsumm, aes(RRID.Age.2, TotalRecombCount, col = Sex2)) +
  geom_point(alpha = 0.2) + 
  stat_smooth() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Age (years)", y = "Autosomal Crossover Count", col = "Sex") +
  ggtitle("C") +
  theme(axis.text.x  = element_text (size = 12, colour = "black"),
        axis.text.y  = element_text (size = 12, colour = "black"),
        strip.text.x = element_text (size = 14, hjust = 0),
        axis.title.y = element_text (size = 12, vjust = 2),
        axis.title.x = element_text (size = 12),
        strip.background = element_rect(fill = "white"), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 16, 2))

png("figs/1_Sex_Fhat3_Age_corrs.png", width = 8, height = 8, units = "in", res = 300)

multiplot(p1, p3, p2, cols = 2)

dev.off()

# grid.arrange(p1, p2, p3, p4,
#              widths = c(3, 3, 4, 2),
#              layout_matrix = matrix(c(1, 2, 3, 4), nrow = 2, byrow = T))

rm(p1, p2, p3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Linkage Maps                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

newcm_sex <- subset(newcm, select = c(SNP.Name, Chr, Oar3_Pos, cMPosition.Female, cMPosition.Male))
newcm_sex <- melt(newcm_sex, id.vars = c("SNP.Name", "Chr", "Oar3_Pos"))
head(newcm_sex)
newcm_sex$variable <- gsub("cMPosition.", "", newcm_sex$variable)

ggplot(newcm_sex, aes(Oar3_Pos/1e6, value, col = variable)) +
  geom_line(size = 1) +
  facet_wrap(~Chr, scales = "free", ncol = 5) +
  labs(x = "Oar_3.1 Position (Mb)", y = "cM Position", col = "Sex") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 10, colour = "black"),
        axis.text.y  = element_text (size = 10, colour = "black"),
        strip.text.x = element_text (size = 12, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 2),
        axis.title.x = element_text (size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave("figs/2_Linkage_Map.png", width = 8, height = 10)

#~~ Compare the lengths
linkcomp <- newcm %>%
  group_by(Chr = as.factor(Chr)) %>%
  summarise(cMPosition = max(cMPosition),
            cMPosition.Female = max(cMPosition.Female),
            cMPosition.Male = max(cMPosition.Male),
            cMPosition.2016 = max(cMPosition.2016),
            cMPosition.Female.2016 = max(cMPosition.Female.2016),
            cMPosition.Male.2016 = max(cMPosition.Male.2016))

head(linkcomp)

summary(lm(cMPosition ~ cMPosition.2016, data = linkcomp))
summary(lm(cMPosition.Male ~ cMPosition.Male.2016, data = linkcomp))
summary(lm(cMPosition.Female ~ cMPosition.Female.2016, data = linkcomp))

linkcomp <- linkcomp  %>% melt(id.vars = "Chr")

linkcomp <- separate(linkcomp, variable, into = c("variable", "Sex", "Year"), sep = "\\.", remove = T)
linkcomp$Year[which(linkcomp$Sex == 2016)] <- 2016
linkcomp$Sex[which(linkcomp$Sex == 2016)] <- "Both Sexes"
linkcomp$Sex[is.na(linkcomp$Sex)] <- "Both Sexes"
linkcomp$Year[is.na(linkcomp$Year)] <- 2020
linkcomp$Year <- paste0("Run_", linkcomp$Year)
linkcomp$variable <- NULL
linkcomp <- dcast(linkcomp, Chr + Sex ~ Year)

ggplot(linkcomp, aes(Run_2016, Run_2020, col = Sex)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  stat_smooth(method = "lm") +
  theme_bw() +
  scale_colour_manual(values = c("black", "red", "blue")) +
  theme(axis.text.x  = element_text (size = 12, colour = "black"),
        axis.text.y  = element_text (size = 12, colour = "black"),
        strip.text.x = element_text (size = 14, hjust = 0),
        axis.title.y = element_text (size = 12, vjust = 2),
        axis.title.x = element_text (size = 12)) +
  labs(x = "Johnston et al (2016) cM Length", y = "Current Study cM Length")
  scale_x_continuous(breaks = seq(0, 400, 50)) +
  scale_y_continuous(breaks = seq(0, 400, 50))
ggsave("figs/2_Linkage_Map_Correspondence.png", width = 5, height = 4)

#~~ how many more separated by crossovers?

head(newcm)

table(diff(newcm$cMPosition) == 0)
table(diff(newcm$cMPosition.2016) == 0)

write.table(newcm, "doc/tables/Table_S2_Linkage_Map.txt", row.names = F, sep = "\t", quote = F)

rm(linkcomp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. GWAS                                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gwas.results <- read.table("results/6_GWAS_Results_600K.txt", header = T, stringsAsFactors = F, sep = "\t")

gwas.results <- arrange(gwas.results, Cumu)

chrinfo <- NULL

for(j in unique(gwas.results$Chromosome)){
  temp1 <- subset(gwas.results, Chromosome == j)
  temp2 <- data.frame(Chromosome = j,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  chrinfo <- rbind(chrinfo, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

chrinfo$Chr2 <- c(0:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", "", 24, "", "", "X")


gwas.results.hold <- gwas.results

gwas.results.hold$Sex[which(gwas.results.hold$Sex == "Both Sexes")] <- "A. Both Sexes"
gwas.results.hold$Sex[which(gwas.results.hold$Sex == "Females")] <- "B. Females"
gwas.results.hold$Sex[which(gwas.results.hold$Sex == "Males")] <- "C. Males"

gwas.results <- gwas.results.hold[seq(1, nrow(gwas.results.hold), 1000),]

gwas.results <- gwas.results.hold

png("figs/2_SEJ_GWAS_600K.png", width = 6, height = 12, units = "in", res = 300)

ggplot(gwas.results, aes(Cumu,-log10(Pc1df), col = factor(Chromosome %% 2))) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr2) +
  scale_colour_brewer(palette = "Set1") +
  labs(x ="Chromosome", y = expression("-log"[10]*"P")) +
  facet_wrap(~Sex, ncol = 1)

dev.off()

png("figs/2_SEJ_GWAS_600K_PP.png", width = 4, height = 12, units = "in", res = 300)

ggplot(gwas.results, aes(-log10(ExpP),-log10(Pc1df))) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        # strip.text.x = element_text (size = 16, hjust = 0),
        # strip.text.y = element_text (size = 16),
        strip.text = element_text(colour = "white", size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x = expression("Expected -log"[10]*"P"), y = expression("Observed -log"[10]*"P")) +
  facet_wrap(~Sex, ncol = 1) +
  geom_abline(slope = 1, intercept = 0)

dev.off()


gwas.results <- subset(gwas.results, Pc1df < 1.28e-6)
gwas.results <- arrange(gwas.results, Sex, Chromosome, Position) %>%
  subset(Chromosome != 0)


