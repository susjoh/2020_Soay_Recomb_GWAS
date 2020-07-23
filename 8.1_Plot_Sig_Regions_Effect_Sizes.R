
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load in and format data                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(plyr)
library(magrittr)
library(tidyr)
library(biomaRt)
library(reshape2)

source("r/biomaRt_functions.R")

#~~ Read in data

# Top GWAS Hits

gwas.results <- read.table("results/6_GWAS_Results_600K.txt", header = T, stringsAsFactors = F, sep = "\t")
gwas.results$Sex[which(gwas.results$Sex == "Both Sexes")] <- "A. Both Sexes"
gwas.results$Sex[which(gwas.results$Sex == "Females")] <- "B. Females"
gwas.results$Sex[which(gwas.results$Sex == "Males")] <- "C. Males"

tophits <- subset(gwas.results, Pc1df < 1.28e-6 & Chromosome != 0)

load("data/Plates_1to87_QC3.GenABEL.RData")

gwas.results$SNP50 <- ifelse(gwas.results$SNP.Name %in% snpnames(soay87), "SNP50", "SNPHD")

meiogenes <- read.table("results/7_Top_GO_Terms_Meio_v3.txt", header = T, stringsAsFactors = F, sep = "\t", quote = "")
geneorthos <- read.table("results/7_Top_Orthologous_Regions_v3.txt", header = T, stringsAsFactors = F, sep = "\t", quote = "")
meiorthos <- subset(geneorthos, Gene.stable.ID.1 %in% meiogenes$ensembl_gene_id)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract regions around the top hits              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Set a distance from the SNP to pick up genes

tophits <- arrange(tophits, Chromosome, Position)

tophits$Group <- NA
tophits$Group[1] <- 1

for(i in 2:nrow(tophits)) tophits$Group[i] <- ifelse(tophits$Chromosome[i] == tophits$Chromosome[i-1] &
                                                       tophits$Position[i] - tophits$Position[i-1] < 2e6 &
                                                       tophits$Position[i] - tophits$Position[i-1] >= 0, 
                                                     tophits$Group[i-1], tophits$Group[i-1] + 1)

topregions <- tophits %>% group_by(Group) %>% 
  summarise(Chromosome = mean(Chromosome),
            MinPos = min(Position),
            MaxPos = max(Position),
            Pc1df = min(Pc1df)) %>%
  na.omit %>%
  filter(MinPos != 0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Plot Chromosome 6                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

topregions

chr6 <- subset(gwas.results, Chromosome == 6 & Position > 110000000 & Position < 118000000)

#~~ Get gene information

x6 <- getSheepGenes(chrid = 6,
                    chrstart = 110000000,
                    chrstop = 118000000)
x6 <- subset(x6, gene_biotype == "protein_coding")

x6$Meiotic <- ifelse(x6$ensembl_gene_id %in% meiorthos$Gene.stable.ID, "yes", "no")

x6 <- rbind(x6, data.frame(ensembl_gene_id = "RNF212",
                           external_gene_name = "RNF212", 
                           chromosome_name = 6,
                           start_position = 116427802,
                           end_position = 116446904,
                           gene_biotype = "protein_coding",
                           Meiotic = "yes"))
x6 <- arrange(x6, start_position)
x6$Orthos <- NA

for(i in 1:nrow(x6)){
  temp <- subset(geneorthos, Gene.stable.ID == x6$ensembl_gene_id[i] & Gene.name.1 != "")
  if(nrow(temp) > 1) x6$Orthos[i] <- temp$Gene.name.1 %>% toupper %>% unique %>% paste(collapse = ", ")
  rm(temp)
}

subset(x6, Meiotic == "yes")

head(x6)
x6$y <- rep(19:22, length.out = nrow(x6))
x6plot <- subset(x6, select = c(ensembl_gene_id, Meiotic, y, start_position, end_position))
x6plot <-  melt(x6plot, id.vars = c("ensembl_gene_id", "y", "Meiotic"))
x6plot <- subset(x6plot, Meiotic == "yes")


ggplot() +
  geom_rect(data = x6plot, mapping=aes(ymin = 19, ymax = 21, xmin = 109, xmax = 118), fill = "white", col = "black") +
  geom_point(data = chr6, aes(Position/1e6,-log10(Pc1df), col = SNP50), size = 1.5, alpha = 0.6) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  geom_line(data = x6plot, aes(x = value/1e6, y = 20, group = ensembl_gene_id), size = 3) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Chromosome 6 Position (Mb)", y = expression("-log"[10]*"P"), col = "SNP Chip") +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(breaks = seq(0, 150, 1)) +
  facet_wrap(~Sex, ncol = 1) +
  coord_cartesian(xlim = c(110, 117.2), ylim = c(-0.3, 19.9))
ggsave("figs/8_Chr_6_Region.png", height = 8, width = 8)


ggplot() +
  #geom_vline(xintercept = 116.427802, linetype = "dashed", colour = "lightblue") +
  geom_point(data = subset(chr6, Sex == "B. Females"), aes(Position/1e6,-log10(Pc1df), col = SNP50), size = 1.5, alpha = 0.6) +
  geom_line(data = x6plot, aes(x = value/1e6, y = 19, group = ensembl_gene_id), size = 3, colour = "blue") +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Chromosome 6 Position (Mb)", y = expression("-log"[10]*"P"), col = "SNP Chip") +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(breaks = seq(0, 150, 1))


ggsave("figs/8_Chr_6_Region_Females.png", height = 3.5, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Plot Chromosome 7                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

topregions

chr7 <- subset(gwas.results, Chromosome == 7 & Position > 19000000 & Position < 26000000)


#~~ Get gene information

x7 <- getSheepGenes(chrid = 7,
                    chrstart = 19000000,
                    chrstop = 26000000)
x7 <- subset(x7, gene_biotype == "protein_coding")

x7$Meiotic <- ifelse(x7$ensembl_gene_id %in% meiorthos$Gene.stable.ID, "yes", "no")
x7 <- arrange(x7, start_position)
x7$Orthos <- NA

for(i in 1:nrow(x7)){
  temp <- subset(geneorthos, Gene.stable.ID == x7$ensembl_gene_id[i] & Gene.name.1 != "")
  if(nrow(temp) > 1) x7$Orthos[i] <- temp$Gene.name.1 %>% toupper %>% unique %>% paste(collapse = ", ")
  rm(temp)
}

subset(x7, Meiotic == "yes")

ggplot(chr7, aes(Position/1e6,-log10(Pc1df), col = SNP50)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Chromosome 7 Position (Mb)", y = expression("-log"[10]*"P"), col = "SNP Chip") +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(breaks = seq(0, 150, 1)) +
  facet_wrap(~Sex, ncol = 1)
ggsave("figs/8_Chr_7_Region.png", height = 8, width = 8)

head(x7)
x7$y <- rep(14:17, length.out = nrow(x7))
x7plot <- subset(x7, select = c(ensembl_gene_id, Meiotic, y, start_position, end_position))
x7plot <-  melt(x7plot, id.vars = c("ensembl_gene_id", "y", "Meiotic"))
x7plot <- subset(x7plot, Meiotic == "yes")

subset(x7, Meiotic == "yes")

ggplot() +
  geom_rect(data = x7plot, mapping=aes(ymin = 13, ymax = 15, xmin = 18, xmax = 27), fill = "white", col = "black") +
  geom_point(data = chr7, aes(Position/1e6,-log10(Pc1df), col = SNP50), size = 1.5, alpha = 0.6) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  geom_line(data = x7plot, aes(x = value/1e6, y = 14, group = ensembl_gene_id), size = 3) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Chromosome 7 Position (Mb)", y = expression("-log"[10]*"P"), col = "SNP Chip") +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(breaks = seq(0, 150, 1)) +
  facet_wrap(~Sex, ncol = 1) +
  coord_cartesian(xlim = c(19, 26), ylim = c(-0.3, 14.2))
ggsave("figs/8_Chr_7_Region.png", height = 8, width = 8)





ggplot() +
  # geom_vline(xintercept = 20.601666 , linetype = "dashed", colour = "lightblue") +
  geom_vline(xintercept = 21.241831, linetype = "dashed", colour = "lightblue") +
  geom_point(data = subset(chr7, Sex == "B. Females"), aes(Position/1e6,-log10(Pc1df), col = SNP50), size = 1.5, alpha = 0.6) +
  #geom_line(data = x7plot, aes(x = value/1e6, y = 12, group = ensembl_gene_id)) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Chromosome 7 Position (Mb)", y = expression("-log"[10]*"P"), col = "SNP Chip") +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(breaks = seq(0, 150, 1))

ggsave("figs/8_Chr_7_Region_BothSexes.png", height = 3.5, width = 8)


genelist <- rbind(x6, x7)
genelist$y <- NULL
head(genelist)

write.table(genelist, "doc/tables/Gene_List_Table.txt", row.names = F, sep = "\t", quote = F)




