#
# 3. Gene and GO information around the top hits
# SEJ
#

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
tophits <- subset(gwas.results, Pc1df < 1.28e-6 & Chromosome != 0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract genes around the top hits                #
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

temp <- subset(tophits, select = c(SNP.Name, Pc1df))
topregions_withSNP <- join(topregions, temp)


genedist <- 2e6

topgenes <- list()

for(i in 1:nrow(topregions)){
  
  print(paste("Running row", i, "of", nrow(topregions)))
  
  x <- getSheepGenes(chrid = topregions$Chromosome[i],
                     chrstart = topregions$MinPos[i] - genedist,
                     chrstop = topregions$MaxPos[i] + genedist)
  
  if(nrow(x) > 0){
    topgenes[[i]] <- x
  }
  
  rm(x)
}

beepr::beep()

topgenes <- bind_rows(topgenes)
toporthos <- getSheepOrthologues(unique(topgenes$ensembl_gene_id))

write.table(toporthos, "results/7_Top_Orthologous_Regions_v3.txt", row.names = F, sep = "\t", quote = F)
write.table(topgenes, "results/7_Top_Genes_v3.txt", row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract genes around orthologous hits            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

toporthos <- arrange(toporthos, Species, Chromosome.scaffold.name, Gene.start..bp.)

toporthos$Diff <- c(diff(toporthos$Gene.start..bp.), NA)
toporthos$Group <- NA

toporthos$Group [1] <- 1

for(i in 2:nrow(toporthos)){
  toporthos$Group[i] <- ifelse(toporthos$Species[i] == toporthos$Species[i-1] &
                                 toporthos$Chromosome.scaffold.name[i] == toporthos$Chromosome.scaffold.name[i-1] &
                                 toporthos$Diff[i-1] < 1e6,
                               toporthos$Group[i-1], toporthos$Group[i-1] + 1)
  
}


toporthoregions <- toporthos %>% 
  group_by(Species, Group) %>% 
  summarise(Chr = Chromosome.scaffold.name[1],
            MinPos = min(Gene.start..bp.) - genedist,
            MaxPos = max(Gene.end..bp.) + genedist) %>%
  na.omit %>%
  filter(MinPos != 0)

toporthoregions$Biomart <- paste0("ensembl_", toporthoregions$Species)

toportholist <- list()

for(i in 1:nrow(toporthoregions)){
  
  print(paste("Running row", i, "of", nrow(toporthoregions)))
  
  chrid = toporthoregions$Chr[i]
  chrstart = toporthoregions$MinPos[i]
  chrstop = toporthoregions$MaxPos[i]
  
  x <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                            'chromosome_name', 'start_position', 'end_position'), 
             filters = c('chromosome_name','start','end'),
             values = list(chrid, chrstart, chrstop), 
             mart = eval(parse(text = toporthoregions$Biomart[i])), quote = "")
  
  head(x)
  
  
  
  
  if(nrow(x) > 0){
    
    x$chromosome_name <- as.character(x$chromosome_name)
    x$Species <- toporthoregions$Species[i]
    
    toportholist[[i]] <- x
  }
  
  rm(x, chrid, chrstart, chrstop)
}

beepr::beep()

toportholist <- bind_rows(toportholist)

rm(toporthoregions, firstRun, i)

#~~ Integrate the lists...

head(toportholist)
head(toporthos)

names(toportholist) <- c("Gene.stable.ID.1", "Gene.name.1", 
                         "Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Species")
toportholist <- subset(toportholist, !Gene.stable.ID.1 %in% toporthos$Gene.stable.ID.1)
toporthos$Diff <- NULL
toporthos$Group <- NULL

toportholist$Gene.stable.ID <- NA
toportholist$Gene.name <- NA

toportholist$X.id..target.Sheep.gene.identical.to.query.gene <- NA

toporthos <- rbind(toporthos, toportholist)
toporthos <- arrange(toporthos, Species, Chromosome.scaffold.name, Gene.start..bp.)

rm(toportholist)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Get all gene and ortho GO terms            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

allgenes <- data.frame(Species = "oaries", Gene.stable.ID.1 = topgenes$ensembl_gene_id)
allgenes <- rbind(allgenes, toporthos[,c("Species", "Gene.stable.ID.1")])
allgenes <- unique(allgenes)
allgenes$Species <- as.character(allgenes$Species)
allgenes$Gene.stable.ID.1 <- as.character(allgenes$Gene.stable.ID.1)

topGO <- NULL

for(i in as.character(unique(allgenes$Species))){
  print(paste("Running", i))
  x <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description',
                            'phenotype_description', 'go_id', 'name_1006', 'definition_1006',
                            'chromosome_name', 'start_position', 'end_position', 'strand'), 
             filters = c('ensembl_gene_id'),
             values = list(subset(allgenes, Species == i)$Gene.stable.ID.1), 
             mart = eval(parse(text = paste0("ensembl_", i))), quote = "")
  x$Species <- i
  topGO <- rbind(topGO, x)
  rm(x)
}

rm(allgenes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Integrate into a virtual gene list         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# toportho_doublecheck <- subset(toporthos, !is.na(Gene.stable.ID))
# 
# toportho_doublecheck <- toportho_doublecheck %>% add_count(Gene.stable.ID, Species) %>% subset(n > 1)
# 
# toporthos$Species %>% table
# 
# topregions
# 
# for(i in 1:nrow(topregions)){
#   
#   x <- subset(topgenes, chromosome_name == topregions$Chromosome[i] &
#                 end_position >= topregions$MinPos[i]-genedist &
#                 start_position <= topregions$MaxPos[i]+genedist)
#   
#   xortho <- subset(toporthos, Gene.stable.ID %in% x$ensembl_gene_id)
#   
#   
#   for(j in unique(xo$Species)){
#     
#     xorthorange <- subset(xortho, Species == j)
#     x
#   }
#   
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Get all gene and ortho GO terms            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

meioterms <- c("meio", "recombin")


head(topGO)

meiolines <- NULL

for(i in meioterms){
  
  print(i)
  
  print(length(c(grep(i, topGO$external_gene_name),
          grep(i, topGO$phenotype_description),
          grep(i, topGO$name_1006),
          grep(i, topGO$definition_1006))))
  
  meiolines <- c(meiolines,
                 grep(i, topGO$external_gene_name),
                 grep(i, topGO$phenotype_description),
                 grep(i, topGO$name_1006),
                 grep(i, topGO$definition_1006))
}

meiolines <- unique(meiolines)

topGOmeio <- topGO[meiolines,]
table(topGOmeio$external_gene_name, topGOmeio$Species)

write.table(topGO, "results/7_Top_GO_Meio_v3.txt", row.names = F, sep = "\t", quote = F)
write.table(topGOmeio, "results/7_Top_GO_Terms_Meio_v3.txt", row.names = F, sep = "\t", quote = F)



finallist <- toporthos %>% group_by(Gene.stable.ID, Gene.name.1, Species) %>% summarise(Count = n()) %>% data.frame
str(topgenes)
str(finallist)

names(finallist)[1] <- "ensembl_gene_id"

finallist <- join(finallist, topgenes)
finallist <- finallist[,c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name", "Species", "Gene.name.1")]

finallist <- subset(finallist, !is.na(ensembl_gene_id))
finallist <- arrange(finallist, chromosome_name, start_position)










