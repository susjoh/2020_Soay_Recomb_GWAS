library(biomaRt)
library(dplyr)


ensembl_hsapiens   = useEnsembl(biomart="ensembl", dataset=c("hsapiens_gene_ensembl"))
ensembl_oaries     = useEnsembl(biomart="ensembl", dataset=c("oaries_gene_ensembl"))
ensembl_btaurus    = useEnsembl(biomart="ensembl", dataset=c("btaurus_gene_ensembl"))
ensembl_mmusculus  = useEnsembl(biomart="ensembl", dataset=c("mmusculus_gene_ensembl"))

genedesc <- function(gene_names){
  
  temptab <- list()
  
  gene_names <- na.omit(gene_names)
  
  try(temptab[[1]] <- cbind(Species = "hsapiens",
                            getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_hsapiens)))
  try(temptab[[2]] <- cbind(Species = "oaries",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_oaries)))
  try(temptab[[3]] <- cbind(Species = "btaurus",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_btaurus)))
  try(temptab[[4]] <- cbind(Species = "mmusculus",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_mmusculus)))
  
  temptab <- bind_rows(temptab)
  
}

orthodesc <- function(gene_ids){
  
  gene.orthos <- NULL
  
  # mice
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_mmusculus)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "mmmusculus"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_mmusculus))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  # humans
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_hsapiens)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "hsapiens"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_hsapiens))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  # cattle
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_btaurus)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "btaurus"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_btaurus))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  #~~ Format the table:
  
  x <- gene.orthos[,c("sheep_gene_id", "external_gene_name")]
  x$external_gene_name <- toupper(x$external_gene_name)
  x <- unique(x)
  x <- subset(x, external_gene_name != "")
  x
  
}