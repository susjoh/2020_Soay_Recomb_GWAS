library(biomaRt)
library(dplyr)


ensembl_hsapiens   = useEnsembl(biomart="ensembl", dataset=c("hsapiens_gene_ensembl"), host = "www.ensembl.org")
ensembl_oaries     = useEnsembl(biomart="ensembl", dataset=c("oaries_gene_ensembl"), host = "www.ensembl.org")
ensembl_btaurus    = useEnsembl(biomart="ensembl", dataset=c("btaurus_gene_ensembl"), host = "www.ensembl.org")
ensembl_mmusculus  = useEnsembl(biomart="ensembl", dataset=c("mmusculus_gene_ensembl"), host = "www.ensembl.org")
ensembl_rnorvegicus  = useEnsembl(biomart="ensembl", dataset=c("rnorvegicus_gene_ensembl"), host = "www.ensembl.org")

getSheepGenes <- function(chrid, chrstart, chrstop){
  
  x <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                            'chromosome_name', 'start_position', 'end_position', 'gene_biotype'), 
             filters = c('chromosome_name','start','end'),
             values = list(chrid, chrstart, chrstop), 
             mart = ensembl_oaries, quote = "")
  
  x
  
}


getSheepOrthologues <- function(sheep_ensemble_gene_id){
  
  y <- rbind(cbind(Species = "hsapiens",
                   getLDS(attributes= c("ensembl_gene_id", "external_gene_name"), 
                          filters = "ensembl_gene_id",
                          values = sheep_ensemble_gene_id, 
                          mart = ensembl_oaries,
                          attributesL = c("ensembl_gene_id", "external_gene_name", 'chromosome_name', 'start_position', 'end_position', "oaries_homolog_perc_id"),
                          martL=ensembl_hsapiens)),
             
             cbind(Species = "btaurus",
                   getLDS(attributes= c("ensembl_gene_id", "external_gene_name"), 
                          filters = "ensembl_gene_id",
                          values = sheep_ensemble_gene_id, 
                          mart = ensembl_oaries,
                          attributesL = c("ensembl_gene_id", "external_gene_name", 'chromosome_name', 'start_position', 'end_position', "oaries_homolog_perc_id"),
                          martL=ensembl_btaurus)),
             
             cbind(Species = "mmusculus",
                   getLDS(attributes= c("ensembl_gene_id", "external_gene_name"), 
                          filters = "ensembl_gene_id",
                          values = sheep_ensemble_gene_id, 
                          mart = ensembl_oaries,
                          attributesL = c("ensembl_gene_id", "external_gene_name", 'chromosome_name', 'start_position', 'end_position', "oaries_homolog_perc_id"),
                          martL=ensembl_mmusculus)),
             
             cbind(Species = "rnorvegicus",
                   getLDS(attributes= c("ensembl_gene_id", "external_gene_name"), 
                          filters = "ensembl_gene_id",
                          values = sheep_ensemble_gene_id, 
                          mart = ensembl_oaries,
                          attributesL = c("ensembl_gene_id", "external_gene_name", 'chromosome_name', 'start_position', 'end_position', "oaries_homolog_perc_id"),
                          martL=ensembl_rnorvegicus))
             
  )
  
  y           
}
