nazwy <- dir("data/", pattern = "dataset")
nazwy <- nazwy[grep(nazwy, pattern = ".Rdata")]
nazwy <- gsub(nazwy, pattern = ".Rdata", replacement = "")

datasets <- strsplit(nazwy, "_")
maper <- unique(sapply(datasets, "[[", 2))
datasets <- unique(sapply(datasets, "[[", 1))

method <- c("edgeR", "DESeq", "limma1", "limma2")
library(dplyr)
library(purrr)
library(tidyr)

library(tibble)
library(reshape)

lista.dziwnych <- list()
lista.niedziwnych <- list()

for (dataset in datasets)
{
  load(file = paste0("results/difficult_genes_", dataset, ".RData"))
  lista.dziwnych[[dataset]] <- list()
  lista.niedziwnych[[dataset]] <- list()
  for (met in method)
  {
    lista.dziwnych[[dataset]][[met]] <-
      res.dane[[met]]$stats %>%
      rownames_to_column("nazwy") %>%
      filter(status.mappers == "Sig" & status.group == "Sig") %>%
      arrange(pval.map, pval.gr)
    
    lista.niedziwnych[[dataset]][[met]] <-
      res.dane[[met]]$stats %>%
      rownames_to_column("nazwy") %>%
      filter(status.mappers != "Sig" & status.group == "Sig") %>%
      arrange(pval.gr, -pval.map)
  }  
}

source("codes/differential_analysis_functions.R")
for (dataset in datasets)
{
  klas <- unlist(read.table(file = paste0('data/', dataset, '.csv')))
  for (map in maper)
  {
    load(file = paste0("data/", dataset, "_", map, ".Rdata"))
    
    for (met in method)
    {
      
      res.dane1 <- obliczenia.all.methods1(cc, klas, method = met)
      res.dane1 <- res.dane1 %>% arrange(padj)
      roznicujace <- res.dane1$ID[which(res.dane1$padj < 0.05)]
      
      lista.dziw1[[dataset]][[map]][[met]] <- intersect(lista.dziwnych[[dataset]][[met]]$nazwy, roznicujace)
      
      lista.niedziw1[[dataset]][[map]][[met]] <- setdiff(roznicujace, lista.dziwnych[[dataset]][[met]]$nazwy) 
      
    }
  }  
}
save(lista.dziw1, lista.niedziw1, file = "results/lists_of_diificult_and_not_difficult_genes_mappers_separately.RData")
