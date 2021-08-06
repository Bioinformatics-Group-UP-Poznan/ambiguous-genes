wyn_dziwne1 <- list()
wyn_niedziwne1 <- list()
nazwy <- dir("data/", pattern = "dataset")
nazwy <- nazwy[grep(nazwy, pattern = ".Rdata")]
nazwy <- gsub(nazwy, pattern = ".Rdata", replacement = "")

datasets <- strsplit(nazwy, "_")
maper <- unique(sapply(datasets, "[[", 2))
datasets <- unique(sapply(datasets, "[[", 1))

method <- c("edgeR", "DESeq", "limma1", "limma2")


source("codes/klasAll.R")
source("codes/ensemble.R")

# before this code you have to run codes from file "lists_of_difficult_and_not_difficult_genes_mappers_separately.R"
load(file = "results/lists_of_diificult_and_not_difficult_genes_mappers_separately.RData")

for (dataset in datasets)
{
  
  L1 <- unlist(read.table(file = paste0('data/', dataset, '.csv')))
  K <- unique(L1)
  K <- as.character(K)
  
  L1 <- as.character(L1)
  for (i in 1:length(K))
  {
    a <- which(L1 == K[i])
    L1[a] <- i-1
  }
  L1<-as.numeric(L1)
  L<-as.factor(L1)
  
  for (met in method[1])
  {
    for (map in maper)
    {
      load(file = paste0("data/", dataset, "_", map, ".Rdata"))
      y <- DGEList(cc, group = L)
      y <- calcNormFactors(y)
      cc <- cpm(y, normalized.lib.sizes=TRUE)
      set.seed(08022018)
      
      ktore <- match(lista.dziw1[[dataset]][[map]][[met]], rownames(cc))
      if (length(ktore) > 5)
      {
        wyniki_dziwne <- NULL
        for (i in 1:10)
        {
          for(f in c(1/3,1/2, 2/3))
          {
            data<-as.matrix(t(cc[ktore[1:min(floor(length(L)*f),length(ktore))],]))
            misclassified_dziwne <- klasAll1(data, L)
            wyniki_dziwne <- rbind(wyniki_dziwne, cbind(misclassified_dziwne,f = min(floor(length(L)*f),length(ktore)), sym = i))
          }
          cat(i, "dziwne ")
        }
        wyn_dziwne1[[dataset]][[map]][[met]] <- wyniki_dziwne
      } else {
        wyn_dziwne1[[dataset]][[map]][[met]] <- data.frame(true.y = L, yhat = NA, prob = NA, f= min(floor(length(L)*f),length(ktore)), sym = 1)
      }  
      
      set.seed(08022018)
      ktore <- match(lista.niedziw1[[dataset]][[map]][[met]], rownames(cc))
      if (length(ktore) > 5)
      {
        wyniki_niedziwne<-NULL
        for(i in 1:10)
        {
          for(f in c(1/3,1/2, 2/3))
          {
            data<-as.matrix(t(cc[ktore[1:min(floor(length(L)*f),length(ktore))],]))
            misclassified_niedziwne <- klasAll1(data, L)
            wyniki_niedziwne <- rbind(wyniki_niedziwne,cbind(misclassified_niedziwne, f = min(floor(length(L)*f),length(ktore)), sym = i))
          }
          cat(i, "niedziwne ")
        }
        wyn_niedziwne1[[dataset]][[map]][[met]] <- wyniki_niedziwne
      } else {
        wyn_niedziwne1[[dataset]][[map]][[met]] <- data.frame(true.y = L, yhat = NA, prob = NA, f= min(floor(length(L)*f),length(ktore)), sym = 1)
      }  
      cat(map)
    } 
    cat(met)
  }
  cat(dataset)
}

save(wyn_dziwne1, wyn_niedziwne1, file = paste0("results/joint_classifier_edgeR_cpm.RData"))