wyszukiwarka <- function(county1, county2, uklad, p.val=0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  if (any(dim(county1)!=dim(county2))) {
    stop("Tables of counts have to have the same sizes.")
  }
  n <- ncol(county1)
  county <- cbind(county1, county2)
  colnames(county) <- c(paste0("mapper1-",1:n), paste0("mapper2-",1:n))
  isexpr <- rowSums(cpm(county) > 5) >= 2
  county <- county[isexpr,]
  grupy <- as.factor(c(uklad, uklad)) 
  mapery <- as.factor(c(rep("mapper1",n),rep("mapper2",n)))
  
  switch(method,
  edgeR = 
  {      
    design1 <- model.matrix( ~ grupy + mapery)
    y <- DGEList(county)
    y <- calcNormFactors(y)
    y1 <- estimateDisp(y, design1)
    fit <- glmFit(y1, design1)
    lrt1 <- glmLRT(fit, coef=2:(ncol(design1)-1))
    lrt2 <- glmLRT(fit, coef = ncol(design1))
    lib.sizes <- y$samples[,3]
    
    out <- data.frame(counts = county)
    out$adj.pval <- p.adjust(lrt2$table$PValue, method = "BH")
    out$status.group <- ifelse(p.adjust(lrt1$table$PValue, method = "BH")<p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(out$adj.pval <p.val, "Sig", "nonSig")
  },
  limma = 
  {
    design1 <- model.matrix( ~ grupy + mapery)
    lib.sizes <- calcNormFactors(county, method = "TMM")
    voom.data = voom(county, design = design1,
                     lib.size = colSums(county) * lib.sizes)
    voom.data$genes = rownames(county)
    voom.fitlimma1 = lmFit(voom.data, design = design1)
    voom.fitbayes1 = eBayes(voom.fitlimma1)
    res1<-topTable(voom.fitbayes1, coef = 2:(ncol(design1)-1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
    res2<-topTable(voom.fitbayes1, coef = ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
    out <- data.frame(counts = county)
    out$adj.pval <- res2$adj.P.Val
    out$status.group <- ifelse(res1$adj.P.Val < p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(res2$adj.P.Val < p.val, "Sig", "nonSig")
  },
  DESeq = 
  {
    design1 <- data.frame(grupy = grupy, mapery = mapery)
    cdsFull = newCountDataSet(round(county), design1)
    cdsFull = estimateSizeFactors( cdsFull )
    cdsFull = estimateDispersions( cdsFull )
    fit0 = fitNbinomGLMs( cdsFull, count ~  grupy + mapery )
    fit1 = fitNbinomGLMs( cdsFull, count ~  grupy )
    fit2 = fitNbinomGLMs( cdsFull, count ~  mapery )
    lib.sizes <- sizeFactors(cdsFull)
    
    pvalsGLM1 = nbinomGLMTest(fit0, fit1)
    padjGLM1 = p.adjust( pvalsGLM1, method="BH" )
    pvalsGLM2 = nbinomGLMTest(fit0, fit2)
    padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
    out <- data.frame(counts = county)
    out$adj.pval <- padjGLM1
    out$status.group <- ifelse(padjGLM2 < p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(padjGLM1 < p.val, "Sig", "nonSig")
  })
  
  res <- list()
  res[["stats"]] <- out
  res[["lib.sizes"]]  <- lib.sizes  
  res
}


wyszukiwarka5 <- function(county, n.mapper = 3, design, p.val=0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  if (ncol(county) %% n.mapper != 0) {
    stop("Number of columns is not a multiplication of number of mappers.")
  }
  n <- ncol(county)/n.mapper
  colnames(county) <- paste0("mapper",rep(1:n.mapper, each=n),"-",rep(1:n,times=n.mapper))
  isexpr <- rowSums(county > 5) >= n/3
  county <- county[isexpr,]
  grupy <- as.factor(rep(design, n.mapper)) 
  mapery <- as.factor(rep(paste0("mapper",1:n.mapper), each = n))
  n.groups <- length(unique(grupy))
  
  switch(method,
         edgeR = 
         {      
           design1 <- model.matrix( ~ grupy + mapery)
           y <- DGEList(county)
           y <- calcNormFactors(y)
           y1 <- estimateDisp(y, design1)
           fit <- glmFit(y1, design1)
           lrt1 <- glmLRT(fit, coef=2:n.groups)
           lrt2 <- glmLRT(fit, coef =(n.groups+1):ncol(design1))
           lib.sizes <- y$samples[,3]
           
           b1 <- match(rownames(lrt1$table),rownames(county))
           b2 <- match(rownames(lrt2$table),rownames(county))
           out <- data.frame(counts = county)
           out$pval.gr <- p.adjust(lrt1$table$PValue[b1], method = "BH")
           out$pval.map <- p.adjust(lrt2$table$PValue[b2], method = "BH")
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         limma1 = 
         {
           design1 <- model.matrix( ~ grupy + mapery)
           lib.sizes <- calcNormFactors(county, method = "TMM")
           voom.data = voom(county, design = design1,
                            lib.size = colSums(county) * lib.sizes)
           voom.data$genes = rownames(county)
           voom.fitlimma1 = lmFit(voom.data, design = design1)
           voom.fitbayes1 = eBayes(voom.fitlimma1)
           res1<-topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           res2<-topTable(voom.fitbayes1, coef = (n.groups+1):ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           out <- data.frame(counts = county)
           b1 <- match(rownames(res1),rownames(county))
           b2 <- match(rownames(res2),rownames(county))
           out$pval.gr <- res1$adj.P.Val[b1]
           out$pval.map <- res2$adj.P.Val[b2]
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         limma2 = 
         {
           design1 <- model.matrix( ~ grupy + mapery)
           DESeq.cds1 <- newCountDataSet(countData = round(county), conditions = design1)
           DESeq.cds1 <- estimateSizeFactors(DESeq.cds1)
           lib.sizes <- sizeFactors(DESeq.cds1)
           DESeq.cds1 <- estimateDispersions(DESeq.cds1, method = "blind", fitType = "local")
           DESeq.vst1 <- DESeq::getVarianceStabilizedData(DESeq.cds1)
           voom.fitlimma1 = lmFit(DESeq.vst1, design = design1)
           voom.fitbayes1 = eBayes(voom.fitlimma1)
           res1<-topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           res2<-topTable(voom.fitbayes1, coef = (n.groups+1):ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           out <- data.frame(counts = county)
           b1 <- match(rownames(res1),rownames(county))
           b2 <- match(rownames(res2),rownames(county))
           out$pval.gr <- res1$adj.P.Val[b1]
           out$pval.map <- res2$adj.P.Val[b2]
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         DESeq = 
         {
           design1 <- data.frame(grupy = grupy, mapery = mapery)
           cdsFull = newCountDataSet(round(county), design1)
           cdsFull = estimateSizeFactors( cdsFull )
           cdsFull = estimateDispersions( cdsFull )
           fit0 = fitNbinomGLMs( cdsFull, count ~  grupy + mapery )
           fit1 = fitNbinomGLMs( cdsFull, count ~  grupy )
           fit2 = fitNbinomGLMs( cdsFull, count ~  mapery )
           lib.sizes <- sizeFactors(cdsFull)
           
           pvalsGLM1 = nbinomGLMTest(fit0, fit1)
           padjGLM1 = p.adjust( pvalsGLM1, method="BH" )
           pvalsGLM2 = nbinomGLMTest(fit0, fit2)
           padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
           out <- data.frame(counts = county)
           b1 <- match(rownames(fit1),rownames(county))
           b2 <- match(rownames(fit2),rownames(county))
           out$pval.gr <- padjGLM2[b2]
           out$pval.map <- padjGLM1[b1]
           out$status.group <- ifelse(padjGLM2[b2] < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(padjGLM1[b1] < p.val, "Sig", "nonSig")
         })
  
  res <- list()
  res[["stats"]] <- out
  res[["lib.sizes"]]  <- lib.sizes  
  res
}


vis.diff <- function(county1, county2, gene.list)
{
  require(ggplot2, reshape)
  
  if (length(gene.list) > 16) {
    stop("Maximum number of genes to plot is 16.")
  }
  if (any(dim(county1)!=dim(county2))) {
    stop("Tables of counts have to have the same sizes.")
  }
  
  n <- ncol(county1)
  dane <- list()
  dane[["mapper1"]] <- county1[match(gene.list, rownames(county1)),]
  dane[["mapper2"]] <- county2[match(gene.list, rownames(county2)),]
  colnames(dane[[1]]) <- 1:n
  colnames(dane[[2]]) <- 1:n
  res <- melt(dane)
  res <- cbind(res, geny = rep(rownames(county)[match(gene.list, rownames(county1))], 2*n))
  p <- ggplot(data=res, aes(x=variable, y=value, group=L1, color=L1)) +
    geom_line() + theme_bw() + facet_wrap( ~ geny, scales = "free", ncol=4) + 
    theme(legend.position = "bottom", legend.title=element_blank()) +
    ylab("Counts") + xlab("Samples")
  print(p)
}