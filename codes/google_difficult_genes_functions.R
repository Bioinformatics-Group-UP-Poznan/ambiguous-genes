dg.seeker <- function(counts1, counts2, exp.design, p.val = 0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  if (any(dim(counts1)!=dim(counts2))) {
    stop("Tables of counts have to have the same sizes.")
  }
  n <- ncol(counts1)
  counts <- cbind(counts1, counts2)
  colnames(counts) <- c(paste0("mapper1-",1:n), paste0("mapper2-",1:n))
  isexpr <- rowSums(cpm(counts) > 5) >= 2
  counts <- counts[isexpr,]
  groups <- as.factor(c(exp.design, exp.design)) 
  mapers <- as.factor(c(rep("mapper1",n),rep("mapper2",n)))
  
  switch(method,
  edgeR = 
  {      
    design1 <- model.matrix( ~ groups + mapers)
    y <- DGEList(counts)
    y <- calcNormFactors(y)
    y1 <- estimateDisp(y, design1)
    fit <- glmFit(y1, design1)
    lrt1 <- glmLRT(fit, coef=2:(ncol(design1)-1))
    lrt2 <- glmLRT(fit, coef = ncol(design1))
    lib.sizes <- y$samples[,3]
    
    out <- data.frame(counts)
    out$adj.pval <- p.adjust(lrt2$table$PValue, method = "BH")
    out$status.group <- ifelse(p.adjust(lrt1$table$PValue, method = "BH")<p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(out$adj.pval <p.val, "Sig", "nonSig")
  },
  limma = 
  {
    design1 <- model.matrix( ~ groups + mapers)
    lib.sizes <- calcNormFactors(counts, method = "TMM")
    voom.data = voom(counts, design = design1,
                     lib.size = colSums(counts) * lib.sizes)
    voom.data$genes = rownames(counts)
    voom.fitlimma1 = lmFit(voom.data, design = design1)
    voom.fitbayes1 = eBayes(voom.fitlimma1)
    res1<-topTable(voom.fitbayes1, coef = 2:(ncol(design1)-1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
    res2<-topTable(voom.fitbayes1, coef = ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
    out <- data.frame(counts = counts)
    out$adj.pval <- res2$adj.P.Val
    out$status.group <- ifelse(res1$adj.P.Val < p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(res2$adj.P.Val < p.val, "Sig", "nonSig")
  },
  DESeq = 
  {
    design1 <- data.frame(groups, mapers)
    cdsFull = newCountDataSet(round(counts), design1)
    cdsFull = estimateSizeFactors( cdsFull )
    cdsFull = estimateDispersions( cdsFull )
    fit0 = fitNbinomGLMs( cdsFull, count ~  groups + mapers )
    fit1 = fitNbinomGLMs( cdsFull, count ~  groups )
    fit2 = fitNbinomGLMs( cdsFull, count ~  mapers )
    lib.sizes <- sizeFactors(cdsFull)
    
    pvalsGLM1 = nbinomGLMTest(fit0, fit1)
    padjGLM1 = p.adjust( pvalsGLM1, method="BH" )
    pvalsGLM2 = nbinomGLMTest(fit0, fit2)
    padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
    out <- data.frame(counts = counts)
    out$adj.pval <- padjGLM1
    out$status.group <- ifelse(padjGLM2 < p.val, "Sig", "nonSig")
    out$status.mappers <- ifelse(padjGLM1 < p.val, "Sig", "nonSig")
  })
  
  res <- list()
  res[["stats"]] <- out
  res[["lib.sizes"]]  <- lib.sizes  
  res[["no.filtered"]] <- length(isexpr)
  res
}


dg.seeker5 <- function(counts, n.mapper = 3, exp.design, p.val = 0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  if (ncol(counts) %% n.mapper != 0) {
    stop("Number of columns is not a multiplication of number of mappers.")
  }
  n <- ncol(counts)/n.mapper
  colnames(counts) <- paste0("mapper", rep(1:n.mapper, each = n), "-", rep(1:n, times = n.mapper))
  isexpr <- rowSums(counts > 5) >= n/3
  counts <- counts[isexpr,]
  groups <- as.factor(rep(exp.design, n.mapper)) 
  mapers <- as.factor(rep(paste0("mapper", 1:n.mapper), each = n))
  n.groups <- length(unique(groups))
  
  switch(method,
         edgeR = 
         {      
           design1 <- model.matrix( ~ groups + mapers)
           y <- DGEList(counts)
           y <- calcNormFactors(y)
           y1 <- estimateDisp(y, design1)
           fit <- glmFit(y1, design1)
           lrt1 <- glmLRT(fit, coef=2:n.groups)
           lrt2 <- glmLRT(fit, coef =(n.groups+1):ncol(design1))
           lib.sizes <- y$samples[,3]
           
           b1 <- match(rownames(lrt1$table),rownames(counts))
           b2 <- match(rownames(lrt2$table),rownames(counts))
           out <- data.frame(counts = counts)
           out$pval.gr <- p.adjust(lrt1$table$PValue[b1], method = "BH")
           out$pval.map <- p.adjust(lrt2$table$PValue[b2], method = "BH")
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         limma1 = 
         {
           design1 <- model.matrix( ~ groups + mapers)
           lib.sizes <- calcNormFactors(counts, method = "TMM")
           voom.data = voom(counts, design = design1,
                            lib.size = colSums(counts) * lib.sizes)
           voom.data$genes = rownames(counts)
           voom.fitlimma1 = lmFit(voom.data, design = design1)
           voom.fitbayes1 = eBayes(voom.fitlimma1)
           res1<-topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           res2<-topTable(voom.fitbayes1, coef = (n.groups+1):ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           out <- data.frame(counts = counts)
           b1 <- match(rownames(res1),rownames(counts))
           b2 <- match(rownames(res2),rownames(counts))
           out$pval.gr <- res1$adj.P.Val[b1]
           out$pval.map <- res2$adj.P.Val[b2]
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         limma2 = 
         {
           design1 <- model.matrix( ~ groups + mapers)
           DESeq.cds1 <- newCountDataSet(countData = round(counts), conditions = design1)
           DESeq.cds1 <- estimateSizeFactors(DESeq.cds1)
           lib.sizes <- sizeFactors(DESeq.cds1)
           DESeq.cds1 <- estimateDispersions(DESeq.cds1, method = "blind", fitType = "local")
           DESeq.vst1 <- DESeq::getVarianceStabilizedData(DESeq.cds1)
           voom.fitlimma1 = lmFit(DESeq.vst1, design = design1)
           voom.fitbayes1 = eBayes(voom.fitlimma1)
           res1<-topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           res2<-topTable(voom.fitbayes1, coef = (n.groups+1):ncol(design1), number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           out <- data.frame(counts = counts)
           b1 <- match(rownames(res1),rownames(counts))
           b2 <- match(rownames(res2),rownames(counts))
           out$pval.gr <- res1$adj.P.Val[b1]
           out$pval.map <- res2$adj.P.Val[b2]
           out$status.group <- ifelse(out$pval.gr < p.val, "Sig", "nonSig")
           out$status.mappers <- ifelse(out$pval.map < p.val, "Sig", "nonSig")
         },
         DESeq = 
         {
           design1 <- data.frame(groups = groups, mapers = mapers)
           cdsFull = newCountDataSet(round(counts), design1)
           cdsFull = estimateSizeFactors( cdsFull )
           cdsFull = estimateDispersions( cdsFull )
           fit0 = fitNbinomGLMs( cdsFull, count ~  groups + mapers )
           fit1 = fitNbinomGLMs( cdsFull, count ~  groups )
           fit2 = fitNbinomGLMs( cdsFull, count ~  mapers )
           lib.sizes <- sizeFactors(cdsFull)
           
           pvalsGLM1 = nbinomGLMTest(fit0, fit1)
           padjGLM1 = p.adjust( pvalsGLM1, method="BH" )
           pvalsGLM2 = nbinomGLMTest(fit0, fit2)
           padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
           out <- data.frame(counts = counts)
           b1 <- match(rownames(fit1),rownames(counts))
           b2 <- match(rownames(fit2),rownames(counts))
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


vis.diff <- function(counts1, counts2, gene.list)
{
  require(ggplot2, reshape)
  
  if (length(gene.list) > 16) {
    stop("Maximum number of genes to plot is 16.")
  }
  if (any(dim(counts1)!=dim(counts2))) {
    stop("Tables of counts have to have the same sizes.")
  }
  
  n <- ncol(counts1)
  dane <- list()
  dane[["mapper1"]] <- counts1[match(gene.list, rownames(counts1)),]
  dane[["mapper2"]] <- counts2[match(gene.list, rownames(counts2)),]
  colnames(dane[[1]]) <- 1:n
  colnames(dane[[2]]) <- 1:n
  res <- melt(dane)
  res <- cbind(res, geny = rep(rownames(counts)[match(gene.list, rownames(counts1))], 2*n))
  p <- ggplot(data=res, aes(x=variable, y=value, group=L1, color=L1)) +
    geom_line() + theme_bw() + facet_wrap( ~ geny, scales = "free", ncol=4) + 
    theme(legend.position = "bottom", legend.title=element_blank()) +
    ylab("Counts") + xlab("Samples")
  print(p)
}