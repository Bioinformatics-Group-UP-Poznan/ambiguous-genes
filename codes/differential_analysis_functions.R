diffAnaAll <- function(counts, n.maper = 3, exp.design, p.val=0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  require("limma")
  if (ncol(counts) %% n.maper != 0) {
    stop("Number of columns is not a multiplication of number of mappers.")
  }
  n <- ncol(counts)/n.maper
  colnames(counts) <- paste0("mapper",rep(1:n.maper, each=n),"-",rep(1:n,times=n.maper))
  isexpr <- rowSums(counts > 5) >= n/3
  counts <- counts[isexpr,]
  groups <- as.factor(rep(exp.design, n.maper)) 
  mapers <- as.factor(rep(paste0("mapper",1:n.maper), each = n))
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
           b <- match(rownames(lrt1$table),rownames(counts))
           out <- data.frame(ID=rownames(lrt1$table)[b], stat=lrt1$table$LR[b], padj=lrt1$table$PValue[b])
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
           res1 <- topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           b <- match(rownames(res1),rownames(counts))
           if(length(unique(groups)) > 2) out <- data.frame(ID=rownames(res1)[b], stat=res1$F[b], padj=res1$adj.P.Val[b]) else
             out <- data.frame(ID=rownames(res1)[b], stat=res1$B[b], padj=res1$adj.P.Val[b])
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
           b <- match(rownames(res1),rownames(counts))
           if(length(unique(groups)) > 2) out <- data.frame(ID=rownames(res1)[b], stat=res1$F[b], padj=res1$adj.P.Val[b]) else
             out <- data.frame(ID=rownames(res1)[b], stat=res1$B[b], padj=res1$adj.P.Val[b])
         },
         DESeq = 
         {
           design1 <- data.frame(groups, mapers)
           cdsFull = newCountDataSet(round(counts), design1)
           cdsFull = estimateSizeFactors( cdsFull )
           cdsFull = estimateDispersions( cdsFull )
           fit0 = fitNbinomGLMs( cdsFull, count ~  groups + mapers )
           fit2 = fitNbinomGLMs( cdsFull, count ~  mapers )

           pvalsGLM2 = nbinomGLMTest(fit0, fit2)
           padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
           b <- match(rownames(fit2),rownames(counts))
           out <- data.frame(ID = rownames(fit0)[b], stat = fit0$`(Intercept)`[b], padj = padjGLM2[b])
         })
  
  res <- out
  res
}

diffAnaAll1 <- function(data, classes, method = "edgeR")
{
  library(limma)
  library(edgeR)
  library(DESeq)
  library(DESeq2)
  isexpr <- rowSums(data > 5) >= length(classes)/3
  data <- data[isexpr,]
  
  switch(method,
         DESeq = 
         { 
            colData <- data.frame(condition=as.factor(classes))
            rownames(colData) <- colnames(data)
            dds <- DESeqDataSetFromMatrix(countData = round(data), colData = colData, design = ~ condition) 
            dds <- DESeq(dds)
            res <- results(dds)
            ID <- rownames(res)
            res <- as.data.frame(res@listData)
            topTab <- data.frame(ID = ID, stat = res$stat, padj = res$padj)
         },
         edgeR = 
         { 
            y <- DGEList(data, group = classes)
            y <- calcNormFactors(y)
            if (length(unique(classes)) > 2)
            {
              design1 <- model.matrix(~ classes)
              y1 <- estimateDisp(y, design1)
              fit <- glmFit(y1, design1)
              lrt1 <- glmLRT(fit, coef=2:length(unique(classes)))
              b <- match(rownames(lrt1$table),rownames(data))
              topTab <- data.frame(ID=rownames(lrt1$table)[b], stat=lrt1$table$LR[b], padj=lrt1$table$PValue[b])
            } else {
              y <- estimateCommonDisp(y)
              y <- estimateTagwiseDisp(y)
              et <- exactTest(y)
              et <- et$table
              padj <- p.adjust(et[,3], method='BH')
              res <- data.frame(et ,padj)
              b <- match(rownames(data), rownames(res))
              topTab <- data.frame(ID = rownames(res)[b], stat = res$logFC[b], padj = res$padj[b])
            }
         },
         limma1 = 
         {
            #### limma + voom
            nf1 <- calcNormFactors(data, method = "TMM")
            voom.data <- voom(data, design = model.matrix(~factor(classes)),
                              lib.size = colSums(data) * nf1)
            voom.data$genes <- rownames(data)
            voom.fitlimma1 <- lmFit(voom.data, design = model.matrix(~factor(classes)))
            voom.fitbayes1 <- eBayes(voom.fitlimma1)
            res <- topTable(voom.fitbayes1, coef = 2:length(unique(classes)), number=nrow(voom.fitbayes1), adjust.method="BH")
            b <- match(rownames(data), rownames(res))
            if(length(unique(classes)) > 2) porownania <- data.frame(ID = rownames(res)[b], stat = res$F[b], padj = res$adj.P.Val[b]) else
              topTab <- data.frame(ID = rownames(res)[b], stat = res$B[b], padj = res$adj.P.Val[b])
         },
         limma2 = 
         {  
            #### limma + vst
            DESeq.cds1 <- newCountDataSet(countData = round(data), conditions = factor(classes))
            DESeq.cds1 <- estimateSizeFactors(DESeq.cds1)
            DESeq.cds1 <- estimateDispersions(DESeq.cds1, method = "blind", fitType = "local")
            DESeq.vst1 <- DESeq::getVarianceStabilizedData(DESeq.cds1)
            DESeq.vst.fitlimma1 <- lmFit(DESeq.vst1, design = model.matrix(~factor(classes)))
            DESeq.vst.fitbayes1 <- eBayes(DESeq.vst.fitlimma1)
            res <- topTable(DESeq.vst.fitbayes1, coef = 2:length(unique(classes)), adjust.method="BH", number=nrow(data))
            b <- match(rownames(data), rownames(res))
            if(length(unique(classes)) > 2) porownania <- data.frame(ID = rownames(res)[b], stat = res$F[b], padj = res$adj.P.Val[b]) else
              topTab <- data.frame(ID = rownames(res)[b], stat = res$B[b], padj = res$adj.P.Val[b])
         })
  res <- topTab
  res
}
