obliczenia.all.methods <- function(county, n.mapper = 3, design, p.val=0.05, method = "edgeR")
{
  require("edgeR")
  require("DESeq") 
  require("DESeq2")
  require("limma")
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
           b <- match(rownames(lrt1$table),rownames(county))
           out <- data.frame(ID=rownames(lrt1$table)[b], stat=lrt1$table$LR[b], padj=lrt1$table$PValue[b])
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
           res1 <- topTable(voom.fitbayes1, coef = 2:n.groups, number=nrow(voom.fitbayes1), adjust.method="BH", sort.by="none")
           b <- match(rownames(res1),rownames(county))
           if(length(unique(grupy)) > 2) out <- data.frame(ID=rownames(res1)[b], stat=res1$F[b], padj=res1$adj.P.Val[b]) else
             out <- data.frame(ID=rownames(res1)[b], stat=res1$B[b], padj=res1$adj.P.Val[b])
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
           b <- match(rownames(res1),rownames(county))
           if(length(unique(grupy)) > 2) out <- data.frame(ID=rownames(res1)[b], stat=res1$F[b], padj=res1$adj.P.Val[b]) else
             out <- data.frame(ID=rownames(res1)[b], stat=res1$B[b], padj=res1$adj.P.Val[b])
         },
         DESeq = 
         {
           design1 <- data.frame(grupy = grupy, mapery = mapery)
           cdsFull = newCountDataSet(round(county), design1)
           cdsFull = estimateSizeFactors( cdsFull )
           cdsFull = estimateDispersions( cdsFull )
           fit0 = fitNbinomGLMs( cdsFull, count ~  grupy + mapery )
           fit2 = fitNbinomGLMs( cdsFull, count ~  mapery )

           pvalsGLM2 = nbinomGLMTest(fit0, fit2)
           padjGLM2 = p.adjust( pvalsGLM2, method="BH" )
           b <- match(rownames(fit2),rownames(county))
           out <- data.frame(ID = rownames(fit0)[b], stat = fit0$`(Intercept)`[b], padj = padjGLM2[b])
         })
  
  res <- out
  res
}

obliczenia.all.methods1<-function(data, klas, method = "edgeR")
{
  library(limma)
  library(edgeR)
  library(DESeq)
  library(DESeq2)
  isexpr <- rowSums(data > 5) >= length(klas)/3
  data <- data[isexpr,]
  
  switch(method,
         DESeq = 
         { 
            colData <- data.frame(condition=as.factor(klas))
            rownames(colData) <- colnames(data)
            dds <- DESeqDataSetFromMatrix(countData = round(data), colData = colData, design = ~ condition) 
            dds <- DESeq(dds)
            res <- results(dds)
            ID <- rownames(res)
            res <- as.data.frame(res@listData)
            porownania <- data.frame(ID = ID, stat = res$stat, padj = res$padj)
         },
         edgeR = 
         { 
            y <- DGEList(data, group = klas)
            y <- calcNormFactors(y)
            if (length(unique(klas)) > 2)
            {
              design1 <- model.matrix(~ klas)
              y1 <- estimateDisp(y, design1)
              fit <- glmFit(y1, design1)
              lrt1 <- glmLRT(fit, coef=2:length(unique(klas)))
              b <- match(rownames(lrt1$table),rownames(data))
              porownania <- data.frame(ID=rownames(lrt1$table)[b], stat=lrt1$table$LR[b], padj=lrt1$table$PValue[b])
            } else {
              y <- estimateCommonDisp(y)
              y <- estimateTagwiseDisp(y)
              et <- exactTest(y)
              et <- et$table
              padj <- p.adjust(et[,3], method='BH')
              res <- data.frame(et ,padj)
              b <- match(rownames(data), rownames(res))
              porownania <- data.frame(ID = rownames(res)[b], stat = res$logFC[b], padj = res$padj[b])
            }
         },
         limma1 = 
         {
            #### limma + voom
            nf1 <- calcNormFactors(data, method = "TMM")
            voom.data <- voom(data, design = model.matrix(~factor(klas)),
                              lib.size = colSums(data) * nf1)
            voom.data$genes <- rownames(data)
            voom.fitlimma1 <- lmFit(voom.data, design = model.matrix(~factor(klas)))
            voom.fitbayes1 <- eBayes(voom.fitlimma1)
            res <- topTable(voom.fitbayes1, coef = 2:length(unique(klas)), number=nrow(voom.fitbayes1), adjust.method="BH")
            b <- match(rownames(data), rownames(res))
            if(length(unique(klas)) > 2) porownania <- data.frame(ID = rownames(res)[b], stat = res$F[b], padj = res$adj.P.Val[b]) else
            porownania <- data.frame(ID = rownames(res)[b], stat = res$B[b], padj = res$adj.P.Val[b])
         },
         limma2 = 
         {  
            #### limma + vst
            DESeq.cds1 <- newCountDataSet(countData = round(data), conditions = factor(klas))
            DESeq.cds1 <- estimateSizeFactors(DESeq.cds1)
            DESeq.cds1 <- estimateDispersions(DESeq.cds1, method = "blind", fitType = "local")
            DESeq.vst1 <- DESeq::getVarianceStabilizedData(DESeq.cds1)
            DESeq.vst.fitlimma1 <- lmFit(DESeq.vst1, design = model.matrix(~factor(klas)))
            DESeq.vst.fitbayes1 <- eBayes(DESeq.vst.fitlimma1)
            res <- topTable(DESeq.vst.fitbayes1, coef = 2:length(unique(klas)), adjust.method="BH", number=nrow(data))
            b <- match(rownames(data), rownames(res))
            if(length(unique(klas)) > 2) porownania <- data.frame(ID = rownames(res)[b], stat = res$F[b], padj = res$adj.P.Val[b]) else
            porownania <- data.frame(ID = rownames(res)[b], stat = res$B[b], padj = res$adj.P.Val[b])
         })
  res <- porownania
  res
}
