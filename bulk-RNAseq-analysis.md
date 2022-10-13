
# Introduction


Attach the libraries.

```{r warning=FALSE, message=FALSE}
library(SummarizedExperiment)
library(edgeR)
library(org.Hs.eg.db)
```

Load the data.

```{r}
full_brain <- get(load('rse_gene_brain_9_scaled.RData'))
full_colon <- get(load('rse_gene_colon_5_scaled.RData'))
full_heart <- get(load('rse_gene_heart_6_scaled.Rdata'))
```

Create three objects in order to easily retrieve information about the data.

```{r}
gtex_brain <- full_brain[,c(9,10,1)]
gtex_heart <- full_heart[,c(4,5,6)]
gtex_colon <- full_colon[,c(8,9,10)]
```

Find the samples id, these will be useful in GTEx Portal.

```{r}
colData(gtex_brain)$sampid
```
```{r}
[1] "GTEX-13X6K-0011-R6a-SM-5P9K6" "GTEX-X4XX-0011-R2A-SM-3P623" 
[3] "GTEX-T5JC-0011-R11A-SM-5S2RX"

```{r}
colData(gtex_heart)$sampid
```
```{r}
```

```{r}
colData(gtex_colon)$sampid
```
```{r}
```

Transform into dataframe and select the columns.

```{r}
br <- as.data.frame(assay(full_brain))[,c(9,10,1)]
col <- as.data.frame(assay(full_colon))[,c(8,9,10)]
hea <- as.data.frame(assay(full_heart))[,c(4,5,6)]
```

# Preprocessing

## Brain

```{r}
brain <- as.data.frame(full_brain@rowRanges)
head(brain)
```
```{r}
```

Find the total of the initial reads.

```{r}
totalreads_br <- colSums(br)
totalreads_br
```
```{r}
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_brain <- br[row.names(brain[brain$bp_length < 200,]),]
totalreads_shortRnas_brain <- colSums(short_genes_brain)
totalreads_shortRnas_brain
```
```{r}
```

Remove short genes from the dataset.

```{r}
brain <- brain[brain$bp_length >= 200,]
br <- br[rownames(brain),]
```

Find the mitochondrial genes and count their reads.

```{r}
mitho_genes_br <- br[row.names(brain[brain$seqnames == "chrM",]),]
totalreads_mithobr <- colSums(mitho_genes_br)
totalreads_mithobr
```
```{r}
```

Final brain data.

```{r}
brain <- brain[!(brain$gene_id %in% rownames(mitho_genes_br)),]
br <- br[rownames(brain),]
```

## Colon

```{r}
colon <- as.data.frame(full_colon@rowRanges)
head(colon)
```
```{r}
```

Find the total of the initial reads.

```{r}
totalreads_col <- colSums(col)
totalreads_col
```
```{r}
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_colon <- col[row.names(colon[colon$bp_length < 200,]),]
totalreads_shortRnas_colon <- colSums(short_genes_colon)
totalreads_shortRnas_colon
```
```{r}
```

Remove short genes from the dataset.

```{r}
colon <- colon[colon$bp_length >= 200,]
col <- col[rownames(colon),]
```

Find the mitochondrial genes and count their reads.

```{r}
mitho_genes_col <- col[row.names(colon[colon$seqnames == "chrM",]),]
totalreads_mithocol <- colSums(mitho_genes_col)
totalreads_mithocol
```
```{r}
```

Final colon data.

```{r}
colon <- colon[!(colon$gene_id %in% rownames(mitho_genes_col)),]
col <- col[rownames(colon),]
```

## Heart

```{r}
heart <- as.data.frame(full_heart@rowRanges)
head(heart)
```
```{r}
```

Find the total of the initial reads.

```{r}
totalreads_hea <- colSums(hea)
totalreads_hea
```
```{r}
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_heart <- hea[row.names(heart[heart$bp_length < 200,]),]
totalreads_shortRnas_heart <- colSums(short_genes_heart)
totalreads_shortRnas_heart
```
```{r}
```

Remove short genes from the dataset.

```{r}
heart <- heart[heart$bp_length >= 200,]
hea <- hea[rownames(heart),]
```
```{r}
```

Find the mitochondrial genes and count their reads.

```{r}
mitho_genes_hea <- hea[row.names(heart[heart$seqnames == "chrM",]),]
totalreads_mithohea <- colSums(mitho_genes_hea)
totalreads_mithohea
```
```{r}
```

Final heart data.

```{r}
heart <- heart[!(heart$gene_id %in% rownames(mitho_genes_hea)),]
hea <- hea[rownames(heart),]
```

## Join the counts table

```{r}
counts_table <- cbind(br,hea,col)
head(counts_table)
```
```{r}
```

# DE gene analysis

Build DGEList element.

```{r}
DE_object <- DGEList(counts = counts_table)
head(DE_object)
```
```{r}
```

Labeling of the samples. I've decided to call "kolon" because the design matrix is created in alphabetical order.

```{r}
tissue <- as.factor(c("Brain","Brain","Brain","Heart","Heart","Heart","Kolon","Kolon","Kolon"))
DE_object$samples$group <- tissue
head(DE_object)
```
```{r}
```

## Filtering and normalization

Filter out genes with low counts.

```{r}
table(rowSums(DE_object$counts==0)==9)
```
```{r}
```
```{r}
keep.exprs <- filterByExpr(DE_object, group = tissue)
DE_object <- DE_object[keep.exprs,, keep.lib.sizes=FALSE]
```

Normalize the counts with TMM method.

```{r}
DE_object <- calcNormFactors(DE_object, method = "TMM")
head(DE_object) 
```
```{r}
```

Distribution of normalized log-cpm values across samples.

```{r, fig.height=8}
logcpm <- cpm(DE_object, log=TRUE)
boxplot(logcpm, outline = F, cex.axis= 0.5, col=rep(c("red","orange","green"),each=3), ylim = c(-5,20))
legend(0.3,20, c("Brain","Heart","Colon"), fill=c("red","orange","green"))
```
![image_1](https://github.com/lucaprocaccio/Bulk-RNAseq/blob/main/figure-html/image_01.png)

## Design matrix

Design of the linear model. No base case, hence no intercept, we have different tissues.

```{r}
design <- model.matrix(~ 0+group, data = DE_object$samples)
colnames(design) <- levels(DE_object$samples$group)
design
```
```{r}
```

## Exploratory analysis

Plot the samples labeled by group.
MDS is a dimensionality reduction technique in which are employed the fold ratio values among the three samples. The samples near in the space are similar between them.

```{r, fig.height=8, fig.width=12}
logcpm <- cpm(DE_object, log=TRUE)
plotMDS(logcpm, labels=rep(c("Brain","Heart","Colon"),each=3))
```
![image_2](https://github.com/lucaprocaccio/Bulk-RNAseq/blob/main/figure-html/image_02.png)

We plot also the samples labeled by sample id. 

```{r, fig.height=8, fig.width=12}
plotMDS(logcpm, labels= colnames(counts_table), col = rep(c("red","orange","green"),each=3))
legend(3.3,3.3, c("Brain","Heart","Colon"), fill=c("red","orange","green"))
```
![image_3](https://github.com/lucaprocaccio/Bulk-RNAseq/blob/main/figure-html/image_03.png)

Estimate the NB dispersion.

```{r}
DE_object <- estimateDisp(DE_object, design)
plotBCV(DE_object)
```
![image_4](https://github.com/lucaprocaccio/Bulk-RNAseq/blob/main/figure-html/image_04.png)

Common dispersion.

```{r}
DE_object$common.dispersion
```
```{r}
```

## Generalized linear model

Fit to generalized linear model. We use a generalized linear model because we assume a negative binomial distribution and not the normal distribution used in the classical linear model.

```{r}
fit <- glmQLFit(DE_object, design)
```

We make comparisons between tissues, then we use multiple testing and decision tests in order to select the right p-value. I've decided to use the default thresholds FDR=0.05 and no threshold on the log-fold ratio because I obtained too few genes in Heart vs Colon.

### Tissue 1 with tissue 2 (Brain and Heart).

```{r}
qlf.HvsB <- glmQLFTest(fit, contrast=c(-1,1,0))
topTags(qlf.HvsB)
```
```{r}
```

```{r}
FDR_HvsB <- p.adjust(qlf.HvsB$table$PValue, method = "BH")
sum(FDR_HvsB < 0.05)
```
```{r}
```
```{r}
summary(decideTests(qlf.HvsB, p.value = 0.05))

```
```{r}
```

```{r}
deg.HvsB <- topTags(qlf.HvsB, n=20000, adjust.method = "BH", sort.by="PValue", p.value = 0.05)$table
head(deg.HvsB)
```
```{r}
```
```{r}
up_genes_HvsB <- row.names(deg.HvsB[deg.HvsB$logFC > 0,])
down_genes_HvsB <- row.names(deg.HvsB[deg.HvsB$logFC < 0,])
```

### Tissue 1 with tissue 3 (Brain and Colon).

```{r}
qlf.CvsB <- glmQLFTest(fit, contrast=c(-1,0,1))
topTags(qlf.CvsB)
```
```{r}
```

```{r}
FDR_CvsB <- p.adjust(qlf.CvsB$table$PValue, method="BH")
sum(FDR_CvsB< 0.05)
```
```{r}
```
```{r}
summary(decideTests(qlf.CvsB, p.value=0.05))
```

```{r}
deg.CvsB <- topTags(qlf.CvsB, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)$table
head(deg.CvsB)
```
```{r}
```
```{r}
up_genes_CvsB <- row.names(deg.CvsB[deg.CvsB$logFC > 0,])
down_genes_CvsB<- row.names(deg.CvsB[deg.CvsB$logFC < 0,])
```

### Tissue 2 with tissue 3 (Heart and Colon).

```{r}
qlf.CvsH <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.CvsH)
```
```{r}
```

```{r}
FDR_CvsH <- p.adjust(qlf.CvsH$table$PValue, method="BH")
sum(FDR_CvsH< 0.05)
```
```{r}
```
```{r}
summary(decideTests(qlf.CvsH,p.value = 0.05))
```
```{r}
```

```{r}
deg.CvsH<- topTags(qlf.CvsH, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)$table
head(deg.CvsH)
```
```{r}
```
```{r}
up_genes_CvsH <- row.names(deg.CvsH[deg.CvsH$logFC > 0,])
down_genes_CvsH<- row.names(deg.CvsH[deg.CvsH$logFC < 0,])
```

## Extraction of the genes

The genes up_regulated in a tissue with respect to the other are the down-regulated in the second with respect to the first.

```{r}
up_genes_BvsH<-down_genes_HvsB
down_genes_BvsH<-up_genes_HvsB

up_genes_BvsC<-down_genes_CvsB
down_genes_BvsC<-up_genes_CvsB

up_genes_HvsC <- down_genes_CvsH
down_genes_HvsC<-up_genes_CvsH
```

Create the vectors for up and down regulated genes in all the tissues.

```{r}
up.brain<-up_genes_BvsH[up_genes_BvsH %in% up_genes_BvsC]
down.brain<-down_genes_BvsH[down_genes_BvsH %in% down_genes_BvsC]

up.heart<-up_genes_HvsB[up_genes_HvsB %in% up_genes_HvsC]
down.heart<-down_genes_HvsB[down_genes_HvsB %in% down_genes_HvsC]

up.colon<-up_genes_CvsB[up_genes_CvsB %in% up_genes_CvsH]
down.colon<-down_genes_CvsB[down_genes_CvsB %in% down_genes_CvsH]
```

### Question 9.

Find the gene which is over-expressed in 1vs2 and 1vs3 with the lowest FDR.

```{r}
up.complete <- deg.HvsB[up.brain,]
gene_9 <- up.complete[which.min(up.complete$FDR),]
gene_9
```
```{r}
```

Create the .csv files. I merged all the obtained .csv files in one unique .xlsx file.

```{r, eval=FALSE}
write.csv(rownames(deg.CvsB),"genes_DE_CvsB.csv")
write.csv(rownames(deg.HvsB), "genes_DE_HvsB.csv")
write.csv(rownames(deg.CvsH), "genes_DE_CvsH.csv")
write.csv(up.brain, "genes_up_brain.csv")
write.csv(down.brain, "genes_down_brain.csv")
write.csv(up.colon, "genes_up_colon.csv")
write.csv(down.colon, "genes_down_colon.csv")
write.csv(up.heart, "genes_up_heart.csv")
write.csv(down.heart, "genes_down_heart.csv")
```

# Enrichment analysis

Create the vectors to insert in Enrichr.

```{r, eval=FALSE}
ensembl.ID<- gsub("\\..*", "", up.brain)
up.brain.ens <- as.data.frame(mapIds(org.Hs.eg.db,keys=ensembl.ID, keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
up.brain.ens.labels<-up.brain.ens$`mapIds(org.Hs.eg.db, keys = ensembl.ID, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")`

ensembl.ID<- gsub("\\..*", "", up.heart)
up.heart.ens <- as.data.frame(mapIds(org.Hs.eg.db,keys=ensembl.ID, keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
up.heart.ens.labels<-up.heart.ens$`mapIds(org.Hs.eg.db, keys = ensembl.ID, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")`

ensembl.ID<- gsub("\\..*", "", up.colon)
up.colon.ens <- as.data.frame(mapIds(org.Hs.eg.db,keys=ensembl.ID, keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
up.colon.ens.labels<-up.colon.ens$`mapIds(org.Hs.eg.db, keys = ensembl.ID, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")`


write.table(up.brain.ens.labels[is.na(up.brain.ens.labels)==F],"up.brain.labels.txt",col.names = F,row.names = F,quote=F)
write.table(up.heart.ens.labels[is.na(up.heart.ens.labels)==F],"up.heart.labels.txt",col.names = F,row.names = F,quote=F)
write.table(up.colon.ens.labels[is.na(up.colon.ens.labels)==F],"up.colon.labels.txt",col.names = F,row.names = F,quote=F)
```
