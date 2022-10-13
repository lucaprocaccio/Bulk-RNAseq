
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
## [1] "GTEX-13X6K-0011-R6a-SM-5P9K6" "GTEX-X4XX-0011-R2A-SM-3P623" 
## [3] "GTEX-T5JC-0011-R11A-SM-5S2RX"
```

```{r}
colData(gtex_heart)$sampid
```
```{r}
## [1] "GTEX-13VXT-1126-SM-5LU3A" "GTEX-13O3P-1626-SM-5K7X3"
## [3] "GTEX-13N1W-0926-SM-5MR36"
```

```{r}
colData(gtex_colon)$sampid
```
```{r}
## [1] "GTEX-11ZU8-1126-SM-5EQ5K" "GTEX-ZAB5-1426-SM-5HL9D" 
## [3] "GTEX-11ILO-2026-SM-5N9CQ"
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
##                    seqnames     start       end  width strand
## ENSG00000000003.14     chrX 100627109 100639991  12883      -
## ENSG00000000005.5      chrX 100584802 100599885  15084      +
## ENSG00000000419.12    chr20  50934867  50958555  23689      -
## ENSG00000000457.13     chr1 169849631 169894267  44637      -
## ENSG00000000460.16     chr1 169662007 169854080 192074      +
## ENSG00000000938.12     chr1  27612064  27635277  23214      -
##                               gene_id bp_length   symbol
## ENSG00000000003.14 ENSG00000000003.14      4535   TSPAN6
## ENSG00000000005.5   ENSG00000000005.5      1610     TNMD
## ENSG00000000419.12 ENSG00000000419.12      1207     DPM1
## ENSG00000000457.13 ENSG00000000457.13      6883    SCYL3
## ENSG00000000460.16 ENSG00000000460.16      5967 C1orf112
## ENSG00000000938.12 ENSG00000000938.12      3474      FGR
```

Find the total of the initial reads.

```{r}
totalreads_br <- colSums(br)
totalreads_br
```
```{r}
## SRR1397094  SRR817686 SRR2166176 
##   36625433   37946538   33468019
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_brain <- br[row.names(brain[brain$bp_length < 200,]),]
totalreads_shortRnas_brain <- colSums(short_genes_brain)
totalreads_shortRnas_brain
```
```{r}
## SRR1397094  SRR817686 SRR2166176 
##      69732      58186      34110
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
## SRR1397094  SRR817686 SRR2166176 
##    6395981    2952326    2080589
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
##                    seqnames     start       end  width strand
## ENSG00000000003.14     chrX 100627109 100639991  12883      -
## ENSG00000000005.5      chrX 100584802 100599885  15084      +
## ENSG00000000419.12    chr20  50934867  50958555  23689      -
## ENSG00000000457.13     chr1 169849631 169894267  44637      -
## ENSG00000000460.16     chr1 169662007 169854080 192074      +
## ENSG00000000938.12     chr1  27612064  27635277  23214      -
##                               gene_id bp_length   symbol
## ENSG00000000003.14 ENSG00000000003.14      4535   TSPAN6
## ENSG00000000005.5   ENSG00000000005.5      1610     TNMD
## ENSG00000000419.12 ENSG00000000419.12      1207     DPM1
## ENSG00000000457.13 ENSG00000000457.13      6883    SCYL3
## ENSG00000000460.16 ENSG00000000460.16      5967 C1orf112
## ENSG00000000938.12 ENSG00000000938.12      3474      FGR
```

Find the total of the initial reads.

```{r}
totalreads_col <- colSums(col)
totalreads_col
```
```{r}
## SRR1398272 SRR1436583 SRR1420413 
##   38272378   38364393   38984105
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_colon <- col[row.names(colon[colon$bp_length < 200,]),]
totalreads_shortRnas_colon <- colSums(short_genes_colon)
totalreads_shortRnas_colon
```
```{r}
## SRR1398272 SRR1436583 SRR1420413 
##      38675      57902     149571
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
## SRR1398272 SRR1436583 SRR1420413 
##    1044689    4148166    1043438
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
##                    seqnames     start       end  width strand
## ENSG00000000003.14     chrX 100627109 100639991  12883      -
## ENSG00000000005.5      chrX 100584802 100599885  15084      +
## ENSG00000000419.12    chr20  50934867  50958555  23689      -
## ENSG00000000457.13     chr1 169849631 169894267  44637      -
## ENSG00000000460.16     chr1 169662007 169854080 192074      +
## ENSG00000000938.12     chr1  27612064  27635277  23214      -
##                               gene_id bp_length   symbol
## ENSG00000000003.14 ENSG00000000003.14      4535   TSPAN6
## ENSG00000000005.5   ENSG00000000005.5      1610     TNMD
## ENSG00000000419.12 ENSG00000000419.12      1207     DPM1
## ENSG00000000457.13 ENSG00000000457.13      6883    SCYL3
## ENSG00000000460.16 ENSG00000000460.16      5967 C1orf112
## ENSG00000000938.12 ENSG00000000938.12      3474      FGR
```

Find the total of the initial reads.

```{r}
totalreads_hea <- colSums(hea)
totalreads_hea
```
```{r}
## SRR1381693 SRR1467904 SRR1442648 
##   39624745   39167915   37718195 
```

Find genes with length < 200 and count their reads.

```{r}
short_genes_heart <- hea[row.names(heart[heart$bp_length < 200,]),]
totalreads_shortRnas_heart <- colSums(short_genes_heart)
totalreads_shortRnas_heart
```
```{r}
## SRR1381693 SRR1467904 SRR1442648 
##     116328     120191      67338 
```

Remove short genes from the dataset.

```{r}
heart <- heart[heart$bp_length >= 200,]
hea <- hea[rownames(heart),]
```
Find the mitochondrial genes and count their reads.

```{r}
mitho_genes_hea <- hea[row.names(heart[heart$seqnames == "chrM",]),]
totalreads_mithohea <- colSums(mitho_genes_hea)
totalreads_mithohea
```
```{r}
## SRR1381693 SRR1467904 SRR1442648 
##   13635931   23086070    8207035 
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
##                    SRR1397094 SRR817686 SRR2166176 SRR1381693 SRR1467904
## ENSG00000000003.14        500       234         93        133        141
## ENSG00000000005.5           1         1          1          1          4
## ENSG00000000419.12        584       511       1015        490        375
## ENSG00000000457.13        332       461        704        218        178
## ENSG00000000460.16        187       412        569         99         95
## ENSG00000000938.12        241       837         73        146        424
##                    SRR1442648 SRR1398272 SRR1436583 SRR1420413
## ENSG00000000003.14        282        486        371        507
## ENSG00000000005.5           5          4          2          3
## ENSG00000000419.12        758       1068        790        906
## ENSG00000000457.13        137        631        631        823
## ENSG00000000460.16         55        295        334        430
## ENSG00000000938.12        169        305         78        507
```

# DE gene analysis

Build DGEList element.

```{r}
DE_object <- DGEList(counts = counts_table)
head(DE_object)
```
```{r}
## An object of class "DGEList"
## $counts
##                    SRR1397094 SRR817686 SRR2166176 SRR1381693 SRR1467904
## ENSG00000000003.14        500       234         93        133        141
## ENSG00000000005.5           1         1          1          1          4
## ENSG00000000419.12        584       511       1015        490        375
## ENSG00000000457.13        332       461        704        218        178
## ENSG00000000460.16        187       412        569         99         95
## ENSG00000000938.12        241       837         73        146        424
##                    SRR1442648 SRR1398272 SRR1436583 SRR1420413
## ENSG00000000003.14        282        486        371        507
## ENSG00000000005.5           5          4          2          3
## ENSG00000000419.12        758       1068        790        906
## ENSG00000000457.13        137        631        631        823
## ENSG00000000460.16         55        295        334        430
## ENSG00000000938.12        169        305         78        507

## $samples
##            group lib.size norm.factors
## SRR1397094     1 30159720            1
## SRR817686      1 34936026            1
## SRR2166176     1 31353320            1
## SRR1381693     1 25872486            1
## SRR1467904     1 15961654            1
## SRR1442648     1 29443822            1
## SRR1398272     1 37189014            1
## SRR1436583     1 34158325            1
## SRR1420413     1 37791096            1
```

Labeling of the samples. I've decided to call "kolon" because the design matrix is created in alphabetical order.

```{r}
tissue <- as.factor(c("Brain","Brain","Brain","Heart","Heart","Heart","Kolon","Kolon","Kolon"))
DE_object$samples$group <- tissue
head(DE_object)
```
```{r}
## An object of class "DGEList"
## $counts
##                    SRR1397094 SRR817686 SRR2166176 SRR1381693 SRR1467904
## ENSG00000000003.14        500       234         93        133        141
## ENSG00000000005.5           1         1          1          1          4
## ENSG00000000419.12        584       511       1015        490        375
## ENSG00000000457.13        332       461        704        218        178
## ENSG00000000460.16        187       412        569         99         95
## ENSG00000000938.12        241       837         73        146        424
##                    SRR1442648 SRR1398272 SRR1436583 SRR1420413
## ENSG00000000003.14        282        486        371        507
## ENSG00000000005.5           5          4          2          3
## ENSG00000000419.12        758       1068        790        906
## ENSG00000000457.13        137        631        631        823
## ENSG00000000460.16         55        295        334        430
## ENSG00000000938.12        169        305         78        507

## $samples
##            group lib.size norm.factors
## SRR1397094 Brain 30159720            1
## SRR817686  Brain 34936026            1
## SRR2166176 Brain 31353320            1
## SRR1381693 Heart 25872486            1
## SRR1467904 Heart 15961654            1
## SRR1442648 Heart 29443822            1
## SRR1398272 Kolon 37189014            1
## SRR1436583 Kolon 34158325            1
## SRR1420413 Kolon 37791096            1
```

## Filtering and normalization

Filter out genes with low counts.

```{r}
table(rowSums(DE_object$counts==0)==9)
```
```{r}
## FALSE  TRUE 
## 42963  7718
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
## An object of class "DGEList"
## $counts
##                    SRR1397094 SRR817686 SRR2166176 SRR1381693 SRR1467904
## ENSG00000000003.14        500       234         93        133        141
## ENSG00000000419.12        584       511       1015        490        375
## ENSG00000000457.13        332       461        704        218        178
## ENSG00000000460.16        187       412        569         99         95
## ENSG00000000938.12        241       837         73        146        424
## ENSG00000000971.15        916      1385        214        376       1394
##                    SRR1442648 SRR1398272 SRR1436583 SRR1420413
## ENSG00000000003.14        282        486        371        507
## ENSG00000000419.12        758       1068        790        906
## ENSG00000000457.13        137        631        631        823
## ENSG00000000460.16         55        295        334        430
## ENSG00000000938.12        169        305         78        507
## ENSG00000000971.15      25192       1188        698       6813

## $samples
##            group lib.size norm.factors
## SRR1397094 Brain 29995666    1.2134888
## SRR817686  Brain 34874308    0.9992263
## SRR2166176 Brain 31194031    1.1230870
## SRR1381693 Heart 25852298    0.7050194
## SRR1467904 Heart 15942445    1.0031820
## SRR1442648 Heart 29394557    0.8625728
## SRR1398272 Kolon 37157633    0.9665450
## SRR1436583 Kolon 34126218    1.1542488
## SRR1420413 Kolon 37165532    1.0789192
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
##            Brain Heart Kolon
## SRR1397094     1     0     0
## SRR817686      1     0     0
## SRR2166176     1     0     0
## SRR1381693     0     1     0
## SRR1467904     0     1     0
## SRR1442648     0     1     0
## SRR1398272     0     0     1
## SRR1436583     0     0     1
## SRR1420413     0     0     1
## attr(,"assign")
## [1] 1 1 1
## attr(,"contrasts")
## attr(,"contrasts")$group
## [1] "contr.treatment"
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
## [1] 0.3929732
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
## Coefficient:  -1*Brain 1*Heart 
##                        logFC   logCPM        F       PValue          FDR
## ENSG00000277543.1   9.166968 5.780761 253.5775 2.044693e-08 0.0003313721
## ENSG00000183072.9  12.626649 6.083311 238.1193 2.766275e-08 0.0003313721
## ENSG00000136574.17 11.863806 5.995002 199.0226 6.523471e-08 0.0003333109
## ENSG00000141448.8   9.624067 5.124279 197.1823 6.818702e-08 0.0003333109
## ENSG00000077009.13 13.433005 8.749115 194.0774 7.354060e-08 0.0003333109
## ENSG00000139914.6   6.959024 4.865118 187.0840 8.756629e-08 0.0003333109
## ENSG00000105605.7  -9.457744 4.626350 182.2489 9.916241e-08 0.0003333109
## ENSG00000198796.6  11.368934 5.954487 177.8669 1.112984e-07 0.0003333109
## ENSG00000102683.7   7.072443 4.752410 172.2622 1.295279e-07 0.0003448034
## ENSG00000185739.13  9.640852 8.190995 162.9497 1.684441e-07 0.0003570690
```

```{r}
FDR_HvsB <- p.adjust(qlf.HvsB$table$PValue, method = "BH")
sum(FDR_HvsB < 0.05)
```
```{r}
## [1] 4955
```
```{r}
summary(decideTests(qlf.HvsB, p.value = 0.05))

```
```{r}
##        -1*Brain 1*Heart
## Down               2711
## NotSig            19003
## Up                 2244
```

```{r}
deg.HvsB <- topTags(qlf.HvsB, n=20000, adjust.method = "BH", sort.by="PValue", p.value = 0.05)$table
head(deg.HvsB)
```
```{r}
##                        logFC   logCPM        F       PValue          FDR
## ENSG00000277543.1   9.166968 5.780761 253.5775 2.044693e-08 0.0003313721
## ENSG00000183072.9  12.626649 6.083311 238.1193 2.766275e-08 0.0003313721
## ENSG00000136574.17 11.863806 5.995002 199.0226 6.523471e-08 0.0003333109
## ENSG00000141448.8   9.624067 5.124279 197.1823 6.818702e-08 0.0003333109
## ENSG00000077009.13 13.433005 8.749115 194.0774 7.354060e-08 0.0003333109
## ENSG00000139914.6   6.959024 4.865118 187.0840 8.756629e-08 0.0003333109
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
## Coefficient:  -1*Brain 1*Kolon 
##                        logFC    logCPM        F       PValue          FDR
## ENSG00000253293.4  10.091483  4.225605 174.3197 1.224467e-07 0.0008510269
## ENSG00000134871.17  6.474954  9.166218 172.9274 1.271858e-07 0.0008510269
## ENSG00000232814.2   6.290290  3.831022 160.5251 1.807898e-07 0.0008510269
## ENSG00000130635.15  6.140994  7.363288 152.9244 2.271962e-07 0.0008510269
## ENSG00000157103.10 -5.718709  5.978018 147.8553 2.661956e-07 0.0008510269
## ENSG00000131471.6   6.996863  7.132885 143.9021 3.022842e-07 0.0008510269
## ENSG00000234638.1  11.263170 10.147793 140.3377 3.399680e-07 0.0008510269
## ENSG00000143196.4   9.337146  7.135186 139.4288 3.504642e-07 0.0008510269
## ENSG00000163359.15  8.588134  8.774674 139.3316 3.516096e-07 0.0008510269
## ENSG00000134198.9   6.282824  4.916399 137.2473 3.772852e-07 0.0008510269
```

```{r}
FDR_CvsB <- p.adjust(qlf.CvsB$table$PValue, method="BH")
sum(FDR_CvsB< 0.05)
```
```{r}
## [1] 3642
```
```{r}
summary(decideTests(qlf.CvsB, p.value=0.05))
```
```{r}
##        -1*Brain 1*Kolon
## Down               2124
## NotSig            20316
## Up                 1518
```

```{r}
deg.CvsB <- topTags(qlf.CvsB, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)$table
head(deg.CvsB)
```
```{r}
##                        logFC   logCPM        F       PValue          FDR
## ENSG00000253293.4  10.091483 4.225605 174.3197 1.224467e-07 0.0008510269
## ENSG00000134871.17  6.474954 9.166218 172.9274 1.271858e-07 0.0008510269
## ENSG00000232814.2   6.290290 3.831022 160.5251 1.807898e-07 0.0008510269
## ENSG00000130635.15  6.140994 7.363288 152.9244 2.271962e-07 0.0008510269
## ENSG00000157103.10 -5.718709 5.978018 147.8553 2.661956e-07 0.0008510269
## ENSG00000131471.6   6.996863 7.132885 143.9021 3.022842e-07 0.0008510269
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
## Coefficient:  -1*Heart 1*Kolon 
##                         logFC    logCPM        F       PValue          FDR
## ENSG00000277543.1  -11.758788  5.780761 324.6686 6.197817e-09 0.0001484873
## ENSG00000183072.9  -11.602741  6.083311 224.6212 3.659434e-08 0.0004383636
## ENSG00000139914.6   -7.259598  4.865118 198.0326 6.680342e-08 0.0004453052
## ENSG00000136574.17 -10.783022  5.995002 184.0708 9.458819e-08 0.0004453052
## ENSG00000253293.4   11.423120  4.225605 177.0878 1.136389e-07 0.0004453052
## ENSG00000077009.13 -11.890852  8.749115 171.4459 1.324728e-07 0.0004453052
## ENSG00000156885.5  -11.605198  7.862392 169.2639 1.407489e-07 0.0004453052
## ENSG00000164708.5   -7.184849  7.626335 164.1210 1.628433e-07 0.0004453052
## ENSG00000239775.1   -6.714521  7.026682 163.1888 1.672822e-07 0.0004453052
## ENSG00000134571.10 -12.213731 10.637982 150.0877 2.481025e-07 0.0005944040
```

```{r}
FDR_CvsH <- p.adjust(qlf.CvsH$table$PValue, method="BH")
sum(FDR_CvsH< 0.05)
```
```{r}
## [1] 934
```
```{r}
summary(decideTests(qlf.CvsH,p.value = 0.05))
```
```{r}
##        -1*Heart 1*Kolon
## Down                630
## NotSig            23024
## Up                  304
```

```{r}
deg.CvsH<- topTags(qlf.CvsH, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)$table
head(deg.CvsH)
```
```{r}
##                         logFC   logCPM        F       PValue          FDR
## ENSG00000277543.1  -11.758788 5.780761 324.6686 6.197817e-09 0.0001484873
## ENSG00000183072.9  -11.602741 6.083311 224.6212 3.659434e-08 0.0004383636
## ENSG00000139914.6   -7.259598 4.865118 198.0326 6.680342e-08 0.0004453052
## ENSG00000136574.17 -10.783022 5.995002 184.0708 9.458819e-08 0.0004453052
## ENSG00000253293.4   11.423120 4.225605 177.0878 1.136389e-07 0.0004453052
## ENSG00000077009.13 -11.890852 8.749115 171.4459 1.324728e-07 0.0004453052
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
##                       logFC  logCPM        F       PValue          FDR
## ENSG00000105605.7 -9.457744 4.62635 182.2489 9.916241e-08 0.0003333109
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
