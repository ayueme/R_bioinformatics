上一篇文章里面简单学习了一下表达矩阵的提取，顺便探索了一下`SummarizedExperiment`对象。

今天学习下用`TCGAbiolinks`做差异分析。

## 加载R包和数据


```r
rm(list = ls())

library(SummarizedExperiment)
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## Loading required package: GenomicRanges
## Loading required package: stats4
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following object is masked from 'package:grDevices':
## 
##     windows
## Loading required package: GenomeInfoDb
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
library(TCGAbiolinks)
```

首先是查询、下载、整理，但是这一步我们在之前已经做好了，直接加载就可以了！

下载方法：[1.新版TCGA数据库学习：批量下载数据](https://mp.weixin.qq.com/s/m8w1L4N2aXAIers_ZJvp_g)


```r
# 查询
query <- GDCquery(project = "TCGA-COAD",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts"
                    )
# 下载
GDCdownload(query, files.per.chunk = 100) #每次下载100个文件
  
# 整理
GDCprepare(query,save = T,save.filename = "TCGA-COAD_mRNA.Rdata")
```

我们已经下载好了，就直接加载即可：


```r
load("TCGA-mRNA/TCGA-COAD_mRNA.Rdata")
```

通常我们会区分mRNA和lncRNA，所以我们这里只选择mRNA即可，方法在上一篇也说过了，非常简单！*直接对`SummarizedExperiment`对象取子集即可！*


```r
se_mrna <- data[rowData(data)$gene_type == "protein_coding",]

se_mrna
## class: RangedSummarizedExperiment 
## dim: 19962 521 
## metadata(1): data_release
## assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(19962): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(521): TCGA-A6-5664-01A-21R-1839-07
##   TCGA-D5-6530-01A-11R-1723-07 ... TCGA-A6-2683-01A-01R-0821-07
##   TCGA-A6-2683-11A-01R-A32Z-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

## 数据预处理

在进行差异分析前进行一些预处理。


```r
dim(se_mrna)
## [1] 19962   521
```


首先根据spearman相关系数去除异常值：


```r
coad_coroutliers <- TCGAanalyze_Preprocessing(se_mrna,cor.cut = 0.7)
## Number of outliers: 0

dim(coad_coroutliers)
## [1] 19962   521
```

还会生成一张相关图：Array Array Intensity correlation (AAIC)：
![PreprocessingOutput](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/PreprocessingOutput.png)

接下来进行标准化：


```r
# normalization of genes
coadNorm <- TCGAanalyze_Normalization(
    tabDF = coad_coroutliers, 
    geneInfo =  geneInfoHT
)
## I Need about  127 seconds for this Complete Normalization Upper Quantile  [Processing 80k elements /s]
## Step 1 of 4: newSeqExpressionSet ...
## Step 2 of 4: withinLaneNormalization ...
## Step 3 of 4: betweenLaneNormalization ...
## Step 4 of 4: exprs ...
```

会使用`EDASeq`包中的方法：


```r
EDASeq::newSeqExpressionSet
EDASeq::withinLaneNormalization
EDASeq::betweenLaneNormalization
EDASeq::counts
```


```r
dim(coadNorm)
## [1] 19469   521
```


然后过滤掉低表达的基因：


```r
# quantile filter of genes
coadFilt <- TCGAanalyze_Filtering(
    tabDF = coadNorm,
    method = "quantile", 
    qnt.cut =  0.25
)
```

看看一通操作下来后还剩多少基因？


```r
dim(coadFilt)
## [1] 14600   521
```

从最开始的19962变成了现在的14600，过滤掉了5000+基因......

最后是根据肿瘤组织和正常组织进行分组：

这里我们只选择了实体瘤和部分正常组织。如果你想选择更多，只要在`typesample`参数中添加更多类型即可。

可选类型见下图，也是根据TCGA-barcode进行判断的：
![sample type](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-27_19-10-43.png)


```r
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
    barcode = colnames(coadFilt),
    typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
    barcode = colnames(coadFilt), 
    typesample = c("TP")
)
```

所以分组这一步我们还是自己搞定！就是根据barcode的第14,15位数，结合上面那张图判断，毫无难度。
![Snipaste_2022-07-27_19-54-42](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-27_19-54-42.png)


```r
# 小于10就是tumor
samplesTumor <- as.numeric(substr(colnames(coadFilt),14,15))<10
```

## 差异分析

非常简单，支持`edgeR`和`limma`两种方法，当然也可以无缝连接`DESeq2`进行差异分析！


```r
# Diff.expr.analysis (DEA)
coadDEGs <- TCGAanalyze_DEA(
    mat1 = coadFilt[,!samplesTumor], # normal矩阵
    mat2 = coadFilt[,samplesTumor], # tumor矩阵
    Cond1type = "Normal",
    Cond2type = "Tumor",
    fdr.cut = 0.01, 
    logFC.cut = 1,
    pipeline = "edgeR", # limma
    method = "glmLRT"
)
## Batch correction skipped since no factors provided
## ----------------------- DEA -------------------------------
## o 41 samples in Cond1type Normal
## o 480 samples in Cond2type Tumor
## o 14600 features as miRNA or genes
## This may take some minutes...
## ----------------------- END DEA -------------------------------
```

差异分析就做好了！结果非常完美，**同时提供了gene_name和gene_type，也就是说我们一开始不取子集也是可以的~~，最后再取也行！**
![Snipaste_2022-07-27_21-18-40](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-27_21-18-40.png)


## 使用DESeq2进行差异分析

连接`DESeq2`那真是太简单了，无缝衔接！！


```r
library(DESeq2)
```

直接把`SummarizedExperiment`对象传给`DESeqDataSet()`函数即可。

不过需要分组信息，这个需要我们手动制作一下。

### 制作分组信息

其实我们的对象中包含了`sample_type`这一列信息，就在`coldata`中，但是有点过于详细了。


```r
table(colData(se_mrna)$sample_type)
## 
##          Metastatic       Primary Tumor     Recurrent Tumor Solid Tissue Normal 
##                   1                 478                   1                  41
```

我们给它修改一下~


```r
new_type <- ifelse(colData(se_mrna)$sample_type == "Solid Tissue Normal",
                   "Normal",
                   "Tumor")
colData(se_mrna)$sample_type <- new_type

table(colData(se_mrna)$sample_type)
## 
## Normal  Tumor 
##     41    480
```

这样就把分组搞定了！skr!

### 差异分析

然后就是愉快的进行差异分析~


```r
ddsSE <- DESeqDataSet(se_mrna, design = ~ sample_type)
## renaming the first element in assays to 'counts'
ddsSE
## class: DESeqDataSet 
## dim: 19962 521 
## metadata(2): data_release version
## assays(6): counts stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(19962): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(521): TCGA-A6-5664-01A-21R-1839-07
##   TCGA-D5-6530-01A-11R-1723-07 ... TCGA-A6-2683-01A-01R-0821-07
##   TCGA-A6-2683-11A-01R-A32Z-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

先过滤，这里我们简单点，


```r
keep <- rowSums(counts(ddsSE)) >= 50
ddsSE <- ddsSE[keep,]

ddsSE
## class: DESeqDataSet 
## dim: 18820 521 
## metadata(2): data_release version
## assays(6): counts stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(18820): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(521): TCGA-A6-5664-01A-21R-1839-07
##   TCGA-D5-6530-01A-11R-1723-07 ... TCGA-A6-2683-01A-01R-0821-07
##   TCGA-A6-2683-11A-01R-A32Z-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

确定谁和谁比，我们设定*Tumor-Normal*。


```r
ddsSE$sample_type <- factor(ddsSE$sample_type, levels = c("Tumor","Normal"))
```

接下来就可以进行差异分析了：


```r
dds <- DESeq(ddsSE)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 1544 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing
res <- results(dds, contrast = c("sample_type","Tumor","Normal"))
res
## log2 fold change (MLE): sample_type Tumor vs Normal 
## Wald test p-value: sample_type Tumor vs Normal 
## DataFrame with 18820 rows and 6 columns
##                     baseMean log2FoldChange     lfcSE       stat      pvalue
##                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
## ENSG00000000003.15  5198.688       0.456284 0.1445636    3.15628 1.59793e-03
## ENSG00000000005.6     42.064       0.414580 0.3014120    1.37546 1.68989e-01
## ENSG00000000419.13  1743.704       0.849473 0.1134251    7.48929 6.92491e-14
## ENSG00000000457.14   463.909      -0.261646 0.0690276   -3.79046 1.50371e-04
## ENSG00000000460.17   328.546       1.272811 0.0940574   13.53229 1.00837e-41
## ...                      ...            ...       ...        ...         ...
## ENSG00000288658.1   5.317261    -1.45986802  0.274148 -5.3251160 1.00889e-07
## ENSG00000288660.1   4.648268     1.47152585  0.450554  3.2660389 1.09063e-03
## ENSG00000288669.1   0.143343    -0.13840184  1.339247 -0.1033431 9.17691e-01
## ENSG00000288674.1   3.201667     0.00548347  0.174731  0.0313824 9.74965e-01
## ENSG00000288675.1  15.048025     2.01434132  0.174631 11.5348224 8.80688e-31
##                           padj
##                      <numeric>
## ENSG00000000003.15 2.75622e-03
## ENSG00000000005.6  2.13290e-01
## ENSG00000000419.13 3.14116e-13
## ENSG00000000457.14 2.95220e-04
## ENSG00000000460.17 2.77449e-40
## ...                        ...
## ENSG00000288658.1  2.74461e-07
## ENSG00000288660.1  1.92231e-03
## ENSG00000288669.1  9.29645e-01
## ENSG00000288674.1  9.78866e-01
## ENSG00000288675.1  1.20917e-29
```

结果很棒，不过没有gene_symbol了，需要自己添加哦~





