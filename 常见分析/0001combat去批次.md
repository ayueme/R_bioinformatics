今天演示下如何使用不同方法去除批次效应。

本期目录：

[toc]

## 准备数据

我们直接用`easyTCGA`下载结肠癌和直肠癌的转录组基因表达数据。


```r
library(easyTCGA)

# 1行代码搞定一切
getmrnaexpr(c("TCGA-COAD","TCGA-READ"))
```

然后就可以加载数据了。

先使用`tpm`数据做个演示。如果是基因表达芯片数据，也是和`tpm`数据的处理方法一样的。


```r
# 加载表达矩阵和临床信息
load(file = "output_mRNA_lncRNA_expr/TCGA-COAD_TCGA-READ_mrna_expr_tpm.rdata")
load(file = "output_mRNA_lncRNA_expr/TCGA-COAD_TCGA-READ_lncrna_expr_tpm.rdata")
load(file = "output_mRNA_lncRNA_expr/TCGA-COAD_TCGA-READ_clinical.rdata")
```

表达矩阵进行`log2`转换，然后提取临床信息中的样本类型和`project`类型。

`easyTCGA`下载的是最新的官网数据，理论上这种数据是最全的，但是也都是未经整理的，很多临床信息的名字和`xena`这些网站下载的并不一样，需要自己多探索哈。`easyTCGA`包会保存临床信息到csv文件里的，打开一看便知。

TCGA的官网是：https://portal.gdc.cancer.gov/。建议初学者熟悉下官网的数据，再结合第三方网站的数据，一起探索。

我们的临床信息中`project_id`就是我们需要的批次信息，一个是`TCGA-COAD`，一个是`TCGA-READ`。

`sample_type`是样本类型，但是太详细了，我们其实只要`normal/tumor`即可。


```r
exprset <- log2(mrna_expr_tpm+0.1)
exprset_lnc <- log2(lncrna_expr_tpm+0.1)
dim(mrna_expr_tpm)
## [1] 19938   701
table(clin_info$project_id)
## 
## TCGA-COAD TCGA-READ 
##       524       177
table(clin_info$sample_type)
## 
##          Metastatic       Primary Tumor     Recurrent Tumor Solid Tissue Normal 
##                   1                 647                   2                  51
```

可以看到一共有19938个编码基因，701个样本。其中`TCGA-COAD`有524个，`TCGA-READ`有177个。

然后是不同的样本类型，其实太详细了，我们只要分为`normal`和`tumor`即可。


```r
# 修改下样本类型，project_id就是批次信息，不改了
clin_info$sample_type <- ifelse(as.numeric(substr(clin_info$barcode,14,15))<10,"tumor","normal")
table(clin_info$sample_type)
## 
## normal  tumor 
##     51    650
table(clin_info$project_id)
## 
## TCGA-COAD TCGA-READ 
##       524       177

head(clin_info[,c("project_id","sample_type")])
##                              project_id sample_type
## TCGA-D5-6540-01A-11R-1723-07  TCGA-COAD       tumor
## TCGA-AA-3525-11A-01R-A32Z-07  TCGA-COAD      normal
## TCGA-AA-3525-01A-02R-0826-07  TCGA-COAD       tumor
## TCGA-AA-3815-01A-01R-1022-07  TCGA-COAD       tumor
## TCGA-D5-6923-01A-11R-A32Z-07  TCGA-COAD       tumor
## TCGA-G4-6322-01A-11R-1723-07  TCGA-COAD       tumor
```

样本顺序和表达矩阵完全一致，方便各种增删改查！


```r
identical(clin_info$barcode, colnames(exprset))
## [1] TRUE
identical(clin_info$barcode, colnames(exprset_lnc))
## [1] TRUE
```

## 数据探索

`easyTCGA`更加侧重下载和整理数据，只是顺便带了差异分析和批量生存分析，下游的各种分析已经有非常多的R包实现了。

首先我们用非常好用的`tinyarray`探索下还没有进行批次矫正的数据。这里就用mRNA的数据作为示例了。


```r
library(tinyarray)
```

### pca

首先是`pca`可视化，1行代码即可，极大节省你自己写代码的时间。


```r
draw_pca(exp = exprset, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-7](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-7-176440171.png)

可以看到`normal`和`tumor`重叠的部分还是蛮多的，数据质量不是非常好，但基本上还是能分得开，也不算太差。如果是非常好的数据应该是分的很开那种，当然也不是绝对的。

我们换个思路，把分组信息换成批次信息再看一看。


```r
draw_pca(exp = exprset, group_list = factor(clin_info$project_id))
```

![plot of chunk unnamed-chunk-8](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-8-176440171.png)

其实这样看还是可以的，`COAD`和`READ`基本混在一起，没有呈现出明显的批次。

### 层次聚类

探索完`pca`之后，我们再使用聚类分析探索下数据。

因为列名太长了，不方便展示，所以我改成了数字了，非常不推荐哈，这样就不能清晰地知道哪个样本不合群了。

这里的聚类还加了个颜色条，展示`normal/tumor`分组，我觉的是个非常不错的技巧。

参考历史推文：[聚类分析可视化之dendextend](https://mp.weixin.qq.com/s/YPz9pL2fDQ7cFX9_Fyw1yg)


```r
suppressPackageStartupMessages(library(dendextend))

tmp <- exprset
colnames(tmp) <- 1:ncol(tmp)

h.clust <- hclust(dist(scale(t(tmp))))
h.clust <- as.dendrogram(h.clust)

sample_colors <- ifelse(clin_info$sample_type == "tumor","red","green")

# 留足画图空间，防止颜色条显示不出来
par(mar=c(15,1,1,1))
plot(h.clust)
colored_bars(colors = sample_colors, dend = h.clust)
```

![plot of chunk unnamed-chunk-9](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-9-176440171.png)

简单的数据探索过程就到这里，如果你发现有明显的异常样本，可以手动删除，我们这里为了省事，就不做这一步了。过段时间再写个专门的推文讨论这个问题。

## combat

首先用`sva`包的`ComBat`函数去批次：


```r
library(sva)
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-41. For overview type 'help("mgcv-package")'.
## Loading required package: genefilter
## Loading required package: BiocParallel

expr_combat <- ComBat(dat = exprset, batch = clin_info$project_id)
## Found 620 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
## Found2batches
## Adjusting for0covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data

# 查看去除批次后的数据
expr_combat[1:4,1:4]
##        TCGA-D5-6540-01A-11R-1723-07 TCGA-AA-3525-11A-01R-A32Z-07
## MT-CO2                     14.50611                     14.11911
## MT-CO3                     14.43620                     13.98267
## MT-ND4                     14.36872                     13.31916
## MT-CO1                     14.61993                     13.99882
##        TCGA-AA-3525-01A-02R-0826-07 TCGA-AA-3815-01A-01R-1022-07
## MT-CO2                     14.04289                     14.85734
## MT-CO3                     14.20102                     14.67482
## MT-ND4                     13.51694                     14.52215
## MT-CO1                     13.62324                     14.51841
expr_combat <- as.data.frame(expr_combat)

# 再画图看看
draw_pca(exp = expr_combat, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-10](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-10-176440171.png)

效果好像不太明显，和没去除批次效应之前的图基本差不多。


```r
# lncRNA的也做一下，保存下数据
expr_lnc_combat <- ComBat(dat = exprset_lnc, batch = clin_info$project_id)
expr_lnc_combat <- as.data.frame(expr_lnc_combat)
save(expr_combat,expr_lnc_combat,clin_info,file = "step1_output.rdata")
```


`ComBat`里面有个`mod`选项，可以用来指定感兴趣的分组（这里就是normal和tumor），告诉函数不要把本来的分组信息给整没了。

我们再试试：


```r
mod <- model.matrix(~factor(clin_info$sample_type))
expr_combat <- ComBat(dat = exprset, batch = clin_info$project_id,mod = mod)
## Found 620 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
## Found2batches
## Adjusting for1covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
draw_pca(exp = expr_combat, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-12](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-12-176440171.png)

还是差不多哈...

## removeBatchEffect

然后我们再尝试下`limma`包的`removeBatchEffect`函数去批次。使用起来也是一模一样的简单。


```r
library(limma)

mod <- model.matrix(~factor(clin_info$sample_type))
expr_rbe <- removeBatchEffect(exprset, batch = clin_info$project_id, design=mod)

# 继续画图
draw_pca(exp = expr_rbe, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-13](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-13-176440171.png)

好像比`ComBat`差一点。

## count数据

如果是基因表达芯片数据，那么思路和上面一模一样，也是可以用`sva`和`limma`即可。

如果是count数据（测序数据），则可以采用下面的思路。

到底是用`count,fpkm,tpm`哪一种？也不用太纠结，我们已经写过很多推文介绍了：

- [DESeq2差异分析及VST变换的探索](https://mp.weixin.qq.com/s/CBznByKNGwPEKIKM5U0Oyw)
- [count、tpm、fpkm等表达量差异](https://mp.weixin.qq.com/s/aff-AX9aA2tSDa2zbB8ZRQ)
- [批量生存分析(logrank和单因素COX)](https://mp.weixin.qq.com/s/o-gCc_1B9SQmNFrG-I6yAQ)

首先加载数据。


```r
load(file = "output_mRNA_lncRNA_expr/TCGA-COAD_TCGA-READ_mrna_expr_counts.rdata")
identical(clin_info$barcode, colnames(mrna_expr_counts))
## [1] TRUE
```

首先还是`sva`包，用其中的`ComBat_seq`函数即可。专门针对`count`数据。


```r
expr_count_combat <- ComBat_seq(counts = as.matrix(mrna_expr_counts), 
                                batch = clin_info$project_id,
                                group = clin_info$sample_type
                                )
## Found 2 batches
## Using full model in ComBat-seq.
## Adjusting for 1 covariate(s) or covariate level(s)
## Estimating dispersions
## Fitting the GLM model
## Shrinkage off - using GLM estimates for parameters
## Adjusting the data
expr_count_combat[1:4,1:4]
##        TCGA-D5-6540-01A-11R-1723-07 TCGA-AA-3525-11A-01R-A32Z-07
## MT-CO1                       474829                       441211
## MT-ND4                       371029                       249643
## MT-CO2                       201649                       211390
## MT-CO3                       217204                       221701
##        TCGA-AA-3525-01A-02R-0826-07 TCGA-AA-3815-01A-01R-1022-07
## MT-CO1                       139679                       219156
## MT-ND4                       120710                       202674
## MT-CO2                        84905                       126534
## MT-CO3                       106881                       126317
```

画图看看去除批次效应之前和之后的pca：


```r
draw_pca(exp = mrna_expr_counts, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-16](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-16-176440171.png)

```r
draw_pca(exp = expr_count_combat, group_list = factor(clin_info$sample_type))
```

![plot of chunk unnamed-chunk-16](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-16-276440171.png)

差别不是很大，因为我们用的`TCGA-COAD`和`TCGA-READ`并没有很明显的批次效应，所以这几种方法用下来差别都不是很大。

也可以使用`DESeq2`，但是这个包是在做差异分析时顺便帮你把批次效应去除，不能单独去除批次效应。


```r
library(DESeq2)

dds1 <- DESeqDataSetFromMatrix(countData = mrna_expr_counts, 
                               colData = clin_info, 
                               design = ~ sample_type+project_id # 批次效应写在这里即可
                               )
```

后面就是做差异分析的步骤了，就不再演示了。可以参考咱们的历史推文：

[DESeq2差异分析及VST变换的探索](https://mp.weixin.qq.com/s/CBznByKNGwPEKIKM5U0Oyw)
