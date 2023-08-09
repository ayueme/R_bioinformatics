`TCGAbiolinks`可以进行甲基化分析，但是功能不如`ChAMP`强大，甲基化分析还是首推`ChAMP`包。

不过为了了解`TCGAbiolinks`包，里面关于甲基化分析的部分还是要学习一下。

主要是甲基化差异分析，甲基化的一些可视化，甲基化和转录组数据的联合作图。

## 加载数据

我们还是使用之前下载好的TCGA-COAD的甲基化β值矩阵。

数据下载见这篇：[使用TCGAbiolinks批量下载最新版TCGA数据库的各种组学数据！](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247492403&idx=1&sn=dcb748ca56e4e13181dedc310aaef5e2&scene=21#wechat_redirect)


```r
library(TCGAbiolinks)
suppressMessages(library(SummarizedExperiment))
load(file = "G:/tcga/TCGA-dnaMethy/COAD_METHY_beta.Rdata")
```

查看一下数据，现在这个β值矩阵还是一个`SummarizedExperiment`对象：


```r
beta.m <- data
beta.m
## class: RangedSummarizedExperiment 
## dim: 485577 352 
## metadata(1): data_release
## assays(1): ''
## rownames(485577): cg13869341 cg14008030 ... cg11478607 cg08417382
## rowData names(52): address_A address_B ... MASK_extBase MASK_general
## colnames(352): TCGA-D5-6530-01A-11D-1721-05
##   TCGA-AA-3660-11A-01D-1721-05 ... TCGA-A6-5664-01A-21D-1837-05
##   TCGA-D5-6533-01A-11D-1721-05
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

## 甲基化差异分析

β矩阵不能含有缺失值，所以先去除缺失值，使用缺失值插补方法进行插补也是可以的：


```r
beta.m <- subset(beta.m, subset = (rowSums(is.na(assay(beta.m))) == 0))
```

如果你的甲基化矩阵是直接使用`TCGAbiolinks`包整理好的`SummarizedExperiment`对象，那么这个甲基化差异分析就非常简单。

我们需要确定谁和谁进行相比，也就是要创建一个含有分组信息的列。

所有样本的信息都在`SummarizedExperiment`对象中的`colData`中，包括分组信息、生存信息、分期信息等，我们这里只要用到分组信息即可，所以先把`colData`简化一下。

如果你这里都看不懂，可以翻看之前的推文：

[新版TCGA数据库学习：表达矩阵提取（mRNA/lncRNA/counts/tpm/fpkm）](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247492462&idx=1&sn=0773b4ad0ecdd7e3cf206cd611249399&scene=21#wechat_redirect)
[新版TCGAbiolinks包学习：差异分析](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247492794&idx=1&sn=ff9e4aa0c76b0d4adadf0d124f2d90a9&scene=21#wechat_redirect)


```r
# 只要两列信息，其实只要sample_type一列也行
sample_info <- as.data.frame(colData(beta.m))[,c("barcode","sample_type")]

# sample_type改成normal和tumor
sample_info[sample_info == "Solid Tissue Normal"] <- "normal"
sample_info[sample_info != "normal"] <- "tumor"

table(sample_info$sample_type)
## 
## normal  tumor 
##     38    314

# 重新变成coldata
colData(beta.m) <- DataFrame(sample_info)
```

有了`SummarizedExperiment`对象和相应的分组信息就可以进行甲基化差异分析了。


```r
# 非常慢
res <- TCGAanalyze_DMC(data = beta.m,
                       groupCol = "sample_type", # colData中这一列含有分组信息
                       group1 = "normal", # normal组
                       group2 = "tumor" # tumor组
                       )

## Group1:normal
## Group2:tumor
## Calculating the p-values of each probe...
## |=================================================================|100%                       Completed after 20 m 
## Saving volcano plot as: methylation_volcano.pdf
```

查看结果:

```r
head(res)
##            mean.normal mean.tumor mean.normal.minus.mean.tumor
## cg21870274   0.7868787  0.5479701                   0.23890856
## cg16619049   0.2781233  0.2591982                   0.01892506
## cg18147296   0.7219909  0.6847735                   0.03721734
## cg13938959   0.6958773  0.4552416                   0.24063575
## cg12445832   0.4654454  0.2479300                   0.21751535
## cg23999112   0.7941180  0.5310289                   0.26308907
##            p.value.normal.tumor p.value.adj.normal.tumor
## cg21870274         9.732669e-18             2.668610e-16
## cg16619049         2.314279e-02             4.306785e-02
## cg18147296         1.024600e-01             1.566558e-01
## cg13938959         4.389698e-18             1.321596e-16
## cg12445832         2.608577e-18             8.391252e-17
## cg23999112         4.649043e-13             5.043637e-12
##                               status
## cg21870274 Hypermethylated in normal
## cg16619049           Not Significant
## cg18147296           Not Significant
## cg13938959 Hypermethylated in normal
## cg12445832 Hypermethylated in normal
## cg23999112 Hypermethylated in normal
```

可以看到和转录组的差异分析结果差不多，但是都是探针信息，基因信息需要注释才行。

同时也会在当前目录下生产一个PDF格式的火山图：

![image-20220903154948951](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903154948951.png)

保存结果：


```r
save(res, file = "tcga-coad_de_methy.rdata")
```


## 甲基化可视化

使用箱线图可视化不同组别之间的甲基化值。


```r
mdm <- TCGAvisualize_meanMethylation(data = beta.m,
                                     groupCol = "sample_type",
                                     print.pvalue = T
                                     )
```

会在目录下生成一个PDF文件：

![image-20220903154914628](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903154914628.png)

## 甲基化旭日图

可以在一张图中查看转录组数据和甲基化数据的情况。

需要准备转录组差异分析结果和甲基化差异分析结果。

转录组差异分析结果使用之前推文中得到的数据：[新版TCGAbiolinks包学习：差异分析](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247492794&idx=1&sn=ff9e4aa0c76b0d4adadf0d124f2d90a9&scene=21#wechat_redirect)


```r
# 这个数据只是部分差异基因，你可以用所有基因的结果
load(file = "coadDEGs.Rdata")
```

甲基化差异分析结果就用本次的结果。

有了两个差异分析的结果，就可以画旭日图了：


```r
starburst <- TCGAvisualize_starburst(
    met = res, 
    exp = coadDEGs,
    group1 = "normal",
    group2 = "tumor",
    met.platform = "Illumina Human Methylation 450",
    genome = "hg38",
    met.p.cut = 10^-5, 
    exp.p.cut = 10^-5,
    names = TRUE
)
```

会在当前目录下生成一个png格式的旭日图：

![starburst](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/starburst.png)







