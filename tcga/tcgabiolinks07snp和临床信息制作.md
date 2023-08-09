
今天给大家演示下如何用自己的数据完成`maftools`的分析，主要是**snp文件和临床信息的制作**，其实很简单，但是网络上的教程都说的不清楚。

这次我们直接用之前**TCGA-COAD和TCGA-READ合并后的数据演示**，合并教程请看前一篇推文：

[新版TCGA数据库不同癌种的组学数据合并](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247493493&idx=1&sn=7d805d747aeee3de70cf5384457b5c59&scene=21#wechat_redirect)

## 加载数据和R包

因为现在的TCGA数据库不能直接下载4种制作好的maf文件了，需要自己整理，如果你还不知道怎么整理，请看这篇内容：

[TCGA的maf突变文件不能下载了？直接用TCGAbiolinks包搞定！](https://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247493492&idx=1&sn=0e7a53ab9f60cdf9abcdf672ad36c834&scene=21#wechat_redirect)


```r
load(file = "./TCGA-colrectum/colrectal_snp.rdata")

library(maftools)
```

除此之外，我们还要给这个数据添加临床信息，也是只要先加载之前合并后的数据。


```r
load(file = "./TCGA-colrectum/colrectal_clin.rdata")
```

`maftools`包添加临床信息的方式非常简单，**只要两个文件都有`Tumor_Sample_Barcode`这一列且对得上就行**。所以我们还要对snp文件和临床信息进行一些简单的处理。

- 对于两个文件中的`Tumor_Sample_Barcode`这一列，我们只要前12个字符即可
- 临床信息中有一些是Normal的样本，需要去除
- 之选择在snp文件中有的样本


```r
# 只要前12个字符
colrec_snp$Tumor_Sample_Barcode <- substr(colrec_snp$Tumor_Sample_Barcode,1,12)
head(colrec_snp$Tumor_Sample_Barcode)
## [1] "TCGA-D5-6530" "TCGA-D5-6530" "TCGA-D5-6530" "TCGA-D5-6530" "TCGA-D5-6530"
## [6] "TCGA-D5-6530"
```

**如果你是像我这样直接用的`TCGAbiolinks`包下载的数据，那这个临床信息直接包含了`sample_type`这一列，不需要自己根据样本名确定到底是normal还是tumor，十分方便。**


```r
index <- unique(colrec_snp$Tumor_Sample_Barcode)

# 只要肿瘤样本
clin_snp <- clin[!clin$sample_type == "Solid Tissue Normal", ]

# 只要snp文件中有的样本
clin_snp <- clin_snp[clin_snp$patient %in% index, ]

# clin中没有Tumor_Sample_Barcode这一列，直接添加一列
clin_snp$Tumor_Sample_Barcode <- clin_snp$patient
```

这样两个需要的文件就制作好了。


```r
colrec_snp[1:5,1:5]
##   X1 Hugo_Symbol Entrez_Gene_Id Center NCBI_Build
## 1  1        AGRN         375790    BCM     GRCh38
## 2  1       ACAP3         116983    BCM     GRCh38
## 3  1      CALML6         163688    BCM     GRCh38
## 4  1       PRKCZ           5590    BCM     GRCh38
## 5  1      WRAP73          49856    BCM     GRCh38
```


```r
clin_snp[1:5,1:5]
##                                                   barcode      patient
## TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664
## TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530
## TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556
## TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818
## TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660
##                                        sample shortLetterCode
## TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664-01A              TP
## TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530-01A              TP
## TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556-01A              TP
## TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818-01A              TP
## TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660-01A              TP
##                                       definition
## TCGA-A6-5664-01A-21R-1839-07 Primary solid Tumor
## TCGA-D5-6530-01A-11R-1723-07 Primary solid Tumor
## TCGA-AA-3556-01A-01R-0821-07 Primary solid Tumor
## TCGA-AA-3818-01A-01R-0905-07 Primary solid Tumor
## TCGA-AA-3660-01A-01R-1723-07 Primary solid Tumor
```

## 读取数据

两个文件都没有问题，直接读取即可。


```r
colrec.maf <- read.maf(colrec_snp,clinicalData = clin_snp)
## -Validating
## --Removed 41322 duplicated variants
## -Silent variants: 67679 
## -Summarizing
## --Mutiple centers found
## BCM;WUGSC;BCM;WUGSC;BCM;BI--Possible FLAGS among top ten genes:
##   TTN
##   SYNE1
##   MUC16
## -Processing clinical data
## --Annotation missing for below samples in MAF:
##   TCGA-AA-3695
##   TCGA-AA-3967
##   TCGA-AF-2689
##   TCGA-AF-3914
##   TCGA-AG-3891
##   TCGA-AG-3906
##   TCGA-AZ-4681
## -Finished in 7.170s elapsed (6.480s cpu)
```

画图也是没有问题的，关于这个包的使用网络上教程非常多，也很全面，大家直接百度即可，这里就简单演示下。


```r
plotmafSummary(colrec.maf)
```

![plot of chunk unnamed-chunk-8](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-8-151561119.png)

最常见的额瀑布图：


```r
oncoplot(maf = colrec.maf,
         clinicalFeatures = c("ajcc_pathologic_stage","vital_status"),
         top = 30,
         sortByAnnotation=T
         )
```

![image-20220820183535812](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220820183535812.png)

这个图其实就是`Complexheatmap`画出来的，之前的推文里详细介绍了这个R包，感兴趣的可以看一看：

[超详细的R语言热图之complexheatmap系列1](https://mp.weixin.qq.com/s/sIGLjqk_Ug4FfwrzWrXprQ)

[超详细的R语言热图之complexheatmap系列2](https://mp.weixin.qq.com/s/3WA9hoHfktm7ZioGC0imEA)

[韦恩图进阶！complexheatmap包画upset plot](https://mp.weixin.qq.com/s/CI-wTadPj2lLuEGM_9tTYA)


```r
colrec.titv = titv(maf = colrec.maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = colrec.titv)
```

![plot of chunk unnamed-chunk-10](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-10-151561119.png)

比较有意思的棒棒糖图，这个图可以用`trackviewer`包画的更好看，下次介绍。

画出来简单，但是解释过程挺复杂的，涉及到下游氨基酸的变化，需要了解**碱基变化命名规则和氨基酸变化命名规则**才能知道具体意思。


```r
lollipopPlot(colrec.maf,
             gene = "TP53",
             AACol = "HGVSp_Short" # 需要注意这一列
             )
## 8 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.
##    HGNC    refseq.ID   protein.ID aa.length
## 1: TP53    NM_000546    NP_000537       393
## 2: TP53 NM_001126112 NP_001119584       393
## 3: TP53 NM_001126118 NP_001119590       354
## 4: TP53 NM_001126115 NP_001119587       261
## 5: TP53 NM_001126113 NP_001119585       346
## 6: TP53 NM_001126117 NP_001119589       214
## 7: TP53 NM_001126114 NP_001119586       341
## 8: TP53 NM_001126116 NP_001119588       209
## Using longer transcript NM_000546 for now.
```

![plot of chunk unnamed-chunk-11](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-151561119.png)

拷贝数变异肯定也是没有问题的，也是用的之前合并后的数据，然后经过gistic处理，就得到了我们需要的文件，关于gistic这个软件的使用，大家百度即可~


```r
all.lesions <- "./TCGA-colrectum/TCGA_COREAD_results/all_lesions.conf_90.txt"
amp.genes <- "./TCGA-colrectum/TCGA_COREAD_results/amp_genes.conf_90.txt"
del.genes <- "./TCGA-colrectum/TCGA_COREAD_results/del_genes.conf_90.txt"
scores.gis <- "./TCGA-colrectum/TCGA_COREAD_results/scores.gistic"

colrec.gistic = readGistic(gisticAllLesionsFile = all.lesions, 
                           gisticAmpGenesFile = amp.genes, 
                           gisticDelGenesFile = del.genes, 
                           gisticScoresFile = scores.gis, 
                           isTCGA = TRUE)
## -Processing Gistic files..
## --Processing amp_genes.conf_90.txt
## --Processing del_genes.conf_90.txt
## --Processing scores.gistic
## --Summarizing by samples

colrec.gistic
## An object of class  GISTIC 
##           ID summary
## 1:   Samples     611
## 2:    nGenes    2791
## 3: cytoBands      88
## 4:       Amp   97455
## 5:       Del  313297
## 6:     total  410752
```

结果也是毫无问题！

准备自己的数据，就是如此的简单。
