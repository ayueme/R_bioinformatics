
# 使用`tinyarray`包简化你的GEO和TCGA数据分析流程！

## 简介

做生信数据分析的小伙伴一定离不开TCGA和GEO这两个数据库，但是在下载以及整理的数据的过程中就会遇到很多麻烦，比如探针id的转换，各种id的转换，数据的过滤，基本的差异分析，PCA，聚类，批量生存分析等。

真正分析的代码也就那几句，但是整理数据的过程真的是太麻烦了！！

今天介绍的这个`tinyarray`包就可以帮助我们解决这些问题，大幅度简化你的数据清理过程！

这个包是生信技能树团队小洁老师写的R包哟，真的是太棒啦！！

## 下载安装

在线安装：


```r
if(!require(devtools))install.packages("devtools")
if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray",upgrade = F)
```

下载zip包后本地安装：


```r
devtools::install_local("tinyarray-master.zip",upgrade = F,dependencies = T)
```

## 简化GEO数据分析

### 简化GEO数据下载


```r
library(tinyarray)

# 直接用GSE号即可，默认会通过曾老师的GEO中国镜像下载，超级快，不需要fq
gse <- geo_download("GSE38713")


载入需要的程辑包：AnnoProbe
AnnoProbe v 0.1.0  welcome to use AnnoProbe!
If you use AnnoProbe in published research, please acknowledgements:
We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.
54675 probes, 43 samples from 0.923767449 to 15.77695556
```

如果报错可添加`by_annopbrobe = T`从官方途径下载。


```r
class(gse)
length(gse)

[1] "list"
[1] 3
```

下载来的`gse`是一个长度为3的列表，第1个就是表达矩阵，第2个是样本信息（pdata），第3个是GPL信息


```r
exp <- gse$exp
exp[1:4,1:4] # 标准的表达矩阵，行是基因，列是样本

        GSM948550 GSM948551 GSM948552 GSM948553
1007_s_at  9.735874  9.271060  9.885801  9.306324
1053_at    7.403046  7.751990  7.579300  5.986061
117_at     2.670557  2.677983  2.758533  2.720663
121_at     4.118984  4.211329  4.787708  4.141341
```


```r
View(gse$pd)
```


```r
gse$gpl

[1] "GPL570"
```


```r
group_list=c(rep('normal',13),rep('UC',30))
group_list <- factor(group_list,levels = c("normal","UC"))
table(group_list) 
## group_list
## normal     UC 
##     13     30
```



### 简化ID转换

简化GEO数据探针转换：


```r
ids <- AnnoProbe::idmap('GPL570') # 配合AnnoProbe
## Setting options('download.file.method.GEOquery'='auto')
## Setting options('GEOquery.inmemory.gpl'=FALSE)
exp1 <- trans_array(exp, ids)

exp1[1:4,1:4]

41937 of 54675 rownames matched
20188 rownames transformed after duplicate rows removed
       GSM948550 GSM948551 GSM948552 GSM948553
RFC2    7.403046  7.751990  7.579300  5.986061
HSPA6   2.670557  2.677983  2.758533  2.720663
PAX8    4.118984  4.211329  4.787708  4.141341
GUCA1A  2.265967  2.257482  2.270735  2.287381
```

连带着把探针重复的问题也一起解决了，是不是非常方便呢！

### 简化差异分析

一句代码，完成差异分析，火山图，热图，PCA，简直不能更棒！


```r
deg <- get_deg_all(
  exp = exp, # 表达矩阵（探针转换前的）
  group_list = group_list, # 分组信息
  ids = ids, # 探针和基因名的对应表
  logFC_cutoff = 1, # logFC
  scale_before = F, # 是否scale
  cluster_cols = T # 热图聚类
)

'select()' returned 1:many mapping between keys and columns
[1] "313 down genes,410 up genes"
```

`deg`也是一个列表，第1个是差异分析结果，第2个是上下调基因，第3个是3张图。


```r
head(deg[[1]])
```

![image-20220204115124398](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220204115124398.png)

```r
str(deg[[2]])

List of 1
 $ deg:List of 3
  ..$ up  :'data.frame':	410 obs. of  2 variables:
  .. ..$ upgenes : chr [1:410] "MANF" "PSME4" "ANKRD22" "SDF2L1" ...
  .. ..$ upprobes: chr [1:410] "202655_at" "212219_at" "238439_at" "218681_s_at" ...
  ..$ down:'data.frame':	313 obs. of  2 variables:
  .. ..$ downgenes : chr [1:313] "TRPM6" "CDKN2B" "CDKN2B-AS1" "LRRN2" ...
  .. ..$ downprobes: chr [1:313] "221102_s_at" "207530_s_at" "1559884_at" "205154_at" ...
  ..$ diff:'data.frame':	723 obs. of  2 variables:
  .. ..$ diffgenes : chr [1:723] "MANF" "PSME4" "ANKRD22" "SDF2L1" ...
  .. ..$ diffprobes: chr [1:723] "202655_at" "212219_at" "238439_at" "218681_s_at" ...
```



```r
deg[[3]]
```

![image-20220204115152573](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220204115152573.png)



这个包的功能远不止此，还可以：简化TCGA数据的ID转换、多组的差异分析、批量生存分析、快速探索表达矩阵等。更多精彩，欢迎到作者的[Github](https://github.com/xjsun1221/tinyarray)中学习!






## 简化TCGA数据分析

### 表达矩阵的整理

简化id转换，简化分组，一步到位，提取mRNA和lncRNA!


```r
# 以xena下载的BRCA为例，其他也一样可以！
rm(list = ls())

library(tinyarray) # 加载包
## 
## 

survinfo <- data.table::fread("E:/projects/lisy/files/TCGA-BRCA.survival.tsv",data.table = F)

expr <- data.table::fread("E:/projects/lisy/files/TCGA-BRCA.htseq_fpkm.tsv.gz",data.table = F)

rownames(expr) <- expr$Ensembl_ID
expr <- expr[,-1] # 变成标准的表达矩阵格式，行是基因，列是样本

expr[1:4,1:4]
##                    TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A1EU-11A
## ENSG00000242268.2        0.09170787      0.000000000       0.05789928
## ENSG00000270112.3        0.01957347      0.004700884       0.01630174
## ENSG00000167578.15       2.23589760      1.863334319       1.70475318
## ENSG00000273842.1        0.00000000      0.000000000       0.00000000
##                    TCGA-A8-A06X-01A
## ENSG00000242268.2          0.000000
## ENSG00000270112.3          0.000000
## ENSG00000167578.15         1.947481
## ENSG00000273842.1          0.000000
```

分别提取mRNA和lncRNA，并转换id为gene symbol，一步到位！


```r
expr_mrna <- trans_exp(expr,mrna_only = T) # 提取mrna，转换id为gene symbol
## 19712 of genes successfully mapping to mRNA,14805 of genes successfully mapping to lncRNA
expr_mrna[1:4,1:4]
##         TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A1EU-11A TCGA-A8-A06X-01A
## RAB4B          2.2358976       1.86333432       1.70475318        1.9474808
## C12orf5        2.3219445       4.22669935       1.97575523        2.8087572
## RNF44          3.6200560       3.54611744       3.39694310        4.7232701
## DNAH3          0.3370874       0.01601615       0.04145534        0.0023613
```



```r
expr_lncrna <- trans_exp(expr,lncrna_only = T) # 提取lncrna，转换id为gene symbol
## 19712 of genes successfully mapping to mRNA,14805 of genes successfully mapping to lncRNA
expr_lncrna[1:4,1:4]
##               TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A1EU-11A
## RP11-368I23.2       0.09170787      0.000000000       0.05789928
## RP11-742D12.2       0.01957347      0.004700884       0.01630174
## EHD4-AS1            0.08466075      0.000000000       0.46162417
## RP11-166P13.4       0.05394113      0.000000000       0.50697074
##               TCGA-A8-A06X-01A
## RP11-368I23.2       0.00000000
## RP11-742D12.2       0.00000000
## EHD4-AS1            0.08891242
## RP11-166P13.4       0.19295943
```

tcga样本分为癌和癌旁，也是一步到位！


```r
group_list <- make_tcga_group(expr_mrna)

table(group_list)
## group_list
## normal  tumor 
##    113   1104
```

下面接差异分析就非常简单了！就不在演示了。

### 批量生存分析的简化

生存分析需要将表达矩阵和临床信息合并，或者变成顺序一样才行，直接手撸代码也是可以的，但是`tinyarray`包把这个过程简化了！



```r
# 去除重复的tumor样本
expr_mrna_fi <- sam_filter(expr_mrna)
## filtered 13 samples.

# 匹配tcga的表达矩阵和临床信息
match_exp_cl(expr_mrna_fi, survinfo, "_PATIENT")
## match successfully.

names(cl_matched)[4:5] <- c("event","time") # 把生存状态和生存时间名字改一下

identical(rownames(cl_matched), colnames(exp_matched))
## [1] TRUE
```


可以看到现在表达矩阵和临床信息都是1181个样本，而且顺序是一致的！

下面就可以批量做生存分析了。

首先展示下最佳截点的计算方法，比较简单的就是根据基因表达量的中位数作为截断值，但是这样有时会导致p值不显著，所以可以使用其他方法。


```r
# 挑选10个基因做，一共有19712个基因
meta <- cl_matched
exprSet <- exp_matched[1:10,]


point_cut(exprSet,meta)
## $RAB4B
## [1] 1.868292
## 
## $C12orf5
## [1] 1.468261
## 
## $RNF44
## [1] 3.582367
## 
## $DNAH3
## [1] 0.007520641
## 
## $RPL23A
## [1] 7.546291
## 
## $ARL8B
## [1] 5.517558
## 
## $CALB2
## [1] 3.726976
## 
## $MFSD3
## [1] 2.525482
## 
## $PIGV
## [1] 3.45386
## 
## $ZNF708
## [1] 2.11762
```

可以看到每个基因的最佳截点都给你算好了！

批量KM生存分析：


```r
surv_KM(exprSet,meta,cut.point = T) # 使用最佳截点
##        MFSD3        CALB2        RAB4B       RPL23A        DNAH3       ZNF708 
## 2.985878e-11 1.549063e-06 8.240647e-05 6.176400e-03 1.078377e-02 1.764367e-02 
##      C12orf5        RNF44 
## 4.037767e-02 4.395910e-02
```

批量cox生存分析：


```r
surv_cox(exprSet, meta, cut.point = T)
##               coef        se         z            p        HR       HRse
## RAB4B   -0.5685845 0.1460777 -3.892343 9.928061e-05 0.5663265 0.08272767
## C12orf5 -0.4358979 0.2147148 -2.030125 4.234383e-02 0.6466837 0.13885258
## RNF44   -0.3019320 0.1504873 -2.006362 4.481766e-02 0.7393884 0.11126856
## DNAH3    0.6244150 0.2484815  2.512924 1.197353e-02 1.8671534 0.46395306
## RPL23A  -0.4741161 0.1745146 -2.716770 6.592238e-03 0.6224350 0.10862400
## CALB2    0.7304669 0.1551399  4.708441 2.496186e-06 2.0760497 0.32207808
## MFSD3   -1.0713457 0.1684207 -6.361128 2.002771e-10 0.3425472 0.05769205
## ZNF708   0.4487043 0.1901424  2.359833 1.828314e-02 1.5662815 0.29781645
##                HRz          HRp    HRCILL    HRCIUL
## RAB4B    -5.242182 1.586884e-07 0.4253293 0.7540644
## C12orf5  -2.544542 1.094210e-02 0.4245476 0.9850483
## RNF44    -2.342186 1.917117e-02 0.5505257 0.9930420
## DNAH3     1.869054 6.161529e-02 1.1472872 3.0387000
## RPL23A   -3.475889 5.091623e-04 0.4421268 0.8762764
## CALB2     3.340959 8.348949e-04 1.5317308 2.8137988
## MFSD3   -11.395899 0.000000e+00 0.2462411 0.4765192
## ZNF708    1.901444 5.724382e-02 1.0789972 2.2736273
```

批量画箱线图：


```r
boxes <- exp_boxplot(exprSet, color = c("blue","red"))

boxes[[1]]
```

![unnamed-chunk-23-134515972](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-23-134515972.png)

拼图：


```r
patchwork::wrap_plots(boxes, nrow = 2)
```

![unnamed-chunk-24-134515972](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-24-134515972.png)

批量画生存分析图：


```r
surv_plots <- exp_surv(exprSet, meta)

patchwork::wrap_plots(surv_plots, nrow = 3)
```

![image-20220204114708988](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220204114708988.png)


这个包的功能远不止此，还可以：多组的差异分析、快速探索表达矩阵、一句代码画热图等。更多精彩，欢迎到作者的[Github](https://github.com/xjsun1221/tinyarray)中学习!





