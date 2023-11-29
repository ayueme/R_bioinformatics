>本文主要参考官方介绍：https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html

[toc]

## 简介

分子分型一直是生信数据挖掘的热门技能，用于分子分型的算法非常多，比如大家常见的非负矩阵分解、一致性聚类、PCA等，一致性聚类我们在之前也介绍过了：[免疫浸润结果分子分型](https://mp.weixin.qq.com/s/96s_hfBH0HjLvvTfNgTIlQ)

今天给大家介绍一个一站式的分子分型R包：`MOVICS`。

该包与其他分子分型R包最大的不同是它能**同时使用多组学的数据**，普通的分子分型R包只能通过一种组学数据进行分析，比如只能通过mRNA的表达矩阵进行分析。但是这R包它可以同时通过比如说mRNA、lncRNA、甲基化数据、突变数据进行分型。

之外，它还提供了分型之后每个亚型的探索以及每个亚型内的分析。所以说这是一个一站式的包。这个的功能主要分为三个部分，示意图如下：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230924081939905.png)

第一个部分是根据不同的组学数据进行分型。大部分是比较不同的分型。第三个部分是对每个分型进行探索，以及获得每个分型特异性的分子。

每个部分包含的主要函数如下，下面会介绍：

- GET Module: get subtypes through multi-omics integrative clustering
  - getElites(): get elites which are those features that pass the filtering procedure and are used for analyses
  - getClustNum(): get optimal cluster number by calculating clustering prediction index (CPI) and Gap-statistics
  - getalgorithm_name(): get results from one specific multi-omics integrative clustering algorithm with detailed parameters
  - getMOIC(): get a list of results from multiple multi-omics integrative clustering algorithm with parameters by default
  - getConsensusMOIC(): get a consensus matrix that indicates the clustering robustness across different clustering algorithms and generate a consensus heatmap
  - getSilhouette(): get quantification of sample similarity using silhoutte score approach
  - getStdiz(): get a standardized data for generating comprehensive multi-omics heatmap
  - getMoHeatmap(): get a comprehensive multi-omics heatmap based on clustering results

- COMP Module: compare subtypes from multiple perspectives
  - compSurv(): compare survival outcome and generate a Kalan-Meier curve with pairwise comparison if possible
  - compClinvar(): compare and summarize clinical features among different identified subtypes
  - compMut(): compare mutational frequency and generate an OncoPrint with significant mutations
  - compTMB(): compare total mutation burden among subtypes and generate distribution of Transitions and Transversions
  - compFGA(): compare fraction genome altered among subtypes and generate a barplot for distribution comparison
  - compDrugsen(): compare estimated half maximal inhibitory concentration (IC50
  ) for drug sensitivity and generate a boxviolin for distribution comparison
  - compAgree(): compare agreement of current subtypes with other pre-existed classifications and generate an alluvial diagram and an agreement barplot

- RUN Module: run marker identification and verify subtypes
  - runDEA(): run differential expression analysis with three popular methods for choosing, including edgeR, DESeq2, and limma
  - runMarker(): run biomarker identification to determine uniquely and significantly differential expressed genes for each subtype
  - runGSEA(): run gene set enrichment analysis (GSEA), calculate activity of functional pathways and generate a pathway-specific heatmap
  - runGSVA(): run gene set variation analysis to calculate enrichment score of each sample based on given gene set list of interest
  - runNTP(): run nearest template prediction based on identified biomarkers to evaluate subtypes in external cohorts
  - runPAM(): run partition around medoids classifier based on discovery cohort to predict subtypes in external cohorts
  - runKappa(): run consistency evaluation using Kappa statistics between two appraisements that identify or predict current subtypes

该包已发表，使用时记得引用：

- Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020). MOVICS: an R package for multi-omics integration and visualization in cancer subtyping. bioRxiv, 2020.2009.2015.297820. [doi.org/10.1101/2020.09.15.297820]

## 安装

目前该包在github，只能通过以下方式安装，注意安装时最好先安装依赖包，因为这个包的依赖包非常多，安装过程中非常容易失败。对于初学者来说，这个包的安装不是很友好哦~


```r
# 网络安装
devtools::install_github("xlucpu/MOVICS")

# 或者下载到本地安装
devtools::install_local("E:/R/R包/MOVICS-master.zip")
```

## GET Module

### 准备数据

我们先看一下示例数据。


```r
library(MOVICS)
## 
```

使用该包自带数据进行演示，这个自带数据是已经清洗好的。过几天再专门写一篇推文介绍怎么准备这个数据。


```r
# TCGA的乳腺癌数据
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))
```

`brca.tcga`里面是多个组学的数据，比如mRNA、lncRNA、甲基化、突变数据等，还有临床信息，比如生存时间和生存状态以及乳腺癌的PAM50分类。

为了演示，这个数据通过MAD筛选了部分数据：
- 500 mRNAs, 
- 500 lncRNA, 
- 1,000 promoter CGI probes/genes with high variation 
- 30 genes that mutated in at least 3% of the entire cohort. 

注意，这里最重要的一点是：每种组学的数据的**样本数量、名字、顺序应该完全一致**。大家可以自己看一下这些数据是什么样的。


```r
names(brca.tcga)
## [1] "mRNA.expr"   "lncRNA.expr" "meth.beta"   "mut.status"  "count"      
## [6] "fpkm"        "maf"         "segment"     "clin.info"

names(brca.yau)
## [1] "mRNA.expr" "clin.info"


# 提取"mRNA.expr""lncRNA.expr""meth.beta""mut.status"
mo.data   <- brca.tcga[1:4]

# 提取raw count data
count     <- brca.tcga$count

# 提取fpkm data
fpkm      <- brca.tcga$fpkm

# 提取maf
maf       <- brca.tcga$maf

# 提取segmented copy number
segment   <- brca.tcga$segment

# 提取生存信息
surv.info <- brca.tcga$clin.info
```

### 筛选基因（降维）

`getElites`，顾名思义，找出精英，找出最牛逼的，也就是说这个函数可以做一些预处理和筛选工作，可以帮你进行数据准备工作。

主要可以做以下预处理：

- 缺失值插补：直接删除或者knn插补
- 筛选分子：可根据mad, sd, pca, cox, freq(二分类数据)进行筛选

其实这个不是第一步，第一步应该是自己先清洗一下数据，比如表达矩阵先进行log转换等。

下面是一些功能演示，还是非常强大的。

缺失值插补：


```r
# scenario 1: 处理缺失值
tmp       <- brca.tcga$mRNA.expr # get expression data
dim(tmp) # check data dimension
## [1] 500 643

tmp[1,1]  <- tmp[2,2] <- NA # 添加几个NA
tmp[1:3,1:3] # check data
##         BRCA-A03L-01A BRCA-A04R-01A BRCA-A075-01A
## SCGB2A2            NA          1.42          7.24
## SCGB1D2         10.11            NA          5.88
## PIP              4.54          2.59          4.35

elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       na.action = "rm", # 直接删除
                       elite.pct = 1) # 保留100%的数据
## --2 features with NA values are removed.
## missing elite.num then use elite.pct

dim(elite.tmp$elite.dat) 
## [1] 498 643


elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       na.action = "impute", # 使用knn进行插补
                       elite.pct = 1) 
## missing elite.num then use elite.pct

dim(elite.tmp$elite.dat) 
## [1] 500 643

elite.tmp$elite.dat[1:3,1:3] # NA values have been imputed 
##         BRCA-A03L-01A BRCA-A04R-01A BRCA-A075-01A
## SCGB2A2         6.867         1.420          7.24
## SCGB1D2        10.110         4.739          5.88
## PIP             4.540         2.590          4.35
```

使用MAD筛选分子：


```r
# scenario 2: 使用MAD筛选，最大中位差
tmp       <- brca.tcga$mRNA.expr 
elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       elite.pct = 0.1) # 保留MAD前10%的基因
## missing elite.num then use elite.pct

dim(elite.tmp$elite.dat) # 500的10%是50
## [1]  50 643
#> [1]  50 643

elite.tmp <- getElites(dat       = tmp,
                       method    = "sd",
                       elite.num = 100, # 保留MAD前100的基因
                       elite.pct = 0.1) # 此时这个参数就不起作用了
## elite.num has been provided then discards elite.pct.

dim(elite.tmp$elite.dat) 
## [1] 100 643
```

使用PCA筛选分子，需要了解一些关于PCA的基础知识：[R语言主成分分析](https://mp.weixin.qq.com/s/B_RD4pAemEE7pHkr2VfDwA)


```r
# scenario 3: 使用PCA筛选分子
tmp       <- brca.tcga$mRNA.expr # get expression data with 500 features
elite.tmp <- getElites(dat       = tmp,
                       method    = "pca",
                       pca.ratio = 0.95) # 主成分的比例
## --the ratio used to select principal component is set as 0.95
dim(elite.tmp$elite.dat) # get 204 elite (PCs) left
## [1] 204 643
```

使用单因素COX回归筛选分子，也就是对每个分子做单因素cox分析，选择有意义的留下，需要提供生存信息：


```r
# scenario 4: 使用cox筛选分子
tmp       <- brca.tcga$mRNA.expr # get expression data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # 生存信息，列名必须有'futime'和'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # 此时这个参数也是不起作用的
## --all sample matched between omics matrix and survival data.
## 5% 10% 15% 20% 25% 30% 35% 40% 45% 50% 55% 60% 65% 70% 75% 80% 85% 90% 95% 100%

dim(elite.tmp$elite.dat) # get 125 elites
## [1] 125 643

table(elite.tmp$unicox$pvalue < 0.05) # 125 genes have nominal pvalue < 0.05 in 
## 
## FALSE  TRUE 
##   375   125

tmp       <- brca.tcga$mut.status # get mutation data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, 
                       p.cutoff  = 0.05,
                       elite.num = 100) 
## --all sample matched between omics matrix and survival data.
## 7% 13% 20% 27% 33% 40% 47% 53% 60% 67% 73% 80% 87% 93% 100%
dim(elite.tmp$elite.dat) # get 3 elites
## [1]   3 643

table(elite.tmp$unicox$pvalue < 0.05) # 3 mutations have nominal pvalue < 0.05
## 
## FALSE  TRUE 
##    27     3
```

使用突变频率筛选分子，这个是准们用于0/1矩阵这种二分类数据的：


```r
# scenario 5: 使用突变频率筛选
tmp       <- brca.tcga$mut.status # get mutation data 
rowSums(tmp) 
## PIK3CA   TP53    TTN   CDH1  GATA3   MLL3  MUC16 MAP3K1  SYNE1  MUC12    DMD 
##    208    186    111     83     58     49     48     38     33     32     31 
##  NCOR1    FLG   PTEN   RYR2  USH2A  SPTA1 MAP2K4  MUC5B    NEB   SPEN  MACF1 
##     31     30     29     27     27     25     25     24     24     23     23 
##   RYR3    DST  HUWE1  HMCN1  CSMD1  OBSCN   APOB  SYNE2 
##     23     22     22     22     21     21     21     21
elite.tmp <- getElites(dat       = tmp,
                       method    = "freq", # must set as 'freq'
                       elite.num = 80, # 这里是指突变频率
                       elite.pct = 0.1) # 此时该参数不起作用
## --method of 'freq' only supports binary omics data (e.g., somatic mutation matrix), and in this manner, elite.pct and elite.num are used to cut frequency.
## elite.num has been provided then discards elite.pct.

rowSums(elite.tmp$elite.dat) # 只保留在80个及以上样本中突变的基因
## PIK3CA   TP53    TTN   CDH1 
##    208    186    111     83

elite.tmp <- getElites(dat       = tmp,
                       method    = "freq", 
                       elite.pct = 0.2) 
## --method of 'freq' only supports binary omics data (e.g., somatic mutation matrix), and in this manner, elite.pct and elite.num are used to cut frequency.
## missing elite.num then use elite.pct
rowSums(elite.tmp$elite.dat) # only genes that are mutated in over than 0.2*643=128.6 
## PIK3CA   TP53 
##    208    186
```


### 确定最佳亚型数量

根据分子表达量对样本进行分型，分子就是上一步得到的mRNA、lncRNA、miRNA、甲基化矩阵等。

先根据CPI和Gaps-statistics确定分成几个亚型：


```r
optk.brca <- getClustNum(data        = mo.data, # 4种组学数据
                         is.binary   = c(F,F,F,T), #前3个不是二分类的，最后一个是
                         try.N.clust = 2:8, # 尝试亚型数量，从2到8
                         fig.name    = "CLUSTER NUMBER OF TCGA-BRCA")#保存的文件名
## calculating Cluster Prediction Index...
## 5% complete
## 5% complete
## 10% complete
## 10% complete
## 15% complete
## 15% complete
## 20% complete
## 25% complete
## 25% complete
## 30% complete
## 30% complete
## 35% complete
## 35% complete
## 40% complete
## 45% complete
## 45% complete
## 50% complete
## 50% complete
## 55% complete
## 55% complete
## 60% complete
## 65% complete
## 65% complete
## 70% complete
## 70% complete
## 75% complete
## 75% complete
## 80% complete
## 85% complete
## 85% complete
## 90% complete
## 90% complete
## 95% complete
## 95% complete
## 100% complete
## calculating Gap-statistics...
## visualization done...
## --the imputed optimal cluster number is 3 arbitrarily, but it would be better referring to other priori knowledge.
```

![unnamed-chunk-10-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-10-186542957.png)

会自动在当前工作目录下产生一个PDF格式的图片。

函数给出的结果是3，但是考虑到乳腺癌的PAM0分类，我们选择k=5，也就是分成5个亚型。

所以这个确定最佳亚型个数是根据你自己的需要来的哈，灵活调整~

### 根据单一算法分型

确定分成几个亚型之后，可以通过算法进行分型了。提供了非常多的方法，大家常见的非负矩阵分解、异质性聚类等等都提供了。

比如根据贝叶斯方法进行分型：


```r
# perform iClusterBayes (may take a while)
iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                      N.clust     = 5,
                                      type        = c("gaussian","gaussian","gaussian","binomial"),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                                      sdev        = 0.05,
                                      thin        = 3)
## clustering done...
## feature selection done...
```

或者使用统一的函数，自己选择方法即可，两种方法得到的结果完全是一样的：


```r
iClusterBayes.res <- getMOIC(data        = mo.data,
                             N.clust     = 5,
                             methodslist = "iClusterBayes", # 指定算法
                             type        = c("gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
                             n.burnin    = 1800,
                             n.draw      = 1200,
                             prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                             sdev        = 0.05,
                             thin        = 3)
```

返回的结果包含一个`clust.res`对象，它有两列：`clust`列指示样本所属的亚型，`samID`列记录对应的样本名称。对于提供特征选择过程的算法（如iClusterBayes、CIMLR和MoCluster），结果还包含一个`feat.res`对象，存储了这种过程的信息。对于涉及分层聚类的算法（例如COCA、ConsensusClustering），样本聚类的相应树状图也将作为`clust.dend`返回，如果用户想要将它们放在热图中会很有用。

### 同时进行多种分型算法

可以同时根据多种算法进行分型，然后整合它们的结果，得到最终的结果，不是一般的强大：


```r
# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"), # 9种算法
                         N.clust     = 5,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))
## --you choose more than 1 algorithm and all of them shall be run with parameters by default.
## SNF done...
## Clustering method: kmeans
## Perturbation method: noise
## PINSPlus done...
## NEMO done...
## COCA done...
## LRAcluster done...
## end fraction
## clustered
```

```R
## ConsensusClustering done...
## IntNMF done...
## clustering done...
## feature selection done...
## CIMLR done...
## clustering done...
## feature selection done...
## MoCluster done...
```

再把贝叶斯的结果一起加进来，这就是10种算法了：


```r
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# 保存下结果
save(moic.res.list, file = "moic.res.list.rda")
```

### 整合多种分型结果

借鉴了consensus ensembles的想法，实现对多个分型算法结果的整合。

可以画出一个一致性热图：


```r
load(file = "moic.res.list.rda")
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
```

![unnamed-chunk-15-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-15-186542957.png)

结果会保存在当前工作目录中。

### 查看分型结果的质量

除了通过上面的热图查看分型结果，还可以使用Silhouette准则判断分型质量。

以下是解释，来源于网络：

>Silhouette准则是一种用于聚类分析中的评价方法，它通过对每个数据点与其所属簇内其他数据点之间的距离进行比较，来衡量聚类质量的好坏。Silhouette准则可以帮助我们确定最佳的聚类数量，从而提高聚类分析的可靠性和准确性。
>Silhouette准则的计算方法如下：对于每个数据点i，计算它与同簇中其他数据点之间的平均距离ai，以及与最近其他簇中数据点之间的平均距离bi。然后，定义每个数据点的Silhouette系数为：
s(i) = (bi - ai) / max(ai, bi)
>Silhouette系数的取值范围在-1到1之间，其中负值表示数据点更容易被分类到错误的簇中，而正值则表示数据点更容易被正确分类。Silhouette系数的平均值可以用来评估整个聚类的质量，因此，Silhouette准则的目标是最大化Silhouette系数的平均值，从而找到最佳的聚类数量。
>当聚类数量增加时，Silhouette系数的平均值通常会先增加后减少。因此，我们需要找到一个聚类数量，使得Silhouette系数的平均值达到最大值。通常，我们会通过绘制Silhouette图来选择最佳的聚类数量。Silhouette图是一种以Silhouette系数为纵轴，聚类数量为横轴的图表，它可以帮助我们直观地理解聚类的质量。
>在使用Silhouette准则进行聚类分析时，需要注意以下几点：
>1. Silhouette系数只适用于欧氏距离或相关度量，对于其他距离度量可能不适用。
>2. Silhouette系数的计算时间较长，因此在处理大规模数据时需要注意计算效率。
>3. Silhouette系数并不是唯一的评价指标，对于特定的聚类问题可能需要采用其他评价指标。

结果会保存在当前工作目录中：


```r
getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)
```

![unnamed-chunk-16-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-16-186542957.png)

```
## png 
##   2
```

### 多组学分型热图

分型之后，肯定是要对每个组学数据进行热图展示不同亚型的表达量情况。

不过需要做一些准备工作。

- 把甲基化的β值矩阵转换为M值矩阵，作者推荐，这样做展示效果更好；
- 数据标准化，画热图之钱一般都会进行这个操作，其实是通过`scale`进行的，比如把所有数据压缩为[-2,2]，超过2的用2表示，小于-2的用-2表示


```r
# β值矩阵转换为M值矩阵
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# 对数据进行标准化
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation
```

我们这里就用贝叶斯分型的结果进行展示，首先是提取每个组学的结果，然后每个组学中选择前10个分子进行标注：


```r
feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)
```

下面就是画图即可，其实也是借助`complexheatmap`实现的，只不过帮你简化了很多过程，结果会自动保存在当前工作目录下，`MOVICS`的默认出图还是很美观的，可能比你自己画的好看~


```r
# 为每个组学的热图自定义颜色，不定义也可
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")
```

![unnamed-chunk-19-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-19-186542957.png)

上面是贝叶斯方法分型结果的展示，你也可以任选一种，毕竟我们有10种算法。

比如选择COCA法的结果进行展示，也是一模一样的用法，结果会自动保存：


```r
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")
```

![unnamed-chunk-20-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-20-186542957.png)

如果你要展示多个临床信息，也是直接添加即可，注意自定义颜色需要使用`circlize`实现：


```r
# extract PAM50, pathologic stage and age for sample annotation
annCol    <- surv.info[,c("PAM50", "pstage", "age"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  PAM50  = c("Basal" = "blue",
                            "Her2"   = "red",
                            "LumA"   = "yellow",
                            "LumB"   = "green",
                            "Normal" = "black"),
                  pstage = c("T1"    = "green",
                             "T2"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                             "TX"    = "black"))

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.brca$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             annRow        = NULL, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")
```

![unnamed-chunk-21-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-21-186542957.png)


是不是非常牛逼？

到这里第一部分的内容就介绍完了，下面就是探索、比较不同的亚型了。

## COMP Module

这部分就是探索不同亚型了。

### 生存分析

既然分成了5个亚型，肯定是要看看5个亚型的生存情况有没有区别，生存分析走一波：


```r
# survival comparison
surv.brca <- compSurv(moic.res         = cmoic.brca,
                      surv.info        = surv.info,
                      convt.time       = "m", # 把天变成月
                      surv.median.line = "h", 
                      xyrs.est         = c(5,10), # 计算5年和10年生存率
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
## --a total of 643 samples are identified.
## --removed missing values.
## --leaving 642 observations.
```

![unnamed-chunk-22-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-22-186542957.png)


查看结果：


```r
print(surv.brca)
## $fitd
## Call:
## survdiff(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
##     na.action = na.exclude)
## 
##               N Observed Expected (O-E)^2/E (O-E)^2/V
## Subtype=CS1 146       28     14.6    12.269    15.455
## Subtype=CS2 132       17     15.1     0.251     0.317
## Subtype=CS3 107        7     13.4     3.028     3.756
## Subtype=CS4 144        7     17.4     6.218     8.175
## Subtype=CS5 113       16     14.6     0.140     0.176
## 
##  Chisq= 22.3  on 4 degrees of freedom, p= 2e-04 
## 
## $fit
## Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
##     na.action = na.exclude, error = "greenwood", type = "kaplan-meier", 
##     conf.type = "plain")
## 
##       n events median 0.95LCL 0.95UCL
## CS1 146     28    130    67.3      NA
## CS2 132     17    114    83.6     144
## CS3 107      7     NA   102.5      NA
## CS4 144      7    216   113.5      NA
## CS5 113     16     NA    97.2      NA
## 
## $xyrs.est
## Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res)
## 
##                 Subtype=CS1 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##  1825     23      23    0.658  0.0647        0.543        0.798
##  3650      4       4    0.509  0.0830        0.370        0.701
## 
##                 Subtype=CS2 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##  1825     27       7    0.855  0.0541        0.755        0.968
##  3650      3       8    0.421  0.1286        0.231        0.766
## 
##                 Subtype=CS3 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##  1825     19       4    0.851  0.0701        0.724        1.000
##  3650      5       3    0.644  0.1183        0.449        0.923
## 
##                 Subtype=CS4 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##  1825     25       4    0.917  0.0443        0.834            1
##  3650      3       2    0.550  0.2027        0.267            1
## 
##                 Subtype=CS5 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##  1825     25      12    0.821  0.0502        0.728        0.926
##  3650      3       4    0.503  0.1599        0.269        0.938
## 
## 
## $overall.p
## [1] 0.000172211
## 
## $pairwise.p
## 
## 	Pairwise comparisons using Log-Rank test 
## 
## data:  mosurv.res and Subtype 
## 
##     CS1     CS2     CS3     CS4    
## CS2 0.11909 -       -       -      
## CS3 0.01345 0.11909 -       -      
## CS4 0.00055 0.03497 0.72167 -      
## CS5 0.16073 0.87433 0.16073 0.03022
## 
## P value adjustment method: BH
```


### 比较临床特征

查看每个亚型的临床特征，类似于基线资料表。


```r
clin.brca <- compClinvar(moic.res      = cmoic.brca,
                         var2comp      = surv.info, #需要比较的临床信息，行名须是样本名
                         strata        = "Subtype", # 分层变量，这里肯定是分型了
                         factorVars    = c("PAM50","pstage","fustat"), #分类变量名字
                         nonnormalVars = "futime", #非正态的连续性变量
                         exactVars     = "pstage", #需要使用精确概率法的变量
                         doWord        = TRUE, # 自动生成Word文档
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
## --all samples matched.
## Registered S3 methods overwritten by 'proxy':
##   method               from    
##   print.registry_field registry
##   print.registry_entry registry
```

查看结果：


```r
print(clin.brca$compTab)
##                           level                      CS1
## 1                      n                             146
## 2             fustat (%)      0              118 (80.8) 
## 3                             1               28 (19.2) 
## 4  futime (median [IQR])        719.00 [431.50, 1355.00]
## 5              PAM50 (%)  Basal                2 ( 1.4) 
## 6                          Her2               38 (26.0) 
## 7                          LumA               47 (32.2) 
## 8                          LumB               53 (36.3) 
## 9                        Normal                6 ( 4.1) 
## 10            pstage (%)     T1               40 (27.4) 
## 11                           T2               85 (58.2) 
## 12                           T3               12 ( 8.2) 
## 13                           T4                9 ( 6.2) 
## 14                           TX                0 ( 0.0) 
## 15                  age                    58.47 ± 13.75
##                         CS2                      CS3                      CS4
## 1                       133                      107                      144
## 2               115 (86.5)               100 (93.5)               137 (95.1) 
## 3                18 (13.5)                 7 ( 6.5)                 7 ( 4.9) 
## 4  746.50 [432.25, 1420.25] 703.00 [436.00, 1270.50] 848.50 [470.50, 1547.75]
## 5                 0 ( 0.0)                 0 ( 0.0)                 0 ( 0.0) 
## 6                 0 ( 0.0)                 0 ( 0.0)                 0 ( 0.0) 
## 7                78 (58.6)                91 (85.0)               128 (88.9) 
## 8                55 (41.4)                16 (15.0)                 0 ( 0.0) 
## 9                 0 ( 0.0)                 0 ( 0.0)                16 (11.1) 
## 10               35 (26.3)                32 (29.9)                40 (27.8) 
## 11               77 (57.9)                63 (58.9)                72 (50.0) 
## 12               17 (12.8)                10 ( 9.3)                31 (21.5) 
## 13                4 ( 3.0)                 2 ( 1.9)                 1 ( 0.7) 
## 14                0 ( 0.0)                 0 ( 0.0)                 0 ( 0.0) 
## 15            59.32 ± 14.03            59.62 ± 12.40            57.58 ± 12.23
##                         CS5      p    test
## 1                       113               
## 2                97 (85.8)   0.001        
## 3                16 (14.2)                
## 4  742.00 [470.00, 1605.00]  0.691 nonnorm
## 5               109 (96.5)  <0.001        
## 6                 0 ( 0.0)                
## 7                 0 ( 0.0)                
## 8                 0 ( 0.0)                
## 9                 4 ( 3.5)                
## 10               25 (22.1)   0.036   exact
## 11               70 (61.9)                
## 12               14 (12.4)                
## 13                3 ( 2.7)                
## 14                1 ( 0.9)                
## 15            55.67 ± 12.11  0.116
```


### 突变全景图

比较突变频率，也就是突变全景图，也是基于`complexheatmap`画出来的，类似于`maftools`：


```r
# mutational frequency comparison
mut.brca <- compMut(moic.res     = cmoic.brca,
                    mut.matrix   = brca.tcga$mut.status, # 0/1矩阵
                    doWord       = TRUE, # 生成Word文档
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # 保留在至少5%的样本中突变的基因
                    p.adj.cutoff = 0.05, # 保留padj<0.05的基因
                    innerclust   = TRUE, # 在每个亚型中进行聚类
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 2,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
## --all samples matched.
```

![unnamed-chunk-26-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-26-186542957.png)

```r
print(mut.brca)
##   Gene (Mutated)       TMB        CS1        CS2        CS3        CS4
## 1         PIK3CA 208 (32%) 55 (37.7%) 10 ( 7.5%) 77 (72.0%) 64 (44.4%)
## 2           TP53 186 (29%) 83 (56.8%) 12 ( 9.0%)  0 ( 0.0%)  8 ( 5.6%)
## 3            TTN 111 (17%) 41 (28.1%) 14 (10.5%) 14 (13.1%) 22 (15.3%)
## 4           CDH1  83 (13%) 11 ( 7.5%)  3 ( 2.3%) 19 (17.8%) 49 (34.0%)
## 5          GATA3  58 ( 9%)  5 ( 3.4%) 32 (24.1%) 10 ( 9.3%) 11 ( 7.6%)
## 6           MLL3  49 ( 8%)  11 (7.5%)  12 (9.0%)  10 (9.3%)  11 (7.6%)
## 7          MUC16  48 ( 8%) 18 (12.3%)  8 ( 6.0%)  6 ( 5.6%)  8 ( 5.6%)
## 8         MAP3K1  38 ( 6%)  2 ( 1.4%)  6 ( 4.5%) 16 (15.0%) 12 ( 8.3%)
## 9          SYNE1  33 ( 5%)   9 (6.2%)   1 (0.8%)   5 (4.7%)   8 (5.6%)
##          CS5   pvalue     padj
## 1  2 ( 1.8%) 6.90e-42 3.10e-41
## 2 83 (73.5%) 2.06e-63 1.85e-62
## 3 20 (17.7%) 2.03e-03 3.04e-03
## 4  1 ( 0.9%) 5.38e-19 1.61e-18
## 5  0 ( 0.0%) 6.21e-11 1.40e-10
## 6   5 (4.4%) 6.31e-01 6.31e-01
## 7  8 ( 7.1%) 2.05e-01 2.31e-01
## 8  2 ( 1.8%) 3.65e-05 6.57e-05
## 9  10 (8.8%) 3.27e-02 4.20e-02
```


### 比较突变负荷

>毋庸置疑，免疫治疗正在成为现代癌症治疗的支柱。最近的分析将肿瘤基因组景观与抗肿瘤免疫联系起来。特别是，新兴的研究显示，肿瘤特异性基因组损伤与免疫检查点激活以及患者对免疫治疗的反应程度和持续时间有关。这些损伤包括高突变负荷和非整倍体情况。为了定量这些可能影响免疫治疗的基因组改变，MOVICS提供了两个函数来计算总突变负荷（TMB）和基因组改变比例（FGA）。具体而言，TMB指的是在肿瘤基因组中发现的突变数量，而FGA是受复制数增减影响的基因组百分比。这两个属性对遗传学研究人员非常有用，因为它们为他们提供了更深入的关于肿瘤基因组构成的信息。我们从`compTMB()`开始。首先，用于此函数的输入`maf`数据必须至少具有以下10列：


```r
names(maf)
##  [1] "Tumor_Sample_Barcode"   "Hugo_Symbol"            "Chromosome"            
##  [4] "Start_Position"         "End_Position"           "Variant_Classification"
##  [7] "Variant_Type"           "Reference_Allele"       "Tumor_Seq_Allele1"     
## [10] "Tumor_Seq_Allele2"
head(maf)
##   Tumor_Sample_Barcode Hugo_Symbol Chromosome Start_Position End_Position
## 1        BRCA-A1XY-01A       USP24       chr1       55159655     55159655
## 2        BRCA-A1XY-01A      ERICH3       chr1       74571494     74571494
## 3        BRCA-A1XY-01A      KIF26B       chr1      245419680    245419680
## 4        BRCA-A1XY-01A       USP34       chr2       61189055     61189055
## 5        BRCA-A1XY-01A      ANTXR1       chr2       69245305     69245305
## 6        BRCA-A1XY-01A       SCN9A       chr2      166199365    166199365
##   Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1
## 1      Missense_Mutation          SNP                T                 T
## 2      Missense_Mutation          SNP                C                 C
## 3                 Silent          SNP                G                 G
## 4                 Silent          SNP                G                 G
## 5                 Silent          SNP                G                 G
## 6                 Silent          SNP                G                 G
##   Tumor_Seq_Allele2
## 1                 C
## 2                 T
## 3                 T
## 4                 C
## 5                 A
## 6                 A
```

关于肿瘤突变负荷，我之前也介绍过：[1行代码计算肿瘤突变负荷TMB](https://mp.weixin.qq.com/s/TPURe613FXKi1tMHzAcJFA)

这个函数可以一次计算不同亚型的TMB，同时展示碱基突变情况：


```r
# compare TMB
tmb.brca <- compTMB(moic.res     = cmoic.brca,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")
## --67 samples mismatched from current subtypes.
## -Validating
## -Silent variants: 24329 
## -Summarizing
## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
## -Processing clinical data
## --Missing clinical data
## -Finished in 3.230s elapsed (2.030s cpu) 
## Kruskal-Wallis rank sum test p value = 1.97e-23
## post-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:
##     CS1        CS2        CS3        CS4       
## CS2 "4.71e-11" " NA"      " NA"      " NA"     
## CS3 "5.11e-10" "7.79e-01" " NA"      " NA"     
## CS4 "1.68e-10" "5.93e-01" "7.79e-01" " NA"     
## CS5 "3.18e-01" "1.23e-12" "1.27e-11" "1.27e-11"
```

![unnamed-chunk-29-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-29-186542957.png)



查看具体结果：


```r
head(tmb.brca$TMB.dat)
##             samID variants       TMB   log10TMB Subtype
## 570 BRCA-A1EW-01A        9 0.2368421 0.09231426     CS1
## 558 BRCA-A1IO-01A       11 0.2894737 0.11041248     CS1
## 560 BRCA-A0C3-01A       11 0.2894737 0.11041248     CS1
## 541 BRCA-A1NG-01A       14 0.3684211 0.13621975     CS1
## 531 BRCA-A1Y2-01A       15 0.3947368 0.14449227     CS1
## 525 BRCA-A5RY-01A       16 0.4210526 0.15261016     CS1
```


### 比较FGA

FGA在上面解释过了。

`compFGA()`不仅计算FGA，还会在每个亚型中针对每个样本计算特定的gain（FGG）或loss（FGL）。

需要修改列名为以下名字：


```r
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")

head(segment)
##          sample chrom     start       end   value
## 1 BRCA-A090-01A     1   3301765  54730235 -0.1271
## 2 BRCA-A090-01A     1  54730247  57443819 -0.0899
## 3 BRCA-A090-01A     1  57448465  57448876 -1.1956
## 4 BRCA-A090-01A     1  57448951  64426102 -0.1009
## 5 BRCA-A090-01A     1  64426648 106657734 -0.1252
## 6 BRCA-A090-01A     1 106657854 106835667  0.1371
```

然后运行即可：


```r
# compare FGA, FGG, and FGL
fga.brca <- compFGA(moic.res     = cmoic.brca,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")
## --2 samples mismatched from current subtypes.
## 5% 10% 15% 21% 26% 31% 36% 41% 46% 51% 57% 62% 67% 72% 77% 82% 88% 93% 98%
```

![unnamed-chunk-32-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-32-186542957.png)




```r
head(fga.brca$summary)
##           samID       FGA        FGG        FGL Subtype
## 1 BRCA-A03L-01A 0.6217991 0.30867268 0.31312638     CS1
## 2 BRCA-A04R-01A 0.2531019 0.09132014 0.16178176     CS2
## 3 BRCA-A075-01A 0.7007067 0.41444237 0.28626433     CS1
## 4 BRCA-A08O-01A 0.6501287 0.45648145 0.19364725     CS3
## 5 BRCA-A0A6-01A 0.1468893 0.06356488 0.08332444     CS4
## 6 BRCA-A0AD-01A 0.1722214 0.03864521 0.13357618     CS2
```

### 比较药物敏感性

比较不同亚型的药物敏感性。基于`pRRophetic`包。

药敏分析我们之前也介绍过了：

- []()
- []()


```r
# drug sensitivity comparison
drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr=fpkm[,cmoic.brca$clust.res$samID],#确保样本顺序一致
                         drugs=c("Cisplatin", "Paclitaxel"), #选择药物
                         tissueType  = "breast", # 选择组织类型
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 
## --all samples matched.
## --log2 transformation done for expression data.
## Cisplatin done...
## Cisplatin: Kruskal-Wallis rank sum test p value = 1.38e-58
## post-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:
##     CS1        CS2        CS3        CS4       
## CS2 "5.50e-18" " NA"      " NA"      " NA"     
## CS3 "3.90e-13" "9.78e-02" " NA"      " NA"     
## CS4 "5.21e-01" "4.54e-24" "3.17e-18" " NA"     
## CS5 "2.29e-13" "3.89e-34" "2.57e-30" "3.15e-14"
## Paclitaxel done...
## Paclitaxel: Kruskal-Wallis rank sum test p value = 1.91e-13
## post-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:
##     CS1        CS2        CS3        CS4       
## CS2 "8.11e-10" " NA"      " NA"      " NA"     
## CS3 "3.87e-09" "8.56e-01" " NA"      " NA"     
## CS4 "9.12e-11" "6.59e-01" "8.55e-01" " NA"     
## CS5 "2.82e-04" "5.20e-02" "4.92e-02" "1.95e-02"
```

![unnamed-chunk-34-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-34-186542957.png)

![unnamed-chunk-34-286542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-34-286542957.png)

查看具体结果：


```r
head(drug.brca$Cisplatin)
##               Est.IC50 Subtype
## BRCA-A03L-01A 3.802060     CS1
## BRCA-A075-01A 3.834646     CS1
## BRCA-A0AW-01A 3.126744     CS1
## BRCA-A0AZ-01A 3.303054     CS1
## BRCA-A0BC-01A 3.409402     CS1
## BRCA-A0BF-01A 3.358053     CS1
```


### 比较结果的一致性

>目前，许多癌症都有传统的分类方法，评估新亚型与先前分类的一致性对于反映聚类分析的稳健性和确定潜在的新亚型至关重要。为了衡量当前亚型与其他已存在的分类之间的一致性（相似性），MOVICS提供了`compAgree()`函数来计算四个统计量：Rand指数（RI）、调整的互信息（AMI）、Jaccard指数（JI）和Fowlkes-Mallows指数（FM）；所有这些度量范围从0到1，值越大，两个评估之间的相似度就越高。该函数还可以生成一张液态图（alluvial-diagram），以当前亚型作为参考，可视化两个评估与当前亚型的一致性。


```r
# customize the factor level for pstage
surv.info$pstage <- factor(surv.info$pstage, levels = c("TX","T1","T2","T3","T4"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.brca,
                        subt2comp = surv.info[,c("PAM50","pstage")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
## --all samples matched.
```

![unnamed-chunk-36-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-36-186542957.png)


查看4个统计量：


```r
print(agree.brca)
##   current.subtype other.subtype        RI         AMI        JI        FM
## 1         Subtype         PAM50 0.6929744 0.399563845 0.2910887 0.4694389
## 2         Subtype        pstage 0.5506751 0.004951334 0.1565996 0.2884961
```

以上是第二部分的主要功能。

## RUN Module

接下来是第3部分的功能介绍。

### 差异分析

支持`edgeR`、`DESeq2`、`limma`3种差异分析方法。

对于多组数据，比如这个示例是分成5组，会自动进行每一个组和其他4组的比较。


```r
# run DEA with edgeR
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA") # prefix of figure name
## --all samples matched.
## --you choose edger and please make sure an RNA-Seq count data was provided.
## edger of CS1_vs_Others exists and skipped...
## edger of CS2_vs_Others exists and skipped...
## edger of CS3_vs_Others exists and skipped...
## edger of CS4_vs_Others exists and skipped...
## edger of CS5_vs_Others exists and skipped...

# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count, # deseq2也需要count
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA")
## --all samples matched.
## --you choose deseq2 and please make sure an RNA-Seq count data was provided.
## deseq2 of CS1_vs_Others done...
## deseq2 of CS2_vs_Others done...
## deseq2 of CS3_vs_Others done...
## deseq2 of CS4_vs_Others done...
## deseq2 of CS5_vs_Others done...

# run DEA with limma
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA")
## --all samples matched.
## --you choose limma and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided.
## --log2 transformation done for expression data.
## limma of CS1_vs_Others done...
## limma of CS2_vs_Others done...
## limma of CS3_vs_Others done...
## limma of CS4_vs_Others done...
## limma of CS5_vs_Others done...
```

结果会自动保存在当前工作目录下。

### 识别特定分子

>在这个过程中，按照log2FoldChange排序选择差异表达最大的基因作为每个亚型的生物标志物（默认情况下每个亚型选择200个生物标志物）。这些生物标志物应该通过显著性阈值（例如，名义P值<0.05和调整后的P值<0.05），并且不能与其他亚型识别出的生物标志物重叠。

这一步需要使用上一步得到的结果，比如下面使用`edgeR`识别上调的100个基因：


```r
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-BRCA", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
## --all samples matched.
## --log2 transformation done for expression data.
```

![unnamed-chunk-39-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-39-186542957.png)




```r
# check the upregulated biomarkers
head(marker.up$templates)
##     probe class dirct
## 1    PNMT   CS1    up
## 2 AKR1B15   CS1    up
## 3    DLK1   CS1    up
## 4    ACE2   CS1    up
## 5  CRISP3   CS1    up
## 6   ACSM1   CS1    up
```


使用`limma`识别下调的50个基因：


```r
# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "TCGA-BRCA",
                       dirct         = "down",
                       n.marker      = 50, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = fpkm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")
## --all samples matched.
## --log2 transformation done for expression data.
```

![unnamed-chunk-41-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-41-186542957.png)


### GSEA

富集分析我们写过很多了：

[富集分析常见类型](https://mp.weixin.qq.com/s/RtF7DPXYaObiDauIQTnkFg)
[enrichplot富集分析可视化](https://mp.weixin.qq.com/s/1mpoaZqdgymhSsMGFrCP_A)
[GSEA富集分析可视化](https://mp.weixin.qq.com/s/cusiasAAVPBq3uIHP0EKZw)
[Goplot富集分析可视化](https://mp.weixin.qq.com/s/DckdtQcPv48DDLyA6oZQew)
[GseaVis富集分析可视化](https://mp.weixin.qq.com/s/hdGkcemBdRuayA2ySMH3hw)
[simplifyEnrichment的使用示例](https://mp.weixin.qq.com/s/BmROSJCTEzHRj9yiM8rcmA)
[单基因富集分析](https://mp.weixin.qq.com/s/q6nkujgTYlbOQpkENjyyxA)
[GSVA和ssGSEA](https://mp.weixin.qq.com/s/aUEP6XnejtHohaPeeEOMOQ)

这里用的是`C5`进行GSEA：


```r
# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
```

会自动在每个亚型内进行GSEA分析，我们使用`edger`的结果：


```r
# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "edger", # name of DEA method
                   prefix       = "TCGA-BRCA", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")
## --all samples matched.
## GSEA done...
## --log2 transformation done for expression data.
## Estimating GSVA scores for 50 gene sets.
## Estimating ECDFs with Gaussian kernels
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=                                                                     |   2%
  |                                                                                 
  |======================================================================| 100%
## gsva done...
## heatmap done...
```

![unnamed-chunk-43-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-43-186542957.png)

查看第一个亚型的部分结果：


```r
print(gsea.up$gsea.list$CS1[1:6,3:6])
##                                                     setSize enrichmentScore
## GO_CORNIFICATION                                         51      -0.7747559
## GO_KERATINIZATION                                        58      -0.7392665
## GO_ADULT_LOCOMOTORY_BEHAVIOR                             57      -0.7048011
## GO_REGULATION_OF_ACTIVIN_RECEPTOR_SIGNALING_PATHWAY      16      -0.8535591
## GO_HEMIDESMOSOME_ASSEMBLY                                12      -0.8908059
## GO_WALKING_BEHAVIOR                                      21      -0.8153693
##                                                           NES      pvalue
## GO_CORNIFICATION                                    -2.117694 0.001219512
## GO_KERATINIZATION                                   -2.056651 0.001197605
## GO_ADULT_LOCOMOTORY_BEHAVIOR                        -1.956867 0.001200480
## GO_REGULATION_OF_ACTIVIN_RECEPTOR_SIGNALING_PATHWAY -1.928063 0.001351351
## GO_HEMIDESMOSOME_ASSEMBLY                           -1.924928 0.001428571
## GO_WALKING_BEHAVIOR                                 -1.917848 0.001333333
```

查看每个通路中不同亚型的富集分数：


```r
head(round(gsea.up$grouped.es,3))
##                                                          CS1    CS2    CS3
## GO_REGULATION_OF_MONOCYTE_CHEMOTAXIS                   0.392 -0.474 -0.329
## GO_INDOLE_CONTAINING_COMPOUND_METABOLIC_PROCESS        0.191 -0.457 -0.425
## GO_BENZENE_CONTAINING_COMPOUND_METABOLIC_PROCESS       0.153 -0.092  0.142
## GO_POSITIVE_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION   0.358 -0.430 -0.236
## GO_AMINE_CATABOLIC_PROCESS                             0.174 -0.379 -0.204
## GO_REGULATION_OF_CELLULAR_AMINO_ACID_METABOLIC_PROCESS 0.332  0.053 -0.278
##                                                           CS4    CS5
## GO_REGULATION_OF_MONOCYTE_CHEMOTAXIS                    0.260  0.040
## GO_INDOLE_CONTAINING_COMPOUND_METABOLIC_PROCESS         0.144  0.430
## GO_BENZENE_CONTAINING_COMPOUND_METABOLIC_PROCESS        0.154 -0.503
## GO_POSITIVE_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION    0.223 -0.012
## GO_AMINE_CATABOLIC_PROCESS                              0.192  0.132
## GO_REGULATION_OF_CELLULAR_AMINO_ACID_METABOLIC_PROCESS -0.524  0.440
```

也可以使用`edger`的结果：


```r
# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2",
                   prefix       = "TCGA-BRCA",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 
```


### GSVA

基因集变异分析，和上面的GSEA分析一样的流程：


```r
# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
```




```r
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)
## --all samples matched.
## --log2 transformation done for expression data.
## Estimating GSVA scores for 21 gene sets.
## Estimating ECDFs with Gaussian kernels
## 
  |                                                                            
  |                                                                      |   0%       
  |===                                                                   |   5%
  |======================================================================| 100%
## gsva done...
```

![unnamed-chunk-48-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-48-186542957.png)



```r
# check raw enrichment score
print(gsva.res$raw.es[1:3,1:3])
##                    BRCA-A03L-01A BRCA-A04R-01A BRCA-A075-01A
## Adhesion              0.09351280    -0.2125523    0.10362446
## Antigen_Processing    0.05451682    -0.3909617    0.05508321
## B-Cell_Functions     -0.01931423    -0.6120864    0.18301567
```


### NTP预测外部数据

>哦，等等，我们是否忘记了什么？是的，还有一个尚未使用的数据集，让我们看看我们是否可以使用这些亚型特异性的生物标志物来验证外部Yau队列中的当前乳腺癌亚型。在这部分中，我们的核心目的是预测外部数据集中每个样本的可能亚型。在大多数情况下，这是一个多分类问题，并且在外部队列中，识别出的生物标志物可能很难整体匹配，因此使用基于模型的预测算法可能不可靠。因此，MOVICS为验证队列中的亚型预测提供了两种无模型方法。首先，MOVICS切换到最近模板预测（NTP）方法，该方法可以灵活应用于跨平台、跨物种和多类别的预测，而无需进行任何分析参数的优化。唯一需要做的就是生成一个模板文件，幸运的是，这已经准备好了。


```r
# run NTP in Yau cohort by using up-regulated biomarkers
yau.ntp.pred <- runNTP(expr       = brca.yau$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU") 
## --original template has 500 biomarkers and 272 are matched in external expression profile.
## cosine correlation distance
## 682 samples; 5 classes; 39-66 features/class
## serial processing; 1000 permutation(s)...
## predicted samples/class (FDR<0.05)
```

![unnamed-chunk-50-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-50-186542957.png)

```
## 
##  CS1  CS2  CS3  CS4  CS5 <NA> 
##  104   85   90  120  140  143
```



```r
head(yau.ntp.pred$ntp.res)
##     prediction  d.CS1  d.CS2  d.CS3  d.CS4  d.CS5 p.value    FDR
## 107        CS2 0.7344 0.5165 0.7377 0.7471 0.7512   0.001 0.0020
## 109        CS1 0.5816 0.7581 0.7172 0.7308 0.7295   0.001 0.0020
## 11         CS1 0.6527 0.7307 0.7430 0.7557 0.7731   0.001 0.0020
## 110        CS2 0.7619 0.5418 0.7583 0.7655 0.7567   0.001 0.0020
## 111        CS1 0.6806 0.7157 0.7198 0.7889 0.7911   0.007 0.0106
## 113        CS4 0.6785 0.7505 0.7394 0.6389 0.7308   0.007 0.0106
```


在yau队列中做生存分析：


```r
# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = brca.yau$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 
## --a total of 682 samples are identified.
```

![unnamed-chunk-52-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-52-186542957.png)

```r
print(surv.yau)
## $fitd
## Call:
## survdiff(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
##     na.action = na.exclude)
## 
##               N Observed Expected (O-E)^2/E (O-E)^2/V
## Subtype=CS1 136       51     43.8     1.169     1.452
## Subtype=CS2 117       51     38.5     4.037     4.876
## Subtype=CS3 119       33     42.3     2.045     2.517
## Subtype=CS4 159       43     57.3     3.590     4.820
## Subtype=CS5 151       50     46.0     0.351     0.441
## 
##  Chisq= 11.2  on 4 degrees of freedom, p= 0.02 
## 
## $fit
## Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
##     na.action = na.exclude, error = "greenwood", type = "kaplan-meier", 
##     conf.type = "plain")
## 
##       n events median 0.95LCL 0.95UCL
## CS1 136     51     NA      NA      NA
## CS2 117     51    205     109      NA
## CS3 119     33    236      NA      NA
## CS4 159     43    222     191      NA
## CS5 151     50     NA      NA      NA
## 
## $xyrs.est
## [1] "[Not Available]: argument of xyrs.est was not specified."
## 
## $overall.p
## [1] 0.02410534
## 
## $pairwise.p
## 
## 	Pairwise comparisons using Log-Rank test 
## 
## data:  mosurv.res and Subtype 
## 
##     CS1   CS2   CS3   CS4  
## CS2 0.699 -     -     -    
## CS3 0.158 0.059 -     -    
## CS4 0.100 0.039 0.901 -    
## CS5 0.804 0.504 0.244 0.158
## 
## P value adjustment method: BH
```

比较一致性：


```r
# compare agreement in Yau cohort
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = brca.yau$clin.info[, "PAM50", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH PAM50")
## --all samples matched.
```

![unnamed-chunk-54-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-54-186542957.png)



```r
print(agree.yau)
##   current.subtype other.subtype        RI       AMI        JI        FM
## 1         Subtype         PAM50 0.7830213 0.3819096 0.3318526 0.4994425
```


### PAM预测外部数据

>除了NTP方法，MOVICS还提供了另一种无模型方法来预测亚型。具体来说，首先在发现（训练）队列（即TCGA-BRCA队列）中使用PAM（partition around medoids）分类器来训练模型，以预测验证（测试）队列（即BRCA-Yau队列）中患者的亚型，并将验证队列中的每个样本分配给与其质心具有最高的Pearson相关性的亚型标签17。最后，将使用in-group proportion (IGP) 统计量来评估发现队列和验证队列之间获得的亚型的相似性和可靠性。


```r
yau.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = cmoic.brca,
                       test.expr   = brca.yau$mRNA.expr)
## --all samples matched.
## --a total of 7303 genes shared and used.
## --log2 transformation done for training expression data.
## --testing expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.
```

查看IGP：


```r
print(yau.pam.pred$IGP)
##       CS1       CS2       CS3       CS4       CS5 
## 0.4545455 0.6416667 0.5616438 0.7814570 0.9225806
```


### 一致性评价

想要知道在使用发现队列时，NTP或PAM的准确性如何？想要知道不同预测结果的一致性如何？可以使用runKappa()函数来进行评估。


```r
# predict subtype in discovery cohort using NTP
tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 
## --original template has 500 biomarkers and 500 are matched in external expression profile.
## cosine correlation distance
## 643 samples; 5 classes; 100-100 features/class
## serial processing; 1000 permutation(s)...
## predicted samples/class (FDR<0.05)
## 
##  CS1  CS2  CS3  CS4  CS5 <NA> 
##   99  105  138  155  107   39
#> --original template has 500 biomarkers and 500 are matched in external expression profile.
#> cosine correlation distance
#> 643 samples; 5 classes; 100-100 features/class
#> serial processing; 1000 permutation(s)...
#> predicted samples/class (FDR<0.05)
#> 
#>  CS1  CS2  CS3  CS4  CS5 <NA> 
#>   99  105  138  155  107   39

# predict subtype in discovery cohort using PAM
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = cmoic.brca,
                        test.expr   = fpkm)
## --all samples matched.
## --a total of 13771 genes shared and used.
## --log2 transformation done for training expression data.
## --log2 transformation done for testing expression data.
#> --all samples matched.
#> --a total of 13771 genes shared and used.
#> --log2 transformation done for training expression data.
#> --log2 transformation done for testing expression data.

# check consistency between current and NTP-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
```

![unnamed-chunk-58-186542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-58-186542957.png)

```r

# check consistency between current and PAM-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")
```

![unnamed-chunk-58-286542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-58-286542957.png)

```r

# check consistency between NTP and PAM-predicted subtype in validation Yau-BRCA
runKappa(subt1     = yau.ntp.pred$clust.res$clust,
         subt2     = yau.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR YAU")
```

![unnamed-chunk-58-386542957](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-58-386542957.png)


以上就是`MOVICS`的全部内容了，非常丰富，提供分子分型的绝大多数分析，只需要提供正确的数据即可。

有点类似于之前介绍过的免疫分析一站式R包：`IOBR`，因为集成了10种分子分型的方法，但是同时也提供了完备的下游分析函数，所以又和甲基化分析的一站式R包`CHAMP`类似。

