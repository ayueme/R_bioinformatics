前面的推文给大家介绍了3种识别免疫相关`lncRNA`的方法：[免疫相关lncRNA的识别](https://mp.weixin.qq.com/s/jrgZ6brGyrh1cAnW6Ddp3w)

今天再给大家介绍下如何使用`WGCNA`识别免疫相关`lncRNA`，也算是`WGCNA`的实战教程。

`WGCNA`可以看做是一种筛选候选分子的方法，也可以是识别特定性状相关分子的方法，还可以看做是一种降维方法。

`WGCNA`的原理和入门教程我们就不讲了，大家自己学习一下，直接进入实操。

[toc]

## 准备数据

为了方便，我们还是用之前推文中的数据：TCGA-COAD和TCGA-READ的数据，已经进行过批次矫正了。数据处理详情请见：[]()


```r
rm(list = ls())

load(file = "step3_output.rdata")
load(file = "step1_expr_lnc.rdata")

clin_sub <- clin_sub[match(colnames(expr_lnc),rownames(clin_sub)),]
identical(colnames(expr_lnc),rownames(clin_sub))
```

```
## [1] TRUE
```

```r
expr_lnc[1:4,1:4]
```

```
##        TCGA-D5-6540-01A-11R-1723-07 TCGA-AA-3525-01A-02R-0826-07
## MALAT1                     5.841875                     5.579839
## NORAD                      7.840943                     6.780140
## SNHG6                      7.012464                     6.435600
## SNHG29                     6.309729                     8.097017
##        TCGA-AA-3815-01A-01R-1022-07 TCGA-D5-6923-01A-11R-A32Z-07
## MALAT1                     6.041946                     6.108745
## NORAD                      6.638278                     9.316181
## SNHG6                      6.290823                     7.406282
## SNHG29                     7.249777                     6.167633
```

```r
clin_sub[1:4,1:4]
```

```
##                              status  age gender stage
## TCGA-D5-6540-01A-11R-1723-07  Alive  >65   male     I
## TCGA-AA-3525-01A-02R-0826-07  Alive  >65   male   III
## TCGA-AA-3815-01A-01R-1022-07  Alive  >65 female    II
## TCGA-D5-6923-01A-11R-A32Z-07  Alive <=65   male     I
```

### 过滤lncRNA

首先是过滤低质量的基因，官方建议先自己过滤一下，通过均值、方差、绝对中位差都可以，我们这里采用的是绝对中位差筛选前5000个lncRNA用于后续的分析。不建议通过差异分析筛选基因（lncRNA）。


```r
#生信技能树https://mp.weixin.qq.com/s/DDGlnHr0QtXdU58D4sWhQg
library(WGCNA)
```

```r
datExpr0 = t(expr_lnc[order(apply(expr_lnc,1,mad), decreasing = T)[1:5000],])# top 5000 mad genes

gsg <- goodSamplesGenes(datExpr0)
```

```
##  Flagging genes and samples with too many missing values...
##   ..step 1
```

```r
gsg$allOK
```

```
## [1] TRUE
```

返回`TRUE`说明没问题，返回`FALSE`就用以下代码删除不符合条件的lncRNA。


```r
# 不OK就需要筛选
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

gsg <- goodSamplesGenes(datExpr0)
gsg$allOK
```

### 过滤样本

然后通过聚类树查看是否有离群样本需要剔除


```r
sampleTree <- hclust(dist(datExpr0), method = "average")

sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red") # 我这里选择的在200的高度砍一刀
#dev.off()
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604183933516.png)

砍完一刀后我们选择需要留下的样本，我这里选择了第2个cluster，另外的都被我都不要了。


```r
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
```

```
## clust
##   0   1   2 
##   1 638  11
```

```r
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

nGenes
```

```
## [1] 5000
```

```r
nSamples
```

```
## [1] 638
```

`datExpr`就是准备好的转置后的表达矩阵了，接下来就用它进行后续的分析。

### 准备性状数据（临床信息）

然后是准备样本数据（临床信息，性状数据），用于后续计算性状和模块/基因之间的相关性。

如果是分类数据需要都变成数值型，其实还有`WGCNA`包中还介绍了其他方法，以后我们专门再写一篇推文比较下。


```r
datTraits <- clin_sub[keepSamples,] # 筛选样本

#分类数据变为0,1...还是1,2...并没有影响，不信你可以试一下
datTraits$cluster <- as.numeric(factor(datTraits$cluster))
datTraits$status <- as.numeric(datTraits$status)
datTraits$age <- as.numeric(datTraits$age)
datTraits$gender <- as.numeric(datTraits$gender)
datTraits$stage <- as.numeric(datTraits$stage)
datTraits$msi <- as.numeric(datTraits$msi)

str(datTraits)# 全都变成数值型
```

```
## 'data.frame':	638 obs. of  6 variables:
##  $ status : num  1 1 1 1 1 1 2 1 1 1 ...
##  $ age    : num  2 2 2 1 2 2 1 2 2 2 ...
##  $ gender : num  2 2 1 2 2 2 1 1 1 1 ...
##  $ stage  : num  1 3 2 1 3 1 3 3 3 3 ...
##  $ msi    : num  4 1 1 4 4 2 4 4 1 4 ...
##  $ cluster: num  1 1 1 1 1 2 1 2 2 2 ...
```

```r
identical(rownames(datTraits),rownames(datExpr))
```

```
## [1] TRUE
```

接下来是看看样本聚类和临床特征之间的关系：


```r
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE) #作者推荐选signed

# Plot the sample dendrogram and the colors underneath.
#pdf(file = "Plots/sampleClustering_addTraits.pdf", width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184025486.png)

```r
#dev.off()
```

到这里是数据准备的过程，这样我们的表达矩阵数据和临床信息就整理好了。


```r
#保存下数据
save(datExpr,datTraits,file = "wgcna-01-dataInput.rdata")
```

## 网络构建和模块识别

接下来是网络构建和模块识别，作者提供了3种方式：

- 一步法
- 分步法
- block-wise法，适合超大数据

但是不管哪种方法，都是从挑选合适的软阈值开始的，这一步是一样的。

### 挑选软阈值


```r
# 加载数据
rm(list = ls())
load(file = "wgcna-01-dataInput.rdata")
```

软阈值一般要求R^2在0.9以上，最小也要在0.8以上。下面是挑选的代码，基本不用改。


```r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
```

```
## pickSoftThreshold: will use block size 5000.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 5000 of 5000
```

```
## Warning: executing %dopar% sequentially: no parallel backend registered
```

```
##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k.   max.k.
## 1      1    0.289 -1.32          0.965 6.01e+02  5.71e+02 1290.000
## 2      2    0.721 -2.01          0.976 1.19e+02  9.84e+01  481.000
## 3      3    0.792 -2.30          0.971 3.12e+01  2.10e+01  215.000
## 4      4    0.854 -2.27          0.988 9.96e+00  5.25e+00  107.000
## 5      5    0.874 -2.17          0.991 3.68e+00  1.42e+00   57.800
## 6      6    0.878 -2.09          0.972 1.52e+00  4.35e-01   33.000
## 7      7    0.875 -2.05          0.971 6.89e-01  1.41e-01   19.600
## 8      8    0.412 -2.67          0.410 3.36e-01  4.93e-02   12.000
## 9      9    0.416 -2.57          0.415 1.75e-01  1.77e-02    7.620
## 10    10    0.421 -2.49          0.419 9.61e-02  6.59e-03    4.930
## 11    12    0.932 -1.71          0.979 3.33e-02  9.95e-04    2.180
## 12    14    0.407 -2.17          0.385 1.36e-02  1.67e-04    1.090
## 13    16    0.931 -1.55          0.932 6.39e-03  2.99e-05    0.721
## 14    18    0.954 -1.46          0.946 3.39e-03  5.46e-06    0.530
## 15    20    0.390 -1.74          0.221 2.00e-03  1.04e-06    0.449
```

```r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1]+2, -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184139116.png)

结果就是看这两个图，选择达到0.9以上的软阈值，右边的图是连接度，选择逐渐平稳的阈值。

是在搞不清楚就按照下面的方法选：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604160951297.png)

代码也会自动给出一个合适的软阈值：


```r
sft$powerEstimate
```

```
## [1] 4
```

推荐选4，但是看图我们还是选6。

### 一步法构建网络和识别模块

一步法构建网络和识别模块。计算时用到了`BLAS`，所以这一步是可以提速的，R自带的`blas`非常慢，提速可以看这个：[让你的R语言提速100倍](https://mp.weixin.qq.com/s/qt3QJIfP3yzQdBnHHxQtMg)


```r
allowWGCNAThreads()#RStudio不支持enableWGCNAThreads()
```

```
## Allowing multi-threading with up to 24 threads.
```

```r
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,#每个模块最少有多少个基因
                      
                       #根据官方的说明，现在电脑内存普遍都在8G以上，20000也完全没问题！
                       maxBlockSize = 20000, #大内存可以设大一点，这个意思是基因数量超过20000就分批处理
                       reassignThreshold = 0, mergeCutHeight = 0.25,#合并模块的阈值
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM", #这个名字就不改了
                       verbose = 3)
```

```
##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will not use multithreading.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file femaleMouseTOM-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 1 genes from module 4 because their KME is too low.
##      ..removing 1 genes from module 7 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.25
##        Calculating new MEs...
```

```r
class(net)
```

```
## [1] "list"
```

```r
names(net)
```

```
##  [1] "colors"         "unmergedColors" "MEs"            "goodSamples"   
##  [5] "goodGenes"      "dendrograms"    "TOMFiles"       "blockGenes"    
##  [9] "blocks"         "MEsOK"
```

```r
table(net$colors)
```

```
## 
##    0    1    2    3    4 
## 2995 1665  226   66   48
```

一共5个模块，其中0是不属于任何模块的基因（lncRNA），一共5000个lncRNA，2995个不属于任何模块......

下面是的代码是可视化基因（lncRNA）聚类树和模块：


```r
# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184400574.png)

这一步做好了也保存下数据。


```r
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs,moduleLabels,moduleColors,geneTree,
     file = "wgcna-02-networkConstruction-auto.rdata")
```

官方教程写“recutBlockwiseTrees可用于重新计算，比从头计算节省时间”，但没有给示例，而且好像没有人这么做，难道大家用的都不是一步法吗？还是说对一步法的结果都很满意？

下面的代码供参考：


```r
tomfiles <- list.files(path = ".",pattern = "^femaleMouseTOM")
tomfiles

# 还有很多参数可以控制模块数量，大家可以看帮助文档
aa <- recutBlockwiseTrees(datExpr = datExpr, goodSamples = gsg$goodSamples, 
                          goodGenes = gsg$goodGenes,
                          TOMFiles = tomfiles, dendrograms = net$dendrograms,
                          minModuleSize = 100,
                          blocks = net$blocks)
```

### 分步法构建网络

挑选软阈值的部分完全一样，这里就不写了，直接从软阈值选择6开始。


```r
rm(list = ls())
load(file = "wgcna-01-dataInput.rdata")
```

分步法构建网络的第一步，计算连接矩阵：


```r
softPower = 6
adjacency = adjacency(datExpr, power = softPower)
```

第二步，将连接矩阵变为拓扑重叠矩阵（Topological Overlap Matrix (TOM)）：


```r
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
```

```
## ..connectivity..
## ..matrix multiplication (system BLAS)..
## ..normalization..
## ..done.
```

```r
dissTOM = 1-TOM
```

可视化基因TOM矩阵的聚类树：


```r
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.1);
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184451078.png)

聚在一起的说明分子间高度相关。

接下来是切割树，也就是识别模块的过程，通过`cutreeDynamic`函数实现，可以控制的参数非常多，这里就选了`deepSplit`和`minClusterSize`，其他的大家可以自己探索下，但其实没必要都用上。


```r
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, #范围0-4，越大模块越多
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize #每个模块最少基因数量
                            )
```

```
##  ..cutHeight not given, setting it to 1  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
table(dynamicMods)
```

```
## dynamicMods
##    0    1    2    3    4    5    6    7    8 
## 1180 1035  982  608  401  295  192  191  116
```

这个结果明显比一步法识别的模块多，用上的lncRNA也多了。

然后是可视化切割后的聚类树和模块：


```r
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
```

```
## dynamicColors
##     black      blue     brown     green      grey      pink       red turquoise 
##       191       982       608       295      1180       116       192      1035 
##    yellow 
##       401
```

```r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184528561.png)

第三步，合并高度相似的模块。

计算模块的`eigengenes`（其实就是第一主成分），再通过相关性进行聚类，达到量化模块间相似性的目的，方便进行模块合并：


```r
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 # 要根据你自己画出来的图选择
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184552495.png)

上面这个图我选择了在高度为0.25的地方进行切割，根据上面这个图，切割后`blue`、`red`、`brown`就会被合并为一个模块了！大家要根据自己画出来的这个聚类树的图选择合适的标准。

选择切割高度进行合并，注意，这里的 相关性 = 1-切割高度，比如切割高度是0.25，那就是选择相关性大于0.75的模块进行合并。


```r
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.25
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 9 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 8 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 7 module eigengenes in given set.
##    Calculating new MEs...
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 7 module eigengenes in given set.
```

```r
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
```

然后又是一个可视化，把合并和未合并的模块画在一起：


```r
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184634685.png)

可以看到有些模块被合并到一起了，根据上面的那个图也能看出来~


```r
#合并完还剩7个模块，意料之中
table(mergedColors)
```

```
## mergedColors
##     black      blue     green      grey      pink turquoise    yellow 
##       191      1782       295      1180       116      1035       401
```

最后是保存数据：


```r
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "wgcna-02-networkConstruction-stepByStep.rdata")
```

### block-wise 略

## 模块和性状（临床信息）关联

加载之前的数据：


```r
rm(list = ls())
load(file = "wgcna-02-networkConstruction-stepByStep.rdata")
load(file = "wgcna-01-dataInput.rdata")
```

计算模块（使用eigengenes代表）和性状（临床信息）的相关性和P值：


```r
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#datTraits$cluster <- ifelse(datTraits$cluster==2,1,0)#变成0,1和变成1,2没有影响

moduleTraitCor = cor(MEs, datTraits, use = "p")#这个cor是WGCNA::cor，可以用于计算任意两个矩阵的每一列之间的相关性（比如500个lncRNA和1000个mRNA），很实用哦！
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
```

画图表示，这个图应该是出现频率最高的了，做`WGCNA`一般都会出现这个图：


```r
sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
textMatrix[1:6,1:6]
```

```
##      [,1]            [,2]             [,3]            [,4]            
## [1,] "0.0023\n(1)"   "0.018\n(0.6)"   "0.034\n(0.4)"  "0.012\n(0.8)"  
## [2,] "0.073\n(0.07)" "0.011\n(0.8)"   "0.0065\n(0.9)" "0.063\n(0.1)"  
## [3,] "0.079\n(0.05)" "-0.1\n(0.01)"   "0.081\n(0.04)" "0.15\n(2e-04)" 
## [4,] "0.036\n(0.4)"  "-0.054\n(0.2)"  "0.026\n(0.5)"  "0.046\n(0.2)"  
## [5,] "-0.055\n(0.2)" "-0.0021\n(1)"   "0.0099\n(0.8)" "-0.076\n(0.06)"
## [6,] "0.032\n(0.4)"  "-0.098\n(0.01)" "0.014\n(0.7)"  "0.12\n(0.002)" 
##      [,5]            [,6]           
## [1,] "-0.3\n(3e-15)" "0.074\n(0.06)"
## [2,] "0.055\n(0.2)"  "0.065\n(0.1)" 
## [3,] "0.48\n(3e-38)" "0.052\n(0.2)" 
## [4,] "0.31\n(9e-16)" "0.014\n(0.7)" 
## [5,] "-0.00057\n(1)" "-0.019\n(0.6)"
## [6,] "0.22\n(1e-08)" "-0.022\n(0.6)"
```

```r
par(mar = c(9, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184748762.png)

同时展示性状（临床信息）和模块的相关系数和P值。

## 提取感兴趣的基因（lncRNA）

上面的几个步骤其实是`WGCNA`的核心步骤，剩下的其实都是基于上面这个结果的个性化分析。

因为我们最终还是要找到几个重要的分子，所以做了这么大一通分析，说到底就是在不断缩小候选分子范围。

只要你提供的性状信息不同，理论上你可以通过这个方法寻找和任何特性相关的分子。

### 提取某个模块的基因

如果你要提取某个模块的基因，可以直接通过以下代码实现：


```r
module = "blue"
module_genes <- colnames(datExpr)[moduleColors==module]
head(module_genes)
```

```
## [1] "AC135388.1" "AP003170.3" "AL121917.2" "AC114291.1" "AL138831.1"
## [6] "AL138960.1"
```

如果要提取所有模块的基因，那就自己写个循环即可。

但是还可以更进一步，通过计算GS和MM，进一步缩小候选分子范围。

### 计算GS和MM

基因和性状的相关性：gene significance GS
基因和模块的相关性：module membership MM

选择我们感兴趣的性状（临床信息），这里是`cluster`，也就是根据`ssGSEA`免疫浸润的结果而得到的分子亚型，所以只要提取和这个性状信息相关的基因，就得到了免疫相关`lncRNA`！

但是很遗憾我们的结果不太好哈，看上面的热图，并没有和`cluster`显著相关的模块，但是对我们的演示影响不大。


```r
# Define variable weight containing the weight column of datTrait
cluster = as.data.frame(datTraits$cluster) # 我们选cluster
names(cluster) = "cluster"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

# 基因和模块的相关性及P值
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# 基因和性状的相关性，这里是和样本亚型的相关性
geneTraitSignificance = as.data.frame(cor(datExpr, cluster, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cluster), sep="");
names(GSPvalue) = paste("p.GS.", names(cluster), sep="");
```

### 选择GS和MM都很高的分子

我们选择黑色模块，矮子里面拔将军吧......

画个图展示下每个基因的GS和MM的分布情况：


```r
module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184834549.png)

如果这个结果非常好，那这些点应该是大概呈一条直线状的，比如这样的：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604164813996.png)

但是很遗憾我们的数据很不好。

接下来就是提取出GS和MM都很高的基因：


```r
tmp1 <- rownames(geneModuleMembership)[abs(geneModuleMembership[moduleGenes, column])>0.7]
tmp2 <- rownames(geneTraitSignificance)[abs(geneTraitSignificance[moduleGenes, 1])>0.04] #这里太小了，但是我们数据差，没办法！

genes <- unique(intersect(tmp1,tmp2))
length(genes)
```

```
## [1] 104
```

```r
head(genes)
```

```
## [1] "AF124730.1" "AC239584.1" "AC090578.2" "AL121890.5" "AL132655.2"
## [6] "AC068831.5"
```

有了这些基因，后面的各种分析就是大家常见的了，比如各种富集分析、生存分析、构建模型等等。

那加上今天介绍的识别免疫相关lncRNA的方法，就一共有4种方法了！

>这么多方法又有选择困难症了？天真，小孩子才做选择，大人都是全都要！

