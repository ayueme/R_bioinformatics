这篇推文主要探讨下`WGCNA`如何处理分类性状。

`eigengenes`可以代表某个模块，在计算出模块的`eigengenes`后，下一步就是探索`eigengenes`和性状之间的关系，也就是模块和性状之间的关系。

大家见到的比较多的是计算相关性，此时需要性状是数字才行。但是大家的性状有很多分类变量，此时应该如何处理呢？

以下是常规的分类变量处理原则：

- 如果是二分类，只要变为0/1即可（也可以变成1/2，没有影响），或者变成因子型；
- 如果是有序多分类，比如治愈、好转、未愈，这种，可以变成数字1,2,3，或者变成因子型；
- 如果是无序多分类，那么此时需要使用`WGCNA`提供的函数进行处理。

假如我们有一个无序分类变量`x`，它有3组：


```r
library(WGCNA)

x <- rep(c("A","B","C"), each = 3)
x
## [1] "A" "A" "A" "B" "B" "B" "C" "C" "C"
```

我们可以把它变成3组之间**两两比较**的形式，使用的是`binarizeCategoricalVariable()`函数：


```r
out <- binarizeCategoricalVariable(x,includePairwise = T,includeLevelVsAll = F)
data.frame(x, out)
##   x B.vs.A C.vs.A C.vs.B
## 1 A      0      0     NA
## 2 A      0      0     NA
## 3 A      0      0     NA
## 4 B      1     NA      0
## 5 B      1     NA      0
## 6 B      1     NA      0
## 7 C     NA      1      1
## 8 C     NA      1      1
## 9 C     NA      1      1
```

或者变成`1-vs-all`的形式：


```r
out <- binarizeCategoricalVariable(x,includePairwise = F,includeLevelVsAll = T)
data.frame(x, out)
##   x A.vs.all B.vs.all C.vs.all
## 1 A        1        0        0
## 2 A        1        0        0
## 3 A        1        0        0
## 4 B        0        1        0
## 5 B        0        1        0
## 6 B        0        1        0
## 7 C        0        0        1
## 8 C        0        0        1
## 9 C        0        0        1
```

`binarizeCategoricalVariable()`是针对1个变量的，通常我们的性状数据都是包含在1个数据框中的，并且可能同时有多个分类变量，此时可以使用`binarizeCategoricalColumns()`。

比如，对于我们之前用过的`datTraits`这个性状数据，我们假设其中的`stage`和`msi`是无序多分类变量，然后对这两个变量进行转换：


```r
load(file = "../000files/wgcna-02-networkConstruction-stepByStep.rdata")
load(file = "../000files/wgcna-01-dataInput.rdata")

out <- binarizeCategoricalColumns(datTraits,
                                  convertColumns = c("stage","msi"),
                                  includePairwise = T,
                                  includeLevelVsAll = F
                                  )
out[1:4,4:9]
##                              stage.2.vs.1 stage.3.vs.1 stage.4.vs.1
## TCGA-D5-6540-01A-11R-1723-07            0            0            0
## TCGA-AA-3525-01A-02R-0826-07           NA            1           NA
## TCGA-AA-3815-01A-01R-1022-07            1           NA           NA
## TCGA-D5-6923-01A-11R-A32Z-07            0            0            0
##                              stage.5.vs.1 stage.3.vs.2 stage.4.vs.2
## TCGA-D5-6540-01A-11R-1723-07            0           NA           NA
## TCGA-AA-3525-01A-02R-0826-07           NA            1           NA
## TCGA-AA-3815-01A-01R-1022-07           NA            0            0
## TCGA-D5-6923-01A-11R-A32Z-07            0           NA           NA

colnames(out)
##  [1] "status"       "age"          "gender"       "stage.2.vs.1" "stage.3.vs.1"
##  [6] "stage.4.vs.1" "stage.5.vs.1" "stage.3.vs.2" "stage.4.vs.2" "stage.5.vs.2"
## [11] "stage.4.vs.3" "stage.5.vs.3" "stage.5.vs.4" "msi.2.vs.1"   "msi.3.vs.1"  
## [16] "msi.4.vs.1"   "msi.3.vs.2"   "msi.4.vs.2"   "msi.4.vs.3"   "cluster"
```

`datTraits`这个数据在之前的推文里，因为4篇推文都是前后有联系的，所以我都放在这里：

- [批次效应去除之combat和removebatcheffect](https://mp.weixin.qq.com/s/yRUmVTimI9f9itoHWxyYrA)
- [免疫浸润结果分子分型（一致性聚类ConsensusClusterPlus）](https://mp.weixin.qq.com/s/96s_hfBH0HjLvvTfNgTIlQ)
- [免疫相关lncRNA的识别](https://mp.weixin.qq.com/s/jrgZ6brGyrh1cAnW6Ddp3w)
- [WGCNA实战：识别免疫相关lncRNA](https://mp.weixin.qq.com/s/Pr33WscVtNQQaoryxTiJ-Q)

接下来就是计算模块（使用eigengenes代表）和性状（临床信息）之间的相关性和P值：


```r
# 计算模块的eigengenes，也就是第一主成分
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0) # 对列（也就是模块）的顺序重新排序，让相似性大的在一起

# 计算模块和性状的相关系数
# 这个cor是WGCNA::cor，可以计算任意两个矩阵的每列之间的相关性
#（比如500个lncRNA和1000个mRNA），很实用！
moduleTraitCor <- cor(MEs, out, use = "p")

# 计算相关系数的P值
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
```

然后画图就可以了：


```r
sizeGrWindow(10,6)

# 把相关系数和P值放在一起
textMatrix <- paste(signif(moduleTraitCor, 2),
                    "\n(",
                    signif(moduleTraitPvalue, 1), ")", 
                    ep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#textMatrix[1:6,1:6]

par(mar = c(9, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(out),
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

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230917104259080.png)

和没进行转换之前的图形比较一下：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230604184748762.png)

## 参考资料

1. https://peterlangfelder.com/2018/11/25/working-with-categorical-variables/
2. https://www.biostars.org/p/293281/
3. https://support.bioconductor.org/p/111449/#111450

