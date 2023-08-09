关于TCGA的差异分析之前介绍过，不过略微有些不够完整，而且主要是演示的`TCGAbiolinks`这个包，对于`DEseq2`介绍的不够，所以今天专门说一下使用`DEseq2`进行差异分析。

对于TCGA的差异分析，很多初学者很纠结，不知道到底是用counts/tpm/fpkm，到底是用哪个包，我这里给出我的建议，**对于TCGA的差异分析，就用counts，`DEseq2`进行差异分析！**

本文目录：
[toc]

## DESeq2分析整理好的counts数据

表达矩阵的下载和整理这里就不演示了，我们直接使用**1行代码整理好的TCGA-COAD的counts数据**，无缝衔接！

- [1行代码提取6种TCGA表达矩阵和临床信息](https://mp.weixin.qq.com/s/1OBGjUKnGyiALmLafYNPUQ)
- [1行代码提取6种TCGA表达矩阵2.0版](https://mp.weixin.qq.com/s/QFGCtrIeaAIichovw6OBVw)
- [1行代码提取TCGA的6种表达矩阵是有视频教程的](https://mp.weixin.qq.com/s/u6VkBcYqakZkaNXjzNTZcw)


```r
# 加载数据和R包
rm(list = ls())
library(DESeq2)

load(file = "G:/tcga/output_expr/TCGA-COAD_mrna_expr_counts.rdata")
```

我们的1行代码获取的表达矩阵是**提取好的标准表达矩阵的形式**，行是基因，列是样本，可以直接使用，不需要任何修改。

样本分组需要稍作整理。


```r
# 根据第14、15个字符进行分组，01-09是tumor，10-29是normal
group <- ifelse(substr(colnames(mrna_expr_counts),14,15)<10,"tumor","normal")
metadata <- data.frame(sample_id = colnames(mrna_expr_counts),
                       group = group
                       )
table(metadata$group)
## 
## normal  tumor 
##     41    480

head(metadata)
##                      sample_id  group
## 1 TCGA-AA-A03F-01A-11R-A16W-07  tumor
## 2 TCGA-G4-6314-01A-11R-1723-07  tumor
## 3 TCGA-A6-3809-01A-01R-A278-07  tumor
## 4 TCGA-AZ-6605-01A-11R-1839-07  tumor
## 5 TCGA-AZ-6605-11A-01R-1839-07 normal
## 6 TCGA-F4-6569-01A-11R-1774-07  tumor
```

接下来就是`DEseq2`进行差异分析的流程了。

首先构建DDS，需要提供3个参数，表达矩阵-直接使用我们的1行代码得到的表达矩阵即可，无需任何修改。
`colData`是样本名和样本分组组成的数据框，`design`是包含分组信息的列。


```r
## 首先构建DDS

dds1 <- DESeqDataSetFromMatrix(countData = mrna_expr_counts, 
                               colData = metadata, 
                               design = ~ group) 
```

接下来时过滤掉表达量低的基因，这一步，可做可不做，因为`DEseq2`在计算结果会自动进行过滤！做了有以下好处：减少内存占用，加快运行速度，画图时减少意外（表达量很低在组间没有差异可能画不出来）。

比较流行的方法是：如果某个基因在一半以上（或者75%）的样本中表达量都是低于10(这个数字也没有标准答案)，那就过滤掉。


```r
# 我这里没做这一步，写出来给大家做个参考
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds1 <- dds1[keep,]
```

真正的差异分析就1行代码而已：


```r
## 差异分析
dds <- DESeq(dds1)

## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 1548 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing

#save(dds,file="output_expr/mRNA_DE_counts.rdata")
```

下面提取结果就好了，我在之前介绍过`limma`进行差异分析的一些注意事项，主要是谁和谁比的问题，搞不好很容易闹乌龙：[limma差异分析，谁和谁比很重要吗？](https://mp.weixin.qq.com/s/vdkDcBzuoqCASts61efjBw)

`DEseq2`不用在一开始指定，在提取结果时指定也可以，使用起来很方便：


```r
# 提取结果，如果你一开始没有用因子level限定组别顺序，这里可以限定
# 添加tidy=T,返回数据框
res <- results(dds, contrast = c("group","tumor","normal")) # 指定tumor比normal
res
## log2 fold change (MLE): group tumor vs normal 
## Wald test p-value: group tumor vs normal 
## DataFrame with 19938 rows and 6 columns
##             baseMean log2FoldChange     lfcSE      stat      pvalue        padj
##            <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
## MT-CO1        356764     -0.8513972  0.146049 -5.829513 5.55894e-09 1.64448e-08
## MT-ND4        347353     -0.1095464  0.176165 -0.621841 5.34046e-01 5.81407e-01
## MT-CO2        238026     -0.0412895  0.153885 -0.268315 7.88457e-01 8.16530e-01
## MT-CO3        230096     -0.4469429  0.144782 -3.087016 2.02176e-03 3.35341e-03
## ACTB          211370     -0.1826209  0.103648 -1.761938 7.80798e-02 1.02550e-01
## ...              ...            ...       ...       ...         ...         ...
## AC084756.2         0             NA        NA        NA          NA          NA
## AL031178.2         0             NA        NA        NA          NA          NA
## AC119733.1         0             NA        NA        NA          NA          NA
## AL451106.1         0             NA        NA        NA          NA          NA
## AC006486.3         0             NA        NA        NA          NA          NA

# 按照p值排序
resOrdered <- res[order(res$pvalue),]
DEG <- na.omit(resOrdered) # 去掉缺失值
DEG
## log2 fold change (MLE): group tumor vs normal 
## Wald test p-value: group tumor vs normal 
## DataFrame with 18353 rows and 6 columns
##            baseMean log2FoldChange     lfcSE         stat       pvalue
##           <numeric>      <numeric> <numeric>    <numeric>    <numeric>
## CDH3       3130.674        5.97400  0.150985      39.5669  0.00000e+00
## KRT80      1343.214        6.77220  0.184900      36.6263 1.09154e-293
## ETV4       3531.664        5.31186  0.152688      34.7889 3.57919e-265
## ESM1        233.411        5.76247  0.180466      31.9311 9.89966e-224
## FOXQ1      2204.014        6.23267  0.199737      31.2044 9.29639e-214
## ...             ...            ...       ...          ...          ...
## MARCHF6 4176.387150    2.06911e-04 0.0796146  0.002598912     0.997926
## GRXCR2     0.830639    8.49323e-04 0.3678967  0.002308591     0.998158
## PDE4B    534.259863   -4.51214e-04 0.2007717 -0.002247397     0.998207
## ATL2    1915.827954   -7.23339e-05 0.0753469 -0.000960012     0.999234
## FN3KRP  1471.955462   -6.17016e-05 0.0761668 -0.000810084     0.999354
##                 padj
##            <numeric>
## CDH3     0.00000e+00
## KRT80   1.00165e-289
## ETV4    2.18963e-261
## ESM1    4.54221e-220
## FOXQ1   3.41233e-210
## ...              ...
## MARCHF6     0.998144
## GRXCR2      0.998316
## PDE4B       0.998316
## ATL2        0.999288
## FN3KRP      0.999354
```

而且结果会提醒你到底是谁比谁：`log2 fold change (MLE): group tumor vs normal`

到这里差异分析其实就做好了，接下来你可以根据`padj`和`log2FoldChange`选择合适的基因，我们就不再演示了。

下面探索下vst标准化后的数据。

## VST探索

如果是差异分析，别纠结，就用counts，使用`DESeq`进行差异分析，后续的生存分析、相关性分析、火山图、热图、PCA、聚类等分析，可以统统使用vst标准化后的数据，当然你也可以选择log2(tpm+1)。

你看这篇cell的文章用的就是vst后的数据：
![doi:10.1016/j.cell.2018.03.052.](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221229100729382.png)

关于vst这种方法的参考文献也放在下面：

- variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010)

vst标准化后的数据有多种获取方式，可以从最开始的`dds1`提取，可以从运行`DESeq`后的`dds`提取，也可以从表达矩阵直接开始，3种结果完全一样！


```r
# 3种方法完全一样
vsd1 <- assay(vst(dds1))
vsd2 <- assay(vst(dds))
vsd <- vst(as.matrix(mrna_expr_counts))
vsd[1:6,1:3]

identical(vsd,vsd1)
identical(vsd1,vsd2)
```

可以画个箱线图看看vst标准化后的数据表达情况，经过vst转换后的表达矩阵表达量还是很好的：


```r
boxplot(vsd[,1:30])
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230109191108422.png)

看看这个表达量，很棒，这也是为什么vst标准化后的数据可以直接进行后续分析(不需要log)的原因。

### 热图

下面画个热图看看。


```r
library(pheatmap)

# 这里的`dds`需要是运行`DESeq`之后的
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

# 这里dds/dds1都行
df <- as.data.frame(colData(dds)[,"group"])
rownames(df) <- colnames(mrna_expr_counts)

pheatmap(vsd[select,],
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_col = df
         )
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230109191137312.png)

TCGA数据由于正常样本比较少，所以直接画出的热图大都不太好看，看不出来明显的四象限，大家可以自己去GEO下载其他测序数据试试看。

### 火山图

再画个火山图看看。


```r
library(ggplot2)
library(ggrepel)

tmp.volcano <- as.data.frame(DEG)
tmp.volcano$type <- ifelse(tmp.volcano$log2FoldChange > 2 & tmp.volcano$padj < 0.01, "up",
                        ifelse(tmp.volcano$log2FoldChange < -2 & tmp.volcano$padj < 0.01, "down", "not-sig"))
tmp.volcano$gene <- rownames(tmp.volcano)

ggplot(tmp.volcano, aes(log2FoldChange, -log10(padj)))+
  geom_point(aes(color=type))+
  scale_color_manual(values = c("blue","black","red"))+
  geom_hline(yintercept = -log10(0.01),linetype=2)+
  geom_vline(xintercept = c(-2,2), linetype=2)+
  geom_text_repel(data = subset(tmp.volcano, abs(log2FoldChange) > 4), 
                  aes(label=gene),col="black",alpha = 0.8)+
  theme_bw()
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-10-163832021.png)


### PCA

然后是主成分分析PCA的可视化，`DESeq`自带这个功能：


```r
plotPCA(vst(dds), intgroup="group")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-163832021.png)

如果你嫌丑，可以返回数据，自己画：


```r
library(ggplot2)

pcaData <- plotPCA(vst(dds), intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-12-163832021.png)

如果你想添加各种元素，比如置信椭圆、箭头等，可以参考之前的推文自己做PCA然后画图：

- [R语言主成分分析](http://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247495269&idx=1&sn=830f994a751c210e5bf9a5a8758b7d4b&chksm=facadee2cdbd57f41c1e96f2bce962394a3f79a63a67314e6e57b160f34053e79692b0a42745#rd)
- [R语言主成分分析可视化(颜值高，很详细)](http://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247495271&idx=1&sn=23a1967958e17c4d212f69b6056a27dc&chksm=facadee0cdbd57f6b3f567f6ef7192cfde52d648643a647f96093de17c4cc5e3eb27cac78b4b#rd)
- [R语言PCA可视化3D版](http://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247495599&idx=1&sn=a6b11c466e6aeb99fc4cce113e01dc6a&chksm=facadf28cdbd563e786399b0cd5f79bdfdde52a4b7719fd743eea25faf9fc1635d470325268d#rd)
- [使用R语言美化PCA图](http://mp.weixin.qq.com/s?__biz=MzUzOTQzNzU0NA==&mid=2247483965&idx=1&sn=f7217eba8b7aac7402a53fb2b1a36cc6&chksm=fac932bacdbebbac7ae1a875c076bf29ada73b86d75e42c6eb200ffc1a27fed436fbfd248737#rd)

## 参考资料

1. DEseq2官方文档





