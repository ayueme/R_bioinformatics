今天给大家介绍GSEA，基因集富集分析，的可视化！主要是带你**详细了解如何自定义GSEA的经典图形**。

之前的推文已经带大家了解了富集分析的常见类型以及如何使用`clusterprofiler`进行富集分析，如何使用`enrichplot`进行可视化：

- [富集分析常见类型](https://mp.weixin.qq.com/s/RtF7DPXYaObiDauIQTnkFg)
- [enrichplot可视化富集分析结果](https://mp.weixin.qq.com/s/1mpoaZqdgymhSsMGFrCP_A)

不过在上期内容中，我们主要是以`ORA`为例，演示了如何进行可视化，虽然`enrichplot`的画图函数对于`ORA`和`GSEA`都是支持的（参考上期内容），但是我们并没有过多介绍`GSEA`富集分析的可视化内容，主要是因为部分函数是专门针对`GSEA`的。

所以本期内容会详细介绍`GSEA`可视化及如何自定义。

本期目录：

[toc]

## 准备数据

用`gse87466`这个GEO的数据做演示，下载整理的过程这次就不演示了。数据可以免费在粉丝QQ群获取。


```r
library(easyTCGA)
load(file = "G:/easyTCGA_test/gse87466.Rdata")
```

这是一个炎症性肠病的数据集，一共108个样本，21个normal，87个uc（ulcerative colitis）。

探针注释我已经提前做好了，但是有一些探针对应多个symbol，为了方便我这里直接删掉了：


```r
exprSet[1:4,1:4]
##                                           GSM2332098 GSM2332099 GSM2332100
## IGK@ /// IGKC                               13.86197   13.76880   13.95740
##                                             13.95740   13.92619   13.79664
## IGL@                                        13.73797   13.61266   13.86197
## IGH@ /// IGHA1 /// IGHA2 /// LOC100126583   13.79664   13.16844   13.76880
##                                           GSM2332101
## IGK@ /// IGKC                               13.95740
##                                             13.86197
## IGL@                                        13.76880
## IGH@ /// IGHA1 /// IGHA2 /// LOC100126583   13.73797
group <- factor(group_list,levels = c("normal","UC"))
table(group)
## group
## normal     UC 
##     21     87
```

首先对这个数据做下差异分析，也是用`easyTCGA`包，1行代码即可，基因芯片数据也是支持的，并且它会自动检测需不需要进行log2转换，如果是`count`矩阵，会自动使用`DESeq2`、`limma`、`edgeR`进行差异分析，如果不是，会自动进行`wilcoxon`和`limma`的差异分析：


```r
library(easyTCGA)

diff_res <- diff_analysis(exprset = exprSet
                          , group = group
                          , is_count = F
                          )
## log2 transform not needed
## => Running limma
## => Running wilcoxon test
## => Analysis done.

# limma的结果
diff_limma <- diff_res$deg_limma

# 多个gene symbol的直接删除，方便演示
diff_limma <- diff_limma[!grepl("/",diff_limma$genesymbol),]
head(diff_limma)
##               logFC   AveExpr         t      P.Value    adj.P.Val        B
## SLC6A14    5.024103  9.413107  21.56440 4.104849e-41 8.514279e-37 82.58182
## LOC389023 -3.550396  5.541681 -21.01057 4.054400e-40 4.204818e-36 80.36199
## SLC23A1   -2.473180  5.649224 -17.88487 3.378001e-34 2.335550e-30 67.08748
## DUOX2      4.911030  9.916299  17.37129 3.569259e-33 1.850839e-29 64.78265
## DPP10     -1.910958  3.991413 -16.98863 2.113068e-32 7.304876e-29 63.04259
## TIMP1      2.125930 11.402645  16.88534 3.425860e-32 1.015131e-28 62.56956
##           genesymbol
## SLC6A14      SLC6A14
## LOC389023  LOC389023
## SLC23A1      SLC23A1
## DUOX2          DUOX2
## DPP10          DPP10
## TIMP1          TIMP1
```


## GSEA富集分析

富集分析首选`clusterProfiler`，没有之一！简单，好用！

`clusterProfiler`为我们提供了非常好用的ID转换函数，这里的**ID转换**和上面说的**探针注释**并不是一回事：


```r
library(clusterProfiler)
## 
## clusterProfiler v4.6.2  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
## 
## If you use clusterProfiler in published research, please cite:
## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
## 
## Attaching package: 'clusterProfiler'
## The following object is masked from 'package:stats':
## 
##     filter

gene_entrezid <- bitr(geneID = diff_limma$genesymbol
                         , fromType = "SYMBOL" # 从symbol
                         , toType = "ENTREZID" # 转成ENTREZID
                         , OrgDb = "org.Hs.eg.db"
                         )
## 
## 'select()' returned 1:many mapping between keys and columns
head(gene_entrezid)
##    SYMBOL ENTREZID
## 1 SLC6A14    11254
## 3 SLC23A1     9963
## 4   DUOX2    50506
## 5   DPP10    57628
## 6   TIMP1     7076
## 7    LCN2     3934
```

富集分析最好用`ENTREZID`进行，关于多种不同的ID，在曾老师的书中都有详细介绍，强烈推荐初学者一定要看：[生信初学者基础知识资源推荐](https://mp.weixin.qq.com/s/T-C2xXbpyICC90TgLIJoSQ)

做GSEA分析对数据格式有要求，之前也说过，需要是一个有序的数值型向量，其名字是基因的ID


```r
gene_entrezid <- merge(gene_entrezid,diff_limma,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist,decreasing = T)

head(genelist)
##     4314    11254    50506     1673     1116     6279 
## 5.123666 5.024103 4.911030 4.608619 4.552790 4.256463
```

我们使用`msigdbr`包从msigdb数据库下载人类的`C5`注释集，大家常用的GO、KEGG的数据其实都是包括在msigdb数据库中的。


```r
library(msigdbr)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
## # A tibble: 6 × 2
##   gs_name                                          entrez_gene
##   <chr>                                                  <int>
## 1 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS       60496
## 2 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS       10840
## 3 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS      160428
## 4 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS        4522
## 5 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS       25902
## 6 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS      441024
```

然后是进行`GSEA`分析：


```r
gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 seed = 456
                 )
## preparing geneSet collections...
## GSEA analysis...
## leading edge analysis...
## done...
```

富集分析得到的结果是一个对象，关于这个对象包括那些东西，如何对它进行各种操作，我们在之前的推文都介绍过了，这里就不多说了~

如何查看某个条目下的所有基因名字，很简单，也是不断的取子集操作：


```r
# 第一个条目的所有基因
gsea_res[[gsea_res$ID[[1]]]]
##  [1] "1673"   "5068"   "2919"   "5967"   "5968"   "1670"   "1671"   "2920"  
##  [9] "3553"   "6373"   "4283"   "10563"  "6374"   "6372"   "6283"   "2921"  
## [17] "6280"   "1755"   "3627"   "5266"   "718"    "725"    "4057"   "931"   
## [25] "629"    "3426"   "6278"   "28461"  "1604"   "6347"   "54210"  "3495"  
## [33] "1380"   "1378"   "1191"   "722"    "1880"   "5199"   "5648"   "5788"  
## [41] "1235"   "9308"   "717"    "3569"   "3500"   "710"    "11005"  "5196"  
## [49] "716"    "715"    "1236"   "64127"  "5079"   "940"    "3507"   "28639" 
## [57] "10417"  "2219"   "5452"   "124976" "2213"   "5450"   "10578"  "28912" 
## [65] "3514"   "6590"   "3119"   "6036"   "3853"   "7124"   "4068"   "5919"  
## [73] "4049"   "6406"   "8547"   "3605"   "5473"   "3458"   "729230" "6480"  
## [81] "240"    "3123"   "4239"
```

下面进入今天的正题，可视化！


```r
library(enrichplot)
library(ggplot2)
```

`enrichplot`中包含超多种可视化方法，可以前一篇推文，我们今天主要介绍专门针对`GSEA`结果的可视化。

## 峰峦图

通过函数`ridgeplot`函数实现。展示核心富集基因的表达倍数变化的分布情况。

X轴是核心基因的表达倍数变化，正值表示上调，负值表示下调。


```r
ridgeplot(gsea_res,
          showCategory = 20,
          fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
          )+
  theme(axis.text.y = element_text(size=8))
## Picking joint bandwidth of 0.212
```

![plot of chunk unnamed-chunk-10](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-10-177988517.png)

## gseadist

展示基因集的logFC分布


```r
ids <- gsea_res@result$ID[10:15]

gseadist(gsea_res,
         IDs = ids,
         type="density" # boxplot
         )+
  theme(legend.direction = "vertical")
```

![plot of chunk unnamed-chunk-11](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-177988517.png)

## gsearank

展示基因的排序以及富集分数的变化。


```r
gsearank(gsea_res,
         geneSetID = 1 # 要展示的基因集
         )
```

![plot of chunk unnamed-chunk-12](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-12-177988517.png)

这个函数还可以直接返回画图数据:


```r
aa <- gsearank(gsea_res, 1, title = gsea_res[1, "Description"],output = "table")

head(aa)
##   gene rank in geneList running ES core enrichment
## 1 1673                4 0.02593614             YES
## 2 5068                8 0.04789110             YES
## 3 2919                9 0.06921595             YES
## 4 5967               10 0.09013150             YES
## 5 5968               11 0.11055815             YES
## 6 1670               19 0.12821970             YES
```

## gseaplot

接下来着重介绍`gseaplot`和`gseaplot2`函数。

`gseaplot`函数可以画两个图：`ES`或者`ranked-gene-list`，通过参数`by`设置，默认是两个图都画出来，如果`by="runningScore"`，则是画出`ES`的图，如果是`by = "preranked"`,则是画出ranked gene list的图，


```r
p <- gseaplot(gsea_res, geneSetID = 1, by = "runningScore", 
         title = gsea_res$Description[1])
p
```

![plot of chunk unnamed-chunk-14](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-14-177988517.png)

单独画的时候这两个图都是`ggplot2`对象，可以使用所有`ggplot2`语法修改图形。


```r
p+theme(plot.title = element_text(size = 8,color="red"))
```

![plot of chunk unnamed-chunk-15](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-15-177988517.png)





单独画logfc标准化（也可以选择其他度量方式）之后进行排序的图形：


```r
p <- gseaplot(gsea_res, geneSetID = 1, by = "preranked", 
         title = gsea_res$Description[1])
p
```

![plot of chunk unnamed-chunk-17](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-17-177988517.png)

可以直接使用`ggplot2`语法进行修改：


```r
p+theme(plot.title = element_text(size = 10,color="blue"))
```

![plot of chunk unnamed-chunk-18](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-18-177988517.png)

如果两个子图都画的话返回的是一个`ggplist`对象，此时如果要修改图形细节，可以使用取子集的方法提取其中的子图形，此时的子图形是`ggplot`对象，又可以使用`ggplot2`语法修改了。


```r
#直接加theme返回null,因为是gglist，不是ggplot object
p <- gseaplot(gsea_res,geneSetID = 1,title = gsea_res$Description[1])
p
```

![plot of chunk unnamed-chunk-19](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-19-177988517.png)

此时的标题可能太大了，需要改小一点，可以通过以下方法进行：


```r
#取子集进行修改
p[[1]] <- p[[1]]+theme(plot.title = element_text(size = 6))
p
```

![plot of chunk unnamed-chunk-20](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-20-177988517.png)

所有的细节都支持使用这种方式进行修改！因为它本质上还是`ggplot`对象，只不过是拼图而已。

## gseaplot2

### 基本图形

GSEA富集分析的主要可视化图形还是以下这种，通过`gseaplot2`实现：


```r
# 默认subplots = 1:3，把3个图放一起
gseaplot2(gsea_res,geneSetID = 1,title = "title",
          subplots = 1:3,
          base_size = 10)
```

![plot of chunk unnamed-chunk-21](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-21-177988517.png)

这张图是由3部分组合而成的，3个部分由参数`subplots`控制，所以也不是`ggplot-object`，而是`gglist`，所以如果你要修改其中细节，也是要通过取子集的方法进行。

下面是第一个子图：


```r
gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
```

![plot of chunk unnamed-chunk-22](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-22-177988517.png)

下面是第1和第2个子图拼一起：


```r
gseaplot2(gsea_res, geneSetID = 1, subplots = 1:2)
```

![plot of chunk unnamed-chunk-23](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-23-177988517.png)

第3个子图是什么样的应该就不用演示了~

支持通过提取子图进行自定义修改：


```r
#把entrezid变为symbol
gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db","ENTREZID")

p <- gseaplot2(gsea_res_symbol,geneSetID = 1,
               title = gsea_res_symbol$Description[1])

p[[1]] <- p[[1]]+
  theme(title = element_text(color = "red"))

p
```

![plot of chunk unnamed-chunk-24](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-24-177988517.png)

### 展示多条通路

可以同时画多条通路。


```r
tmp <- as.data.frame(gsea_res_symbol)
colnames(tmp)
##  [1] "ID"              "Description"     "setSize"         "enrichmentScore"
##  [5] "NES"             "pvalue"          "p.adjust"        "qvalue"         
##  [9] "rank"            "leading_edge"    "core_enrichment"
head(tmp,2)
##                                                          ID
## GOBP_HUMORAL_IMMUNE_RESPONSE   GOBP_HUMORAL_IMMUNE_RESPONSE
## GOBP_ADAPTIVE_IMMUNE_RESPONSE GOBP_ADAPTIVE_IMMUNE_RESPONSE
##                                                 Description setSize
## GOBP_HUMORAL_IMMUNE_RESPONSE   GOBP_HUMORAL_IMMUNE_RESPONSE     226
## GOBP_ADAPTIVE_IMMUNE_RESPONSE GOBP_ADAPTIVE_IMMUNE_RESPONSE     431
##                               enrichmentScore      NES pvalue     p.adjust
## GOBP_HUMORAL_IMMUNE_RESPONSE        0.7441821 2.629284  1e-10 7.277857e-09
## GOBP_ADAPTIVE_IMMUNE_RESPONSE       0.6841391 2.548253  1e-10 7.277857e-09
##                                     qvalue rank                   leading_edge
## GOBP_HUMORAL_IMMUNE_RESPONSE  5.821805e-09 1251  tags=37%, list=8%, signal=34%
## GOBP_ADAPTIVE_IMMUNE_RESPONSE 5.821805e-09 2125 tags=46%, list=13%, signal=41%
```

比如我们通过查找发现自己想要把第4,5,6个通路画在一起，只要给`geneSetID`提供参数即可：


```r
#提供多个通路
p <- gseaplot2(gsea_res,geneSetID = 4:6)
p
```

![plot of chunk unnamed-chunk-26](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-26-177988517.png)

现在这个通路名字很长，显示不全，如果你对`ggplot`很熟悉你就能知道上面这张图的通路名字在这里很明显就是一个图例，所以我们可以取子图，然后对子图进行操作，把图例放到图形上方：


```r
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
```

![plot of chunk unnamed-chunk-27](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-27-177988517.png)

当然你也可以通过修改图例标签实现，详情请见：[ggplot2修改图例详细解读](https://mp.weixin.qq.com/s/OGiOV2ve43gr1G0J-ahYzQ)


```r
p <- gseaplot2(gsea_res,geneSetID = 4:6)
p[[1]] <- p[[1]]+scale_color_hue(labels=c("aa","bb","cc"))
p
```

![plot of chunk unnamed-chunk-28](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-28-177988517.png)

当然一些比较简单的参数y叔已经给你准备好了，比如`base_size`控制整体字体大小，`color`改变颜色映射。


```r
p <- gseaplot2(gsea_res,geneSetID = 4:6,
               base_size = 10,
               color = c("#E495A5", "#86B875", "#7DB0DD")
               )
p
```

![plot of chunk unnamed-chunk-29](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-29-177988517.png)

但是如果你需要个性化的出图或者为了发文章，那肯定是需要自己DIY一番的喽~


```r
p <- gseaplot2(gsea_res,geneSetID = 4:6)
p[[1]] <- p[[1]]+scale_color_viridis_d(labels=c("lalala","heiheihei","dadada"))+
  geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)+
  theme(legend.position = "top")
p[[2]] <- p[[2]]+scale_color_viridis_d()
p[[3]] <- p[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
p
```

![plot of chunk unnamed-chunk-30](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-30-177988517.png)

所有的一切都可以自定义！这么多`ggplot2`的语法我是怎么知道的呢？很简单，两本`ggplot2`说明书：《R数据可视化手册》和《ggplot2数据分析与图形艺术》，买一本认真看一遍你就懂了！

### 展示P值信息

除此之外还可以显示`pvalue`信息，但是很遗憾不能显示`NES`~


```r
#有些通路名字很长，表格会显示不出来，注意调整图形宽度
p <- gseaplot2(gsea_res, geneSetID = 4:6, 
          pvalue_table = TRUE # 显示pvalue
          )
p
```

![plot of chunk unnamed-chunk-31](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-31-177988517.png)

你会发现通路名字不是很长，表格有时会显示不出来，但是此时通路名字变成了`annotate`而不是图例，所以就不能像上面修改图例那样修改这里的通路名字了！


```r
# 这段代码没啥用了~
p[[1]] <- p[[1]]+
  theme(legend.position = "top",
        legend.direction = "vertical"
        )
p
```

![plot of chunk unnamed-chunk-32](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-32-177988517.png)

这里的通路名字并不是图例，所以上面的代码并没有达到把通路名字移走的效果~

这种情况下通路名字和`p_value table`其实是`annotate`。

`gseaplot2`函数中对这个`annotate`的处理如下代码所示：


```r
if (pvalue_table) {
        pd <- x[geneSetID, c("Description","pvalue", "p.adjust")]
        # pd <- pd[order(pd[,1], decreasing=FALSE),]
        rownames(pd) <- pd$Description

        pd <- pd[,-1]
        # pd <- round(pd, 4)
        for (i in seq_len(ncol(pd))) {
            pd[, i] <- format(pd[, i], digits = 4)
        }
        tp <- tableGrob2(pd, p.res)

        p.res <- p.res + theme(legend.position = "none") +
            annotation_custom(tp,
                              xmin = quantile(p.res$data$x, .5),
                              xmax = quantile(p.res$data$x, .95),
                              ymin = quantile(p.res$data$runningScore, .75),
                              ymax = quantile(p.res$data$runningScore, .9))
    }
```

那这样的话其实我们可以自己添加，绝对更加个性化，想怎么弄都行。

如果你一定要用`ggplot2`的默认颜色，可以通过以下方式获取，`scales`作为`ggplto2`扩展包，功能十分实用，我们之前也详细介绍过：[实用R包scales包介绍](https://mp.weixin.qq.com/s/MS-c_tXoBwaF6yu3_62R4g)


```r
library(scales)
hex <- hue_pal()(3)
hex
## [1] "#F8766D" "#00BA38" "#619CFF"
```

这个`annotate`其实是一个`gtable`对象，我们可以通过`gridExtra`包实现对它的精细化控制，然后把它加到图形中即可。


```r
library(gridExtra)

# 选择4,5,6条通路
x <- gsea_res_symbol
geneSetID <- 4:6

# 提取NES，P值等信息
pd <- x[geneSetID, c( "NES","pvalue", "p.adjust")]
pd <- pd[order(rownames(pd), decreasing=FALSE),]
for (i in seq_len(ncol(pd))) {pd[, i] <- format(pd[, i], digits = 4)}

# 通过修改table的主题来修改表格细节
tt <- ttheme_minimal(base_size = 10,
  core=list(#bg_params = list(fill = NA, col=NA),
            fg_params=list(col=c("#F8766D","#00BA38","#619CFF"))
            )
  )

tp <- tableGrob(pd,rows = NULL,theme = tt)

# 修改表格每个格子的宽度和高度
#tp$widths <- unit(rep(1.2,ncol(tp)), "cm")
tp$heights <- unit(rep(0.4,nrow(tp)),"cm") # cell height
```

这样一个表格就画好了，可以画出来看一下先：


```r
plot(tp)
```

![plot of chunk unnamed-chunk-36](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-36-177988517.png)

就是我需要的效果！

有了这个东西，我们取子集，然后用`ggplot2`语法DIY即可：


```r
p <- gseaplot2(gsea_res, geneSetID = 4:6)
p[[1]] <- p[[1]]+
  annotation_custom(tp,
                    xmin = 10000,
                    xmax = 14000,
                    ymin = 0.4,
                    ymax = 0.8
                    )+
  theme(plot.title = element_text(size = 5),
        legend.position = "top",
        legend.direction = "vertical"
        )
p
```

![plot of chunk unnamed-chunk-37](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-37-177988517.png)

很好，就是我需要的最终效果，其他颜色、大小这些就不再演示了~

除此之外，`ggpp`包也能达到类似的效果，下面是一个简单的演示：


```r
x <- gsea_res_symbol
geneSetID <- 4:6

pd <- x[geneSetID, c( "NES","pvalue", "p.adjust")]
rownames(pd) <- NULL
for (i in seq_len(ncol(pd))){pd[, i] <- format(pd[, i], digits = 4)}
```

通过`ggpp`中的`annotate`函数实现：


```r
library(ggpp)
## 
## Attaching package: 'ggpp'
## The following object is masked from 'package:ggplot2':
## 
##     annotate

p <- gseaplot2(gsea_res, geneSetID = 4:6)
p[[1]] <- p[[1]]+
  annotate("table", x = 12000, y = 0.4, label = pd,
           size=3,table.theme = ttheme_gtminimal
           )+
  theme(plot.title = element_text(size = 5),
        legend.position = "top",
        legend.direction = "vertical"
        )
p
```

![plot of chunk unnamed-chunk-39](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-39-177988517.png)

但是这样的问题是并不能知道哪条通路的值，因为全是黑色，这个包还有一个`geom_table`函数，理论上可以为不同的行映射不同的颜色，大家可以自己探索下~

### 展示基因名字

在通路中添加想要展示的基因呢？通过`geom_gsea_gene`函数即可。


```r
#选择基因
symbol <- gsea_res_symbol[[gsea_res_symbol$ID[[1]]]]
head(symbol)
## [1] "DEFB4A" "REG3A"  "CXCL1"  "REG1A"  "REG1B"  "DEFA5"
length(symbol)
## [1] 83

# 随便选5个
g <- sample(symbol,5)
g
## [1] "C1R"   "C4BPB" "IL6"   "LTA"   "CCL2"
```

添加到图形中即可：


```r
p <- gseaplot(gsea_res_symbol, 1, by='runningScore') 

p+geom_gsea_gene(g)
```

![plot of chunk unnamed-chunk-41](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-41-177988517.png)

还是可以通过取子集的方式修改其中的子图形：


```r
p <- gseaplot(gsea_res_symbol, 1) 
p
```

![plot of chunk unnamed-chunk-42](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-42-177988517.png)

```r
p[[2]] <- p[[2]]+geom_gsea_gene(g, geom = ggplot2::geom_label)
p
```

![plot of chunk unnamed-chunk-42](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-42-277988517.png)

`geom_gsea_gene`函数当然也是支持`gseaplot2`的：


```r
library(ggrepel)

p <- gseaplot2(gsea_res_symbol, geneSetID = 6)
p[[1]] <- p[[1]]+geom_gsea_gene(g, geom=geom_text_repel)+
  theme(legend.position = "top",
        legend.direction = "vertical"
        )
p
```

![plot of chunk unnamed-chunk-43](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-43-177988517.png)

多条通路当然也没有问题：


```r
# 每条通路随便展示3个基因
g <- sample(symbol,3)

p <- gseaplot2(gsea_res_symbol, geneSetID = 4:6)
p[[1]] <- p[[1]]+geom_gsea_gene(mapping = aes(color= Description), 
                                g,
                                geom=geom_text_repel)+
  theme(legend.position = "top",
        legend.direction = "vertical"
        )
p
```

![plot of chunk unnamed-chunk-44](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-44-177988517.png)

可以看到有些基因是重复的，因为所有通路都是展示的这几个基因，可以通过分别提取子集再添加基因的方式展示不同通路中的不同基因：


```r
g11 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[1]]],5)
g22 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[2]]],5)
g33 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[3]]],5)
desc <- gsea_res_symbol$Description[1:3]

p <- gseaplot2(gsea_res_symbol, geneSetID = 1:3)

p[[1]] <- p[[1]]  + 
    geom_gsea_gene(mapping=aes(colour = Description), g11, geom=geom_text_repel, geneSet=desc[1]) + 
    geom_gsea_gene(mapping=aes(colour = Description), g22, geom=geom_text_repel, geneSet=desc[2]) +
    geom_gsea_gene(mapping=aes(colour = Description), g33, geom=geom_text_repel, geneSet=desc[3])
p
```

![plot of chunk unnamed-chunk-45](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-45-177988517.png)

今天的内容就到这里，后面会继续给大家介绍富集分析可视化。

