上一次推文已经给大家介绍了常见的富集分析类型以及如何使用全能的R包`clusterprofiler`实现，详情请见：[富集分析常见类型](https://mp.weixin.qq.com/s/RtF7DPXYaObiDauIQTnkFg)

相必大家也知道`clusterprofiler`的富集结果可以画出很多漂亮的图，其中`enrichplot`是专门用于对接`clusterprofiler`的富集结果的可视化R包。

但是大家有没有这样的困惑呢？

- `ehrichplot`可以画哪些图形？
- 每个绘图函数可以对接哪些富集分析的结果？
- 其中一些参数怎么用？

今天就先给大家详细介绍`enrichplot`包。

当然了，富集分析的可视化还有许多其他R包，我们会在后续的推文中继续介绍。

本期目录：

[toc]

## 准备数据

用`gse87466`这个GEO的数据做演示，下载整理的过程这次就不演示了。数据可在我们的粉丝QQ群免费下载。


```r
load(file = "G:/easyTCGA_test/gse87466.Rdata")
```

这是一个炎症性肠病的数据集，一共108个样本，21个normal，87个uc（ulcerative colitis）。


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
                          , is_count = F # 不是count数据
                          , logFC_cut = 0 # 可以直接筛选结果
                          , pvalue_cut = 1 
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

因为接下来会同时演示`ORA`和`GSEA`两种富集分析，所以我们把筛选后的差异基因用于`ORA`分析，所有的基因用于`GSEA`分析。

选取`logFC > 1` & `adj.P.Val<0.01` 的基因作为差异基因进行后续的`ORA`分析：


```r
deg_genes <- diff_limma[abs(diff_limma$logFC)>1 & diff_limma$adj.P.Val<0.01,]
deg_genes <- deg_genes$genesymbol


length(deg_genes)
## [1] 1192
head(deg_genes)
## [1] "SLC6A14"   "LOC389023" "SLC23A1"   "DUOX2"     "DPP10"     "TIMP1"
```

1192个差异基因等下用于`ORA`富集分析。

然后准备下`GSEA`需要的格式。

富集分析最好用`ENTREZID`进行，关于多种不同的ID，在曾老师的书中都有详细介绍，强烈推荐初学者一定要看：[生信初学者基础知识资源推荐](https://mp.weixin.qq.com/s/T-C2xXbpyICC90TgLIJoSQ)


```r
suppressMessages(library(clusterProfiler))

gene_entrezid <- bitr(geneID = diff_limma$genesymbol
                         , fromType = "SYMBOL" # 从symbol
                         , toType = "ENTREZID" # 转成ENTREZID
                         , OrgDb = "org.Hs.eg.db"
                         )
## 
## 'select()' returned 1:many mapping between keys and columns

gene_entrezid <- merge(gene_entrezid,diff_limma,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist,decreasing = T)

head(genelist)
##     4314    11254    50506     1673     1116     6279 
## 5.123666 5.024103 4.911030 4.608619 4.552790 4.256463
```

这样GSEA需要的数据也准备好了。


## 富集分析

富集分析首选`clusterProfiler`，没有之一！简单，好用！

富集分析最好用`ENTREZID`进行，但其实不转换也可以进行，富集分析时会给你转换，你只要指定类型即可，这里是因为`enrichGO`富集分析会借助`Org`注释包进行，里面含有多种不同的基因ID，它可以自动帮你进行转换，如果没有使用`Org`注释包的富集分析函数就只能用`ENTREZID`。

首先进行`ORA`:


```r
ora_res <- enrichGO(gene = deg_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "SYMBOL",#这里指定ID类型
                   ont = "ALL", # "BP", "MF", "CC" 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   minGSSize = 10,# 最少的基因数量
                   maxGSSize = 500, # 最大的基因数量
                   readable = T # 把ENTREZID转换为SYMBOL
                   )

head(ora_res,3)
##            ONTOLOGY         ID                              Description
## GO:0050900       BP GO:0050900                      leukocyte migration
## GO:0097530       BP GO:0097530                    granulocyte migration
## GO:0002237       BP GO:0002237 response to molecule of bacterial origin
##            GeneRatio   BgRatio       pvalue     p.adjust       qvalue
## GO:0050900   84/1002 398/18903 2.247121e-28 1.191424e-24 8.456271e-25
## GO:0097530   52/1002 158/18903 1.187840e-27 3.148963e-24 2.235014e-24
## GO:0002237   77/1002 360/18903 1.792733e-26 2.871742e-23 2.038253e-23

class(ora_res)
## [1] "enrichResult"
## attr(,"package")
## [1] "DOSE"
```

这个结果是一个`enrichResult`对象，

下面进行`GSEA`富集分析：


```r
gsea_res <- gseGO(gene = genelist,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500
                   )
## preparing geneSet collections...
## GSEA analysis...
## leading edge analysis...
## done...

head(gsea_res,3)
##            ONTOLOGY         ID
## GO:0006959       BP GO:0006959
## GO:0002250       BP GO:0002250
## GO:0002460       BP GO:0002460
##                                                                                                                          Description
## GO:0006959                                                                                                   humoral immune response
## GO:0002250                                                                                                  adaptive immune response
## GO:0002460 adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
##            setSize enrichmentScore      NES pvalue     p.adjust       qvalue
## GO:0006959     227       0.7406039 2.653306  1e-10 4.840719e-09 3.654586e-09
## GO:0002250     439       0.6876408 2.594300  1e-10 4.840719e-09 3.654586e-09
## GO:0002460     286       0.6947694 2.536881  1e-10 4.840719e-09 3.654586e-09
##            rank                   leading_edge
## GO:0006959 1251  tags=36%, list=8%, signal=33%
## GO:0002250 2125 tags=47%, list=13%, signal=42%
## GO:0002460 2243 tags=49%, list=14%, signal=43%

class(gsea_res)
## [1] "gseaResult"
## attr(,"package")
## [1] "DOSE"
```

这个结果是`gseaResult`对象。

有了这两个结果，我们就可以演示`enrichplot`的用法了。

不知道大家注意到没有，`clusterProfiler`主要就是进行`ORA`和`GSEA`富集分析，

- 如果是进行`ORA`，那么结果就是`enrichResult`对象，
- 如果是进行`GSEA`，那么结果就是`gseaResult`对象。
- （还有一个例外，`compareCluster`，它的结果是`compareClusterResult`对象）

`enrichplot`可以专门对接`clusterProfiler`的输出结果，目前一共包括十几个主要的绘图函数，主要就是针对`enrichResult`和`enrichResult`的(还有`compareClusterResult`)，不过并不是每种图形都可以同时对接多个对象。

我们今天主要以`ORA`的富集结果可视化为例，因为对于同时支持`enrichResult`和`enrichResult`多种结果的函数来说，它们的用法基本一模一样！并且一些函数是专门针对`GSEA`结果的，我们会在后面的推文中继续介绍。

`enrichplot`是基于`ggplot2`的，所以所有的`ggplot2`特性都是支持的。

## 条形图

通过`barplot`实现，此函数只能对接`enrichResult`对象，所以`GSEA`的结果它是画不出来的哦~

用于展示最重要的或者你感兴趣的条目的富集结果，比如富集到的基因个数、条目名字，P值等信息。


```r
library(enrichplot)

p1 <- barplot(ora_res, 
        showCategory=10 # 展示多少条目
        ,x = "Count" # X轴展示那个变量，默认Count,也可以是GeneRatio
        ,label_format = 30 # 默认对名字超过30个字符的进行折叠
        ,font.size = 12 # 字体大小
        ,title = "Bar plot for ORA"
        ) 
p2 <- barplot(ora_res, 
        showCategory=10 
        ,x = "GeneRatio" 
        ) 

cowplot::plot_grid(p1,p2)
```

![plot of chunk unnamed-chunk-8](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-8-177986737.png)

我们在做GO的ORA时，同时进行了"BP", "MF", "CC"三种类别的分析，所以可以根据这个变量进行分面展示。


```r
barplot(ora_res, showCategory=10
        ,split = "ONTOLOGY" # 分面，GO ORA特有
        ) +
  facet_grid(ONTOLOGY~., scale="free") # ggplot2的分面语法
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230624170453848.png)

`ggplot2`可以非常方便的进行各种自定义可视化，可以参考历史推文：

- [ggplot2坐标轴修改详细教程](https://mp.weixin.qq.com/s/VhraGHo-HQ-mmcIwNn7XMQ)
- [ggplot2图例修改详细教程](https://mp.weixin.qq.com/s/OGiOV2ve43gr1G0J-ahYzQ)
- [ggplot2分面图形大改造](https://mp.weixin.qq.com/s/kgFy52uxM_z_vI0C6s2Rgg)

因为`barplot`函数是对接`enrichResult`对象的，所以其他类型的`ORA`分析也是可以直接用的，详情可参考上一篇推文：[富集分析常见类型](https://mp.weixin.qq.com/s/RtF7DPXYaObiDauIQTnkFg)

GSEA的条形图怎么办呢？提取数据自己画即可。

## gene-concept网络图

通过函数`cnetplot`实现，这个函数可以同时对接`enrichResult`、`gseaResult`、`compareClusterResult`3种结果，非常强大。

网络图可以展示不同的条目富集到了哪些基因，还可以用连接线连起来，内容比起条形图和气泡图更加丰富。

首先看下默认的出图及参数含义，有一些参数可能在接下来的版本中被移除(因为有了更灵活的其他参数可以替代)，所以这里就不介绍了。


```r
# 默认
p1 <- cnetplot(ora_res
         , showCategory = 5 # 也可以直接写条目名字
         , layout = "kk" #网络形状，’star’, ’circle’, ’gem’, ’dh’, ’graphopt’, ’grid’, ’mds’, ’randomly’, ’fr’, ’kk’, ’drl’ or ’lgl’
         , node_label = "all" #显示哪些节点的标签，’category’, ’gene’, ’all’(默认), ’none’
         , shadowtext = "all" #哪些节点标签需要添加阴影，’category’, ’gene’, ’all’(默认), ’none’
         # 控制节点和连线的颜色
         , color.params = list(foldChange = NULL 
                               , edge = FALSE #根据富集的不同条目上色
                               , category = "#E5C494" #条目节点颜色
                               , gene = "#B3B3B3" #基因节点颜色
                               )
         
         # 控制标签和节点的大小
         , cex.params = list(category_node = 1
                             , gene_node = 1
                             , category_label = 1
                             , gene_label = 1)
         
         # 控制哪些节点和连线高亮显示
         , hilight.params = list(category = NULL
                                 , alpha_hilight = 1
                                 , alpha_no_hilight = 0.3)

         )
p1
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230624170629146.png)

下面我们看看不同的形状，并且只显示条目的名字，按照基因变化倍数给基因上色，按照富集到的不同条目给连线上色：

注意这里ORA的`foldChange`参数需要格式，和做`GSEA`的格式一样，不过不需要排序也可以。


```r
#对于ORA需要自己准备下 foldchange
foldChange <- diff_limma$logFC
names(foldChange) <- diff_limma$genesymbol

p2 <- cnetplot(ora_res, showCategory = 3,
               layout = "star", node_label = "category"
               , color.params = list(foldChange = foldChange #把基因的倍数变化映射给基因节点的颜色
                               , edge = T 
                               , category = "red" 
                               #, gene = "yellow" # 因为设置了foldChange，所以这里没用了
                               )
               )
## Scale for size is already present.
## Adding another scale for size, which will replace the existing scale.

p2
```

![plot of chunk unnamed-chunk-11](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-177986737.png)

上图因为设置了foldChange，所以不能统一基因颜色了，这里演示下基因颜色怎么改成统一的，并演示下如何控制条目和基因标签的大小：


```r
p3 <- cnetplot(ora_res, 
               showCategory = 2, layout = "circle"
               , color.params = list(gene = "yellow" # 基因节点的颜色
                                     )
               # 控制标签和节点的大小
               , cex.params = list(category_node = 2
                             , gene_node = 0.5
                             , category_label = 2
                             , gene_label = 0.5)
               )

p3
```

![plot of chunk unnamed-chunk-12](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-12-177986737.png)

- `showCategory`参数可以接受条目名字组成的字符串向量，
- 高亮显示1个或者多个条目及基因，


```r
cate <- c("nucleotide receptor activity","peptide binding")

p4 <- cnetplot(ora_res, showCategory = cate, layout="graphopt"
               # 控制哪些节点和连线高亮显示
               , hilight.params = list(category = "nucleotide receptor activity"
                                 , alpha_hilight = 0.8 # 高亮显示的颜色深度
                                 , alpha_no_hilight = 0.3)
               )

p4
```

![plot of chunk unnamed-chunk-13](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-13-177986737.png)

对于GSEA的富集分析也是一样的支持：


```r
# 把ENTREZID改为symbol，方便演示
gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db",keyType = "ENTREZID")

cnetplot(gsea_res_symbol, showCategory=3
         , color.params = list(foldChange = genelist)#这里可以直接用genelist，因为格式符合要求
         )
## Scale for size is already present.
## Adding another scale for size, which will replace the existing scale.
```

![plot of chunk unnamed-chunk-14](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-14-177986737.png)

由于`showCategory`可以使用字符创型向量，富集分析的结果也是非常简单地可以进行取子集操作，所以我们也可以单独展示上下调的通路。


```r
names(gsea_res@result)
##  [1] "ONTOLOGY"        "ID"              "Description"     "setSize"        
##  [5] "enrichmentScore" "NES"             "pvalue"          "p.adjust"       
##  [9] "qvalue"          "rank"            "leading_edge"    "core_enrichment"
```

对结果取子集，上一篇推文介绍过的，非常简单的`dplyr`数据操作：


```r
#先把entrezid变成symbol
gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db")

#选择上调的通路，这里选了前6个
up <- gsea_res_symbol %>% 
  arrange(desc(NES)) %>% 
  slice(1:6)

up <- up$Description
```


```r
cnetplot(gsea_res_symbol, showCategory=up,
         color.params = list(foldChange = genelist)
         )
## Scale for size is already present.
## Adding another scale for size, which will replace the existing scale.
```

![plot of chunk unnamed-chunk-17](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-17-177986737.png)

除此之外，`cnetplot`还支持一个`list`类型，在我的探索下，发现它的作用是可视化两个富集结果组成的列表，列表内必须指定元素名字：


```r
gg1 <- deg_genes[1:500]
gg2 <- deg_genes[500:1000]

# 做两次富集分析
ggg1 <- enrichGO(gg1,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL")
ggg2 <- enrichGO(gg2,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL")

# 结果放进list，并给定名字
cnetplot(list(aa = ggg1, bb = ggg2)
         #,node_label = "category"
         )
```

![plot of chunk unnamed-chunk-18](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-18-177986737.png)

此时子节点是条目名字！那么此时上面那么多控制参数还有哪些可用呢？留给需要的人自己探索~

## 气泡图

通过函数`dotplot`实现，和`barplot`函数很像，只不过是增加了点的大小这个映射，可以多展示一列变量。

同时支持`enrichResult`、`gseaResult`、`compareClusterResult`3种结果。


```r
dotplot(ora_res, showCategory=20)
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230624182347850.png)

分面：


```r
dotplot(ora_res, showCategory=10
        ,split = "ONTOLOGY"
        ) +
  facet_grid(ONTOLOGY~., scale="free")
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230624182415441.png)

和`barplot`用法基本一样，我们就不再过多演示了。

## emapplot网络图

`emapplot`也是一种网络图，不过可以把相似的条目聚集到一起，便于识别不同的功能模块。

这个函数同样也是可以直接使用`enrichResult`、`gseaResult`、`compareClusterResult`3种结果。

不过在使用前，必须用`pairwise_termsim`函数添加相似性矩阵才行~


```r
ora_pt <- pairwise_termsim(ora_res)

#默认参数画图
emapplot(ora_pt
         ,showCategory = 30
         ,color = "p.adjust" # 映射给条目颜色 ’pvalue’, ’p.adjust’ or ’qvalue’
         ,shadowtext = TRUE # 显示标签阴影
         ,repel = FALSE # 解决标签重叠问题
         ,node_label = "category" #显示谁的标签  ’category’,’group’,’all', ’none’
         
         #形状控制参数
         ,layout.params = list(layout = NULL #和cnetplot的形状参数一样，后面也是，就不多说了
                               ,coords = NULL #控制位置，需要含2列的data.frame,x是x轴坐标,y是y轴坐标
                               )
         
         #控制连线
         ,edge.params = list(show = TRUE # 是否显示连线
                             , min = 0.2) #判断两个节点是否相似的阈值
         
         #大小控制参数
         ,cex.params = list(category_node = 1 # 节点大小
                            , category_label = 1 #节点标签大小
                            , line = 1 # 线的粗细
                            #, pie2axis #饼图大小
                            #, label_group #分组标签大小
                            )
         
         #控制高亮，和cnetplot同，不再多说
         ,hilight.params = list(category = NULL
                                ,alpha_hilight = 1
                                ,alpha_no_hilight = 0.3)
         
         #控制聚类
         ,cluster.params = list(cluster = FALSE #对条目聚类
                                ,method = stats::kmeans #聚类方法
                                ,n = NULL #聚类个数
                                , legend = FALSE #显示聚类图例
                                ,label_style = "shadowtext"
                                ,label_words_n = 4 #聚类标签个数
                                ,label_format = 30 #最长的字符数，控制折叠
                                )
         )
```

![plot of chunk unnamed-chunk-21](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-21-177986737.png)

由于形状、大小、高亮我们在上面的`cnetplot`中演示过了，所以这里就不再演示了~

简单展示下聚类相关的参数。


```r
library(ggplot2)

emapplot(ora_pt,
         #node_label = "none",
         cluster.params = list(cluster = T
                                ,method = stats::kmeans 
                                ,n = 5 
                                , legend = T 
                                ,label_style = "shadowtext"
                                ,label_words_n = 4 
                                ,label_format = 30 
                                )
         )+
  theme(legend.position = "top",
        legend.direction = "vertical"
        )
```

![plot of chunk unnamed-chunk-22](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-22-177986737.png)

这张图就很有意思了，如果你用过`cytoscape`做过富集分析，那你一定对这种样式很熟悉。它给你增加了置信椭圆，一眼就能看出几个模块。

## 有向无环图

通过`goplot`实现，支持`enrichResult`、`gseaResult`

注意在做富集分析时要指定`ont = "BP"`，不然画不出来~


```r
ora_bp <- enrichGO(gene = deg_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "SYMBOL",
                   ont = "BP")

# 去除冗余的结果
ora_sim <- simplify(ora_bp)
```

画图：


```r
#ora_bp也可以
goplot(ora_sim
       ,showCategory = 10
       ,color = "p.adjust"
       ,layout = "sugiyama"
       ,geom = "text" #f ’label’ or ’text’
       )
```

![plot of chunk unnamed-chunk-24](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-24-177986737.png)


## 热图

热图展示的信息和`cnetplot`展示的信息差不多，通过函数`heatplot`实现，同时支持`enrichResult`和`gseaResult`。

不过参数就比`cnetplot`函数少的多了，并且基本上都在上面的函数中介绍过了，这里也就不再重复解释了。


```r
#对于ORA需要自己准备下 foldchange
foldChange <- diff_limma$logFC
names(foldChange) <- diff_limma$genesymbol

#自己准备下pvalue，也是一样的格式
pval <- diff_limma$adj.P.Val
names(pval) <- diff_limma$genesymbol

#选择想要展示的条目，适合基因数目比较少的
cate <- c("nucleotide receptor activity",
          "organic cation transmembrane transporter activity",
          "polyol transmembrane transporter activity",
          "purinergic nucleotide receptor activity",
          "steroid dehydrogenase activity",
          "oxidoreductase activity, acting on the CH-NH2 group of donors",
          "glycosphingolipid binding","organic acid:sodium symporter activity",
          "prostanoid receptor activity")

#画图即可
heatplot(ora_res
         ,showCategory = cate
         ,symbol = "dot" # rect dot
         ,foldChange = foldChange
         ,pvalue = pval
         ,label_format = 30
         )
```

![plot of chunk unnamed-chunk-25](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-25-177986737.png)

## ssplot降维图形

通过函数`ssplot`实现，支持`enrichResult`、`gseaResult`、`compareClusterResult`3种结果。

可以把富集分析结果进行降维，画出降维之后的图形。

这个函数也有非常多参数，其中一些是控制大小、控制标签等，和上面介绍的函数中即将被移除的参数名字一样，会在接下来的版本中被移除，所以就不多介绍了。


```r
#需要计算相似性
ora_pt <- pairwise_termsim(ora_res)

ssplot(ora_pt
       ,showCategory=30
       ,drfun = NULL #降维用的算法 stats::cmdscale (the default),vegan::metaMDS, or ape::pcoa
       )
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230622150759895.png)

## 聚类树

把富集结果通过聚类树的性质展示，使用函数`treeplot`实现，支持`enrichResult`、`gseaResult`、`compareClusterResult`3种结果.

这个图形展示的信息有点类似`emapplot`网络图，可以通过**聚类**发现功能近似的条目，有利于发现功能模块。

很多参数已经在上面几个函数介绍过了，这里还是简单一点说。

不知道大家发现了吗，可以聚类、计算距离的函数都是需要提前用`pairwise_termsim`处理一下的哦。


```r
library(ggplot2)

#计算相似性
go_pt <- pairwise_termsim(ora_res)

treeplot(go_pt
         ,showCategory = 30
         ,color = "p.adjust" #pvalue, p.adjust, qvalue
         ,label_format = NULL
         ,fontsize = 4
         
         #添加颜色框
         ,hilight.params = list(hilight = T #是否添加颜色框
                                , align = "both" # 颜色框的对齐方式’none’,’left’,’right’,’both’
                                )
         
         #控制各个元素之间的距离
         ,offset.params = list(bar_tree = rel(1)
                               ,tiplab = rel(1)
                               ,extend = 0.3
                               ,hexpand = 0.1
                               )
         
         #聚类相关参数，见前面的解释
         ,cluster.params = list(method = "ward.D"
                                , n = 5
                                , color =  NULL #控制颜色框的颜色，c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
                                , label_words_n = 4
                                , label_format = 30
                                ),
         )
```

![plot of chunk unnamed-chunk-27](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-27-177986737.png)

## upsetplot

通过`upsetplot`的形式展示多个条目之间的交集情况，可以看出哪些基因富集在哪些条目中。

通过函数`upsetplot`函数实现，同时支持`enrichResult`、`gseaResult`两种结果。

这个函数的参数比较少，对于ORA结果，会展示不同条目之间重叠的基因：


```r
upsetplot(ora_res, n=10)
```

![plot of chunk unnamed-chunk-28](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-28-177986737.png)

对于GSEA结果，除了展示不同条目间基因的交集外，还可以通过箱线图展示不同条目的倍数变化（foldchange）：

我们这里对`GSEA`的结果取个子集，因为原结果我不喜欢，取子集的方法也在之前的推文介绍过了：

- [富集分析常见类型]()


```r
gsea_sub <- gsea_res[gsea_res$setSize<200,asis=T]

#通路名字太长了，没有折叠
upsetplot(gsea_sub,n=10) 
```

![plot of chunk unnamed-chunk-29](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-29-177986737.png)

关于`upsetplot`，在之前的历史推文中介绍非常多，可以参考：

- [韦恩图进阶之UpSetR（1）](https://mp.weixin.qq.com/s/Ikwp4uHtdXy-KNSlpvNgQg)
- [韦恩图进阶之UpSetR（2）](https://mp.weixin.qq.com/s/W2kRqCVxtlPc-JWJU-0ePA)
- [韦恩图进阶之UpSetR（3）](https://mp.weixin.qq.com/s/kT4igXA9_nH4-8QTLyTEiA)
- [韦恩图进阶之UpSetR（4）](https://mp.weixin.qq.com/s/lTGLtGvOkP7l3c7FiWujHg)
- [韦恩图进阶之ComplexHeatmap](https://mp.weixin.qq.com/s/CI-wTadPj2lLuEGM_9tTYA)
- [韦恩图进阶之ggupset](https://mp.weixin.qq.com/s/czbGBdCArjCm9ouuHTJbJw)
- [最强大的upset plot 之ComplexUpset（上）](https://mp.weixin.qq.com/s/aisqPyMZ4-UDpPqyDfohTw)
- [最强大的upset plot 之ComplexUpset（下）](https://mp.weixin.qq.com/s/XYIs52AXBMaARhYBhRm2dw)

---

今天的介绍就到这里，主要是以`ORA`结果为例，展示了`enrichplot`强大的绘图功能，无缝对接`clusterprofiler`的富集结果。

但是关于`compareClusterResult`的结果也没有介绍，因为有一些参数是专门针对它的，所以放到后面再介绍。



