`TCGAbiolinks`是一个大而全的R包，常见的分析都能做，比如差异分析、富集分析、生存分析等等。上次学习了差异分析，今天学习下富集分析和生存分析。

但其实这些分析我们都是把数据保存好用其他包做的~

## 富集分析

如果大家要做富集分析，那肯定是首选`clusterprofiler`了，不过今天这个还是要学习下。

还是用上篇推文得到的`coadDEGs`继续进行富集分析。

[新版TCGAbiolinks包学习：差异分析](https://mp.weixin.qq.com/s/0SLQOZRkZ4hOQY1ETnQRUA)


```r
load(file = "coadDEGs.Rdata")
```

在`TCGAbiolinks`里进行富集分析很简单，就一句代码搞定。


```r
library(TCGAbiolinks)

Genelist <- coadDEGs$gene_name # gene_symbol

# 进行GO和KEGG分析
ansEA <- TCGAanalyze_EAcomplete(
    TFname = "TCGAbiolinks enrichment analysis",
    RegulonList = Genelist
)
## [1] "I need about  1 minute to finish complete  Enrichment analysis GO[BP,MF,CC] and Pathways... "
## [1] "GO Enrichment Analysis BP completed....done"
## [1] "GO Enrichment Analysis MF completed....done"
## [1] "GO Enrichment Analysis CC completed....done"
## [1] "Pathway Enrichment Analysis completed....done"


# 富集分析结果可视化
TCGAvisualize_EAbarplot(
    tf = rownames(ansEA$ResBP), 
    GOBPTab = ansEA$ResBP,
    GOCCTab = ansEA$ResCC,
    GOMFTab = ansEA$ResMF,
    PathTab = ansEA$ResPat,
    nRGTab = Genelist, 
    nBar = 10
)
## png 
##   2
```

然后就可以得到一张条形图：
![Snipaste_2022-07-28_17-08-33](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-28_17-08-33.png)


## 生存分析

可以使用之前保存好的数据，也可以直接下载使用，临床数据不大，一般都能下载下来~


```r
# 临床数据可以像这样下载，直接就是一个数据框
clin.coad <- GDCquery_clinic("TCGA-COAD", "clinical")

dim(clin.coad)
## [1] 461  70
```

会直接得到这样一个数据框，非常方便：
![Snipaste_2022-07-28_18-09-10](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-28_18-09-10.png)

然后就可以进行生存分析了，比如我们用`gender`作为分组变量：


```r
TCGAanalyze_survival(
    data = clin.coad,
    clusterCol = "gender",
    main = "TCGA Set\n COAD",
    height = 10,
    width=10
)
## File saved as: survival.pdf
```

结果会得到这样一个图：
![Snipaste_2022-07-28_18-05-56](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-07-28_18-05-56.png)

肯定是不如自己用别的包画的好看~

也可以使用基因表达量进行分组，然后进行差异分析，只要把基因表达量数据加进去即可~

但是如果这样做的话，我们肯定是选择使用`survival`包做，比`TCGAanalyze_survival()`更加灵活好用~



除此之外，还可以进行火山图、热图、PCA图的绘制以及甲基化的一些简单分析，但是相比于它下载和准备的数据的功能，其他功能太弱了，都是对于其他包的封装，对于`TCGAbiolinks`的这些分析感兴趣的可以自行去官网学习哦~







