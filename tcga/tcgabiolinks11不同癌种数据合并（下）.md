前几天在我们的QQ交流群中有小伙伴提出合并LUAD和LUSC数据时报错：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230218121311618.png)

使用的是我推文中的教程：[新版TCGA数据库不同癌种的组学数据合并](https://mp.weixin.qq.com/s/0hcQ1m_9l1TtvXgEG20F5Q)

上面这篇教程是用COAD和READ演示的，没有任何问题，但是LUAD和LUSC就出问题了，吓得我赶紧去试了一下，果然是报错的。

所以这里给大家提供另一种更加简单的方法，只要下载时同时写上LUAD和LUSC即可，整个过程交给`TCGAbiolinks`包。代码如下：

```R
# 加载R包
library(TCGAbiolinks)

## 合并，同时写TCGA-LUAD和TCGA-LUSC即可
query <- GDCquery(project = c("TCGA-LUAD","TCGA-LUSC"),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts"
                  )
# 下载,网络不好的这一步经常报错，建议使用我们推文中的其他方法
GDCdownload(query, files.per.chunk = 100) #每次下载100个文件

# 整理
GDCprepare(query,save = T,
           save.filename = "G:/tcga/lusc_luad测试/TCGA_LUAD_LUSC_mRNA.Rdata")

```

其中，下载这一步，网络不好的经常报错，建议使用我们推文中的其他方法：

[可能是最适合初学者的TCGA官网下载和表达矩阵整理教程](https://mp.weixin.qq.com/s/rbnWvstRsfhbi9il-qSYpQ)

整个过程非常丝滑，没有任何问题：

```R
> GDCprepare(query,save = T,save.filename = "G:/tcga/lusc_luad测试/TCGA_LUAD_LUSC_mRNA.Rdata")

|=============|100%   Completed after 30 s 

Starting to add information to samples
 => Add clinical information to samples
 => Adding TCGA molecular information from marker papers
 => Information will have prefix 'paper_' 
luad subtype information from:doi:10.1038/nature13385
lusc subtype information from:doi:10.1038/nature11404
Available assays in SummarizedExperiment : 
  => unstranded
  => stranded_first
  => stranded_second
  => tpm_unstrand
  => fpkm_unstrand
  => fpkm_uq_unstrand
=> Saving file: G:/tcga/lusc_luad测试/TCGA_LUAD_LUSC_mRNA.Rdata
=> File saved
class: RangedSummarizedExperiment 
dim: 60660 1153 
metadata(1): data_release
assays(6): unstranded stranded_first ... fpkm_unstrand
  fpkm_uq_unstrand
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
  ENSG00000288674.1 ENSG00000288675.1
rowData names(10): source type ... hgnc_id havana_gene
colnames(1153): TCGA-78-7156-01A-11R-2039-07
  TCGA-44-6774-01A-21R-1858-07 ...
  TCGA-52-7809-01A-21R-2125-07 TCGA-85-7843-01A-11R-2125-07
colData names(95): barcode patient ...
  paper_Homozygous.Deletions paper_Expression.Subtype
```

保存后的数据是一个`summarisedExperiment`对象，这个对象里面包含了表达矩阵（包含counts、tpm、fpkm）和临床信息等多种信息，可以参考我们之前的推文了解它：

[新版TCGAbiolinks包学习：表达矩阵提取（mRNA/lncRNA/counts/tpm/fpkm）](https://mp.weixin.qq.com/s/wI0_GyVl5LiKAjX5C3f-NQ)

同理，这种方法也可以用于结直肠癌的合并。

easy！