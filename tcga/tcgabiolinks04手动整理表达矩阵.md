很多人因为网络原因不能使用`TCGAbiolinks`这个神包下载TCGA的RNA-seq数据，只能通过浏览器访问GDC TCGA的官网进行下载，而下载后得到的是一个个文件夹，对于如何整理成一个表达矩阵也是很麻烦的。

今天给大家介绍一个简单点的方法，使用`TCGAbiolinks`包整理你通过浏览器官网下载的rna-seq数据。

通常大家通过浏览器下载后会得到下面的这种很多个文件夹：
![Snipaste_2022-08-03_17-51-57](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-08-03_17-51-57.png)

每个文件夹里是一个样本的表达量数据，tsv格式的：
![Snipaste_2022-08-03_17-52-14](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-08-03_17-52-14.png)

这时候你可以通过之前介绍过的方法得到表达矩阵：新版TCGA数据库表达矩阵整理

但是这个方法对于新手还是不够友好，尤其是根据Json文件匹配数据时，但是TCGA表达量数据又是很常用的，这个操作还是很高频的需求。

前几天学习`TCGAbiolinks`包时意外发现，**即使是手动下载的数据，只要构建合适的路径，也是可以通过`GDCprepare()`函数进行整理从而简单的得到表达矩阵的！**

`TCGAbiolinks`包下载的表达量数据的文件路径是有规律的，如果你没有特别指明，通常是位于**GDCdata\TCGA-COAD\harmonized\Transcriptome_Profiling\Gene_Expression_Quantification**这个路径下的。

这个包下载数据就是三板斧操作，`query`,`download`,`prepare`，而且最后`GDCprepare()`需要的还是`GDCquery()`得到的对象，因此我们完全可以通过构建一个适合它的路径，让`GDC_prepare()`帮我们整理成表达矩阵！

比如我上面的各个样本文件夹的路径在我的电脑中是这样的：**G:\tcga\GDCdata\TCGA-COAD\harmonized\Transcriptome_Profiling\Gene_Expression_Quantification**，我的`get_expr.R`脚本是放在**G:\tcga**这个路径下的。

脚本内容如下：

```{r}
library(TCGAbiolinks)

## =============================================================
## ______  ___  ____   ___                                        
##   ||   |    |      |   | |    o  __  |   o  _         __         
##   ||   |    | ___  |___| |__  | |  | |   | | | | |_/ |__         
##   ||   |___ |____| |   | |__| | |__| |__ | | |_| | \  __|       
## ------------------------------------------------------------
## Query, download & analyze - GDC                  
## Version:2.25.2
## ==============================================================


# 查询这一步是需要的！即使网在栏，这一步应该可以成功的...
query <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts"
                  )
# 下载这一步就不用了，我们是通过官网手动下载的~
# GDCdownload(query, files.per.chunk = 100) #每次下载100个文件
  
# 整理
GDCprepare(query,save = T,save.filename = "example.rdata")

##|===============================================================================|100%   ##                   Completed after 1 m 
##Starting to add information to samples
## => Add clinical information to samples
## => Adding TCGA molecular information from marker papers
## => Information will have prefix 'paper_' 
##coad subtype information from:doi:10.1038/nature11252
##Available assays in SummarizedExperiment : 
##  => unstranded
##  => stranded_first
##  => stranded_second
##  => tpm_unstrand
##  => fpkm_unstrand
##  => fpkm_uq_unstrand
##=> Saving file: example.rdata
##=> File saved
```

这样我们的数据就整理好了：
![Snipaste_2022-08-03_18-22-19](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-08-03_18-22-19.png)

下次使用直接load即可：

```{r}
rm(list = ls())
load(file = "example.rdata")

se <- data
se

class: RangedSummarizedExperiment 
dim: 60660 521 
metadata(1): data_release
assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1 ENSG00000288675.1
rowData names(10): source type ... hgnc_id havana_gene
colnames(521): TCGA-A6-5664-01A-21R-1839-07 TCGA-D5-6530-01A-11R-1723-07 ...
  TCGA-A6-2683-01A-01R-0821-07 TCGA-A6-2683-11A-01R-A32Z-07
colData names(107): barcode patient ... paper_vascular_invasion_present paper_vital_status
```

**这个`se`就是我们之前介绍过的SummarizedExperiment对象，你可以对它进行各种操作，得到counts矩阵、tpm矩阵、fpkm矩阵都是小事一桩，犹如探囊取物一般简单流畅！** 详情可参考之前的推文。

>关于TCGA表达矩阵提取，告诉我，你还有哪里搞不定！？