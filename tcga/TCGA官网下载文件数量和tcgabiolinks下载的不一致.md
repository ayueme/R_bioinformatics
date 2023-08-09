最近有小伙伴问我TCGA的表达矩阵整理问题，用到了我的一篇推文中的教程：

TCGA官网下载的数据也可以用TCGAbiolinks包搞定，只需2行代码！

但是总是遇到以下报错：

```R
# 查询这一步是需要的！即使网在烂，这一步应该可以成功的...
query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts"
                  )

# 下载这一步就不用了，我们是通过官网手动下载的~
#GDCdownload(query, files.per.chunk = 100) #每次下载100个文件

# 整理，网友在这一步遇到了一下报错！！！
GDCprepare(query,save = T,save.filename = "tcga_read.rdata")

Error in GDCprepare(query, save = T, save.filename = "tcga_read.rdata") : 
  I couldn't find all the files from the query. Please check if the directory parameter is right or `GDCdownload` downloaded the samples.
```

看这个报错提示是文件不全，让检查文件路径，在确定了文件路径和代码、网络都没有问题后，我觉得非常神奇！

理论上是不应该的呀！这个包就是用的官方的API下载的，不应该和官网直接下载的数据量不一样啊！

于是我赶紧检查了一下。

首先是看这个`query`一共查到了几个文件：

```R
tmp <- query$results[[1]]

# 查看查询到的文件夹名字这一列
head(tmp$id)

## [1] "00f55a16-0ee5-4939-8efb-de34e68d4ccd" "229e0c80-ada5-4fd3-8e93-a9bc1fac11a4"
## [3] "4b9b8b25-96e3-4667-8315-124711dcc1e0" "0ef12067-7e43-4d08-9374-a961430dd5ab"
## [5] "b9c1d14a-a169-4e73-a9c6-e884b005a160" "984e14e8-3272-4101-a1cb-81056eec7f8c"

# 一共177个
length(tmp$id)

## 177
```
也就是说我们通过代码的方式查询到了`TCGA-READ`一共177个文件！但是！网友下载的是91个！

我也赶紧去官网点点点看了一下，竟然也是91个！

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912205545012.png)

太神奇了，难道是`TCGAbiolinks`包出问题了吗？？？

冷静思考之后，我把网页中**Primary Site**中的打勾去掉了，然后就一切归于平静：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912210137181.png)

只要不选择**Primary Site**中的选项，就和`TCGAbiolinks`包下载的数据完全一样！又试了其他几个癌种，都是一样的了！

果然我还是太年轻，没见过世面啊！

**解决了这个小小的问题后，大家又可以愉快的只用2行代码解决表达矩阵的整理问题了！**

接下来还是使用官网网页下载，然后自己新建指定文件路径，就可以用2行代码搞定表达矩阵了：

```R
# 查询
query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts"
                  )

# 整理
GDCprepare(query,save = T,save.filename = "tcga_read.rdata")

|==========================================|100%    Completed after 14 s 
Starting to add information to samples
 => Add clinical information to samples
 => Adding TCGA molecular information from marker papers
 => Information will have prefix 'paper_' 
read subtype information from:doi:10.1038/nature11252
Available assays in SummarizedExperiment : 
  => unstranded
  => stranded_first
  => stranded_second
  => tpm_unstrand
  => fpkm_unstrand
  => fpkm_uq_unstrand
=> Saving file: tcga_read.rdata
=> File saved
class: RangedSummarizedExperiment 
dim: 60660 177 
metadata(1): data_release
assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1
  ENSG00000288675.1
rowData names(10): source type ... hgnc_id havana_gene
colnames(177): TCGA-AF-3911-01A-01R-1736-07 TCGA-DY-A1DC-01A-31R-A155-07 ...
  TCGA-AF-2692-11A-01R-A32Z-07 TCGA-EI-6882-01A-11R-1928-07
colData names(107): barcode patient ... paper_vascular_invasion_present
  paper_vital_status
```

全程不到1分钟即可完成，里面包含了fpkm/tpm/counts的表达矩阵、以及超级详细的临床信息！可以参考另一篇推文：超简单的表达矩阵提取。

舒服！