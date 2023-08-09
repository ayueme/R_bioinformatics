众所周知，TCGA数据库改版了！！改的比之前更好用了！

对于常规转录组数据，主要是以下几点改变：
- 下载一次即可获得counts、TPM、FPKM三种类型的表达矩阵，再也不用单独下载了
- 自带gene symbol，不用自己找各种方法转换了
- 自带基因类型，可以直接区分mRNA和lncRNA了

`TCGAbiolinks`不仅是数据下载，它能访问、下载全部的TCGA数据（除了受限制的），用它下载的数据是*最新最全*的！这和直接去GDC官网，使用网页下载的方式是一样的。

除了常规的转录组数据，还包括甲基化数据、SNP数据、突变数据、临床数据等多种数据类型，还能进行数据分析，比如差异分析、生存分析、聚类等，除此之外，它也具有强大的绘图功能，可以直接绘制突变瀑布图等多种图形，是一个全面的TCGA包。

作为官方唯一推荐的专用下载及分析可视化一体的R包：`TCGAbiolinks`，也进行了相应的更新。

而`xena`的数据并不会及时更新，最新的数据还停留在2019年。

因为网络问题一直没怎么学习过这个强大的R包，最近数据更新了，学习下。

## 安装

需要安装版本在2.25.1以上的版本！


```r
# 经典2选1
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
```

>注意：目前`bioconductor`上面的`TCGAbiolinks`还停留在2.24.3版本，你需要安装开发版本哦~

如果你安装不成功，可以下载到本地安装，如果你不会本地安装，请翻看视频：xxxxxxxxxxx

## 使用

对网络有要求！

如果这一步都不成功，建议下面的就别运行了，因为很可能也不会成功。


```r
# 查看TCGA中33种癌症的简称
library(TCGAbiolinks)

projects <- TCGAbiolinks::getGDCprojects()$project_id ##获取癌症名字
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]

projects
##  [1] "TCGA-READ" "TCGA-UCS"  "TCGA-COAD" "TCGA-CESC" "TCGA-PAAD" "TCGA-ESCA"
##  [7] "TCGA-KIRP" "TCGA-PCPG" "TCGA-HNSC" "TCGA-BLCA" "TCGA-STAD" "TCGA-SARC"
## [13] "TCGA-CHOL" "TCGA-LAML" "TCGA-THYM" "TCGA-ACC"  "TCGA-SKCM" "TCGA-LUAD"
## [19] "TCGA-LIHC" "TCGA-KIRC" "TCGA-KICH" "TCGA-DLBC" "TCGA-PRAD" "TCGA-OV"  
## [25] "TCGA-MESO" "TCGA-LUSC" "TCGA-GBM"  "TCGA-UVM"  "TCGA-LGG"  "TCGA-BRCA"
## [31] "TCGA-TGCT" "TCGA-THCA" "TCGA-UCEC"
```

## 批量下载mRNA和lncRNA的数据

需要良好的网络环境，网络不好就别试了。全部数据40+G。


```r
sapply(projects, function(project){
  
  # 查询
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts"
                    )
  # 下载
  GDCdownload(query, method = "api", files.per.chunk = 100) #每次下载100个文件
  
  # 整理
  GDCprepare(query,save = T,save.filename = paste0(project,"_mRNA.Rdata"))
  }
  )
```

如果`query`能成功，但是下载成功，可以通过网页下载后，放在指定目录中，然后再运行`GDCprepare`函数也是可以成功的！

## 批量下载临床数据

也可以使用`GDCquery_clinic()`直接下载。


```r
sapply(projects, function(project){
  
  query <- GDCquery(project = project,
                    data.category = "Clinical", 
                    file.type = "xml")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  saveRDS(clinical,file = paste0(project,"_clinical.rds"))
})
```

使用方法做个小记录，可以通过不同的参数快速获取不同的临床数据：


```r
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

To get the following information please change the clinical.info argument
=> new_tumor_events: new_tumor_event 
=> drugs: drug 
=> follow_ups: follow_up 
=> radiations: radiation
```

## 批量下载miRNA


```r
sapply(projects, function(project){
  
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "miRNA Expression Quantification"
                    )
  
  GDCdownload(query)
  
  GDCprepare(query, save = T,save.filename = paste0(project,"_miRNA.Rdata"))
  
})
```


## 批量下载SNP


```r
sapply(projects, function(project){
  
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)
  
  GDCdownload(query)
  
  GDCprepare(query, save = T,save.filename = paste0(project,"_SNP.Rdata"))
  
})
```

## 批量下载CNV


```r
sapply(projects, function(project){
  
  query <- GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",              
    access = "open"
)
  
  GDCdownload(query)
  
  GDCprepare(query, save = T,save.filename = paste0(project,"_CNV.Rdata"))
  
})
```


## 批量下载甲基化

数据太大了，只下载一个COAD的演示一下~

β值矩阵：


```r
coad_methy <- GDCquery(
    project = "TCGA-COAD", 
    data.category = "DNA Methylation", 
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 27" # Illumina Human Methylation 450
    )
GDCdownload(coad_methy)
GDCprepare(coad_methy,save = T,save.filename="COAD_METHY_beta.Rdata")
```

IDAT：


```r
coad_methy <- GDCquery(
    project = "TCGA-COAD", 
    data.category = "DNA Methylation", 
    data.type = "Masked Intensities",
    platform = "Illumina Human Methylation 27", # Illumina Human Methylation 450
    legacy = FALSE
    )
GDCdownload(coad_methy)
GDCprepare(coad_methy,save = T,save.filename="COAD_METHY_idat.Rdata")
```


## 批量下载蛋白质数据


```r
sapply(projects, function(project){
  
  query <- GDCquery(
    project = project,
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification"
    )
  
  GDCdownload(query)
  
  GDCprepare(query, save = T,save.filename = paste0(project,"_protein.Rdata"))
  
})
```

亲测可用，我下载了2天1夜......

除此之外，还有其他数据可用，大家可以去[官网](https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html "TCGAbiolinks官方教程")学习~

