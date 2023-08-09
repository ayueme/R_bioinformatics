之前的2行代码提取表达矩阵由于大家的R语言水平参差不齐，导致很多新手会报错，于是我把前面的代码打包为一个脚本，1行代码就可以了！

脚本已上传到QQ群，需要的小伙伴加群下载即可~

**只需要1行代码就可以获取分别获取mRNA和lncRNA的counts/fpkm/tpm总计6种类型类型的表达矩阵以及临床信息，表达矩阵是标准形式，行是基因，列是样本，行名是gene symbol。**

>**使用这种方法有4个前提条件：**
>
>- `TCGAbiolinks`包的版本必须要在2.25.1以上
>- 需要使用`TCGAbiolinks`下载的数据或者按照这个教程下载的数据：可能是最适合初学者的TCGA下载教程
>- 必须按照这篇教程构建正确的路径：手动下载的TCGA数据也可以用TCGAbiolinks包整理
>- 脚本必须和`GDCdata`放在一个路径下

## 使用方法

加载需要的R包：

```{r,echo=FALSE}
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
```

加载脚本"getTCGAexpr.r"，这个脚本必须和**GDCdata**位于同一个位置。

![脚本位置必须对！](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929132733384.png)

加载这个脚本：

```{r}
source("getTCGAexpr.r")
```

使用函数，需要提供TCGA的癌症简称，比如：TCGA-LUSC。

```{r}
getTCGAexpr(project = "TCGA-LUSC")

##--------------------------------------
##o GDCquery: Searching in GDC database
##--------------------------------------
##Genome of reference: hg38
##--------------------------------------------
##oo Accessing GDC. This might take a while...
##--------------------------------------------
##ooo Project: TCGA-LUSC
##--------------------
##oo Filtering results
##--------------------
##ooo By data.type
##ooo By workflow.type
##----------------
##oo Checking data
##----------------
##ooo Checking if there are duplicated cases
##ooo Checking if there are results for the query
##-------------------
##o Preparing output
##-------------------
##|=====================================================|100%                      ##Completed after 16 s 
##Starting to add information to samples
## => Add clinical information to samples
## => Adding TCGA molecular information from marker papers
## => Information will have prefix 'paper_' 
##lusc subtype information from:doi:10.1038/nature11404
##Available assays in SummarizedExperiment : 
##  => unstranded
##  => stranded_first
##  => stranded_second
##  => tpm_unstrand
##  => fpkm_unstrand
##  => fpkm_uq_unstrand
##=> Saving file: output_expr/TCGA-LUSC_expr.rdata
##=> File saved
```

全程不到一分钟即可！

完成后会在当前目录多出一个**output_expr**文件夹，里面就是6个表达矩阵和临床信息：

![完成后会多出一个文件夹](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929132755500.png)

**output_expr**文件夹里面就是提取好的信息：

![提取好的表达矩阵和临床信息](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929133048502.png)

- `TCGA-LUSC_expr.rdata`：原始的se对象，所有信息都是从这里面提取的；
- `TCGA-LUSC_clinical.rdata`：TCGA-LUSC的临床信息；
- `TCGA-LUSC_lncRNA_expr_counts.rdata`：lncRNA的counts矩阵；
- `TCGA-LUSC_lncRNA_expr_fpkm.rdata`：lncRNA的fpkm矩阵；
- `TCGA-LUSC_lncRNA_expr_tpm.rdata`：lncRNA的tpm矩阵；
- `TCGA-LUSC_mRNA_expr_counts.rdata`：mRNA的counts矩阵；
- `TCGA-LUSC_mRNA_expr_fpkm.rdata`：mRNA的fpkm矩阵；
- `TCGA-LUSC_mRNA_expr_tpm.rdata`：mRNA的tpm矩阵；



表达矩阵示例：

![lncRNA的counts矩阵](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929133408560.png)

![mRNA的counts矩阵](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929133522374.png)



![mRNA的tpm矩阵](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929134624108.png)

![临床信息](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220929133658274.png)