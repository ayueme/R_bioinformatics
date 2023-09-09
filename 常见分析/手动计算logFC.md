前段时间有小伙伴问怎么手动计算`logFC`，今天说一下。

`logFC`是`log fold change`的缩写，也就是log之后的差异倍数。这个差异倍数意思是某个基因在A组表达量的平均值是B组表达量平均值的几倍。

这个东西的计算其实很简单的，就是常规的对数计算而已。

一般来说，我们用`tpm`或者`fpkm`时，通常都会先进行log2处理，在log2处理后的表达矩阵里，如果某个基因在A组表达量是x，在B组表达量是y，那么这个基因的`logFC = x - y`。

>这不是巧合，只是一个很简单的数学公式log(x/y)=log(x)-log(y)

## 准备数据

用`gse87466`这个GEO的数据做演示，下载整理的过程这次就不演示了。数据在粉丝qq群文件，需要的加群下载即可。


```r
library(easyTCGA)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(tidyr)

load(file = "G:/easyTCGA_test/gse87466.Rdata")
```

这是一个炎症性肠病的数据集，一共108个样本，21个`normal`，87个`uc`（ulcerative colitis）。

探针注释我已经提前做好了，但是有一些探针对应多个`symbol`，为了方便我这里直接删掉了：


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

exprSet <- exprSet[!grepl("/",rownames(exprSet)),]

group <- factor(group_list,levels = c("normal","UC"))
table(group)
## group
## normal     UC 
##     21     87
```

## 使用limma做差异分析

首先对这个数据做下差异分析，也是用`easyTCGA`包，1行代码即可，基因芯片数据也是支持的。

- 如果是`count`矩阵，会自动使用`DESeq2`、`limma`、`edgeR`进行差异分析；
- 如果是`tpm`、`fpkm`、`基因表达芯片数据`，它会自动检测需不需要进行log2转换，然后进行`wilcoxon`和`limma`的差异分析：


```r
library(easyTCGA)

diff_res <- diff_analysis(exprset = exprSet
                          , group = group
                          , is_count = F
                          )
## => log2 transform not needed
## => Running limma
## => Running wilcoxon test
## => Analysis done.

# limma的结果
diff_limma <- diff_res$deg_limma

diff_limma <- diff_limma %>% 
  arrange(desc(logFC))

head(diff_limma)
##            logFC   AveExpr         t      P.Value    adj.P.Val        B
## MMP3    5.125542  9.319122 12.754123 2.387453e-23 1.132335e-20 42.56916
## SLC6A14 5.025462  9.414419 21.582496 3.792392e-41 7.554445e-37 82.65877
## DUOX2   4.912367  9.917941 17.377496 3.456356e-33 1.721265e-29 64.81419
## DEFB4A  4.609705  8.710986  8.923018 1.217313e-14 5.214811e-13 22.81705
## CHI3L1  4.554672  9.777758 12.195694 4.281634e-22 1.332658e-19 39.72614
## S100A8  4.259343 10.016995 12.112299 6.602449e-22 1.906098e-19 39.29950
##         genesymbol
## MMP3          MMP3
## SLC6A14    SLC6A14
## DUOX2        DUOX2
## DEFB4A      DEFB4A
## CHI3L1      CHI3L1
## S100A8      S100A8
```

现在有很多文章中直接使用的`wilcoxon`检验，但是它并不能计算`logFC`：


```r
diff_wilc <- diff_res$deg_wilcoxon
head(diff_wilc)
##              pvalue genesymbol
##        9.202462e-03           
## IGL@   1.182949e-06       IGL@
## RPL41  1.708300e-02      RPL41
## RPL24  2.428482e-01      RPL24
## RPL37A 1.903898e-01     RPL37A
## RPS11  8.255636e-03      RPS11
```

## 自己计算logFC

根据前面的理论，我们可以自己计算`logFC`，思路就是分别计算某个基因在两组中的平均表达量，然后直接相减即可。

下面我们用`dplyr`中的`rowwise`操作实现这一过程，当然还有其他方法，选择自己喜欢的即可。


```r
library(dplyr)
library(tidyr)

uc <- colnames(exprSet)[group_list == "UC"]
normal <- colnames(exprSet)[group_list == "normal"]

logfc_df <- exprSet %>% 
  rowwise() %>% 
  mutate(mean_uc=rowMeans(across(all_of(uc))),
         mean_normal=rowMeans(across(all_of(normal))),
         logfc = mean_uc - mean_normal, # 这个就是logFC了
         .keep = "none"
            ) %>% 
  bind_cols(genesymbol = rownames(exprSet)) %>% 
  arrange(desc(logfc))
  
head(logfc_df)
## # A tibble: 6 × 4
## # Rowwise: 
##   mean_uc mean_normal logfc genesymbol
##     <dbl>       <dbl> <dbl> <chr>     
## 1   10.3         5.20  5.11 MMP3      
## 2   10.4         5.37  5.01 SLC6A14   
## 3   10.9         5.95  4.92 DUOX2     
## 4    9.60        5.00  4.60 DEFB4A    
## 5   10.7         6.10  4.56 CHI3L1    
## 6   10.8         6.55  4.29 S100A8
```

可以看到，我们手动计算的这个logfc和上面limma包计算的logFC基本上是一样的（有误差，可以忽略）哦。

## 参考资料

当然是来自于万能的生信菜鸟团啦，2015年的文章了：http://www.bio-info-trainee.com/1209.html

点击**阅读原文**可直接访问参考文章。

