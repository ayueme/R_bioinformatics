前几天看到了生信技能树的推文：[什么情况下我们可以修改基因名字](https://mp.weixin.qq.com/s/vPE1Mvj-GHgktwAU3NWwmg)

里面提到了2个函数很好用：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230516152439283.png)

其实这个需求我知道在小洁老师的R包`tinyarray`里有函数可以实现：[宝藏R包tinyarray：常用图表一键收走](https://mp.weixin.qq.com/s/cdcuwqajMqp1t9Gs0-jEfA)

我还知道果子老师在2018年就介绍过这个技能了：[多个基因在多亚组疾病中的展示](https://mp.weixin.qq.com/s/DOvSutbqotCK5JDHFdjFaw)

感觉这个函数非常实用，于是我也想写一个。

写好之后我把它们加在`easyTCGA`包里了。

目前`easyTCGA`已经完美支持TCGA常见5种组学数据的下载和整理：**gene**、**miRNA**、**Copy Number Variation**、

**simple nucleotide variation**、**DNA methylation**

**批量生存分析也默认支持最佳截点了**。

大家可以去github了解详情：https://github.com/ayueme/easyTCGA

今天主要说下几个可视化小函数。

## 准备数据

以`TCGA-BRCA`为例。


```r
library(easyTCGA)
#getmrnaexpr("TCGA-BRCA") # 下载只要1行代码
load(file = "G:/easyTCGA_test/output_mRNA_lncRNA_expr/TCGA-BRCA_mrna_expr_tpm.rdata")
load(file = "G:/easyTCGA_test/output_mRNA_lncRNA_expr/TCGA-BRCA_clinical.rdata")
```

准备表达数据，分子，和组别。


```r
expr <- log2(mrna_expr_tpm+0.1)
marker <- "CXCL1"
markers <- c("CXCL1","TP53","BRAF","EGFR","CTLA4","VEGFB","NTRK2")
sample_group <- ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumor","normal")
table(sample_group)
## sample_group
## normal  tumor 
##    113   1118
sample_groups <- clin_info$paper_BRCA_Subtype_PAM50 #临床信息非常丰富！
table(sample_groups)
## sample_groups
##  Basal   Her2   LumA   LumB Normal 
##    197     82    571    209     40
```

## plot_gene

首先是`plot_gene`。可以实现：*任意数量基因在任意癌种（TCGA33种其中之一都可以）的任意分组中的表达量箱线图*

1个基因在两个组的表达量：


```r
res <- plot_gene(expr = expr,marker,sample_group)
```

![plot of chunk unnamed-chunk-3](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-3-174794620.png)

多个基因在两个组的表达量：


```r
res <- plot_gene(expr = expr,markers,sample_group)
```

![plot of chunk unnamed-chunk-4](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-4-174794620.png)

1个基因在多个组的表达量：


```r
res <- plot_gene(expr = expr,marker,sample_groups)
```

![plot of chunk unnamed-chunk-5](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-5-174794620.png)

多个基因在多个组的表达量：


```r
res <- plot_gene(expr = expr,markers,sample_groups)
```

![plot of chunk unnamed-chunk-6](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-6-174794620.png)

当然大家对美的追求是无止境的，所以每个画图函数都可以返回画图数据，方便你DIY。


```r
head(res)
##                      sample_id  group markers expression
## 1 TCGA-B6-A0RH-01A-21R-A115-07   Her2   CXCL1 -2.9072513
## 2 TCGA-BH-A1FU-11A-23R-A14D-07   <NA>   CXCL1  2.3879966
## 3 TCGA-BH-A1FU-01A-11R-A14D-07 Normal   CXCL1  3.7425560
## 4 TCGA-AR-A0TX-01A-11R-A084-07   Her2   CXCL1  1.5995082
## 5 TCGA-A1-A0SE-01A-11R-A084-07   LumA   CXCL1 -0.7147754
## 6 TCGA-BH-A1FC-11A-32R-A13Q-07   <NA>   CXCL1  5.1622906
```

## plot_gene_paired

任意基因在某一癌种配对样本中的表达量箱线图；


```r
pair_sam <- plot_gene_paired(expr = expr,marker = markers)
```

![plot of chunk unnamed-chunk-8](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-8-174794620.png)


```r
pair_sam <- plot_gene_paired(expr = expr,marker = "CXCL1")
```

![plot of chunk unnamed-chunk-9](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-9-174794620.png)

如果你不满意也可以自己用返回的数据画图。

## plot_km

根据任意基因的表达量分组，并画出K-M生存曲线（支持最佳截点）

需要先准备下临床数据，需要一个数据框，只含有两列，列名是`time`和`event`，event用1表示死亡，0表示存活。


```r
# 准备临床数据
clin <- clin_info[,c("days_to_last_follow_up","vital_status")]
names(clin) <- c("time","event")
clin$event <- ifelse(clin$event=="Dead",1,0)
```

画图，默认根据最佳截点，否则根据中位数，支持返回数据


```r
res <- plot_KM(exprset = expr, marker = marker, clin = clin,optimal_cut = T)
```

![plot of chunk unnamed-chunk-11](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-174794620.png)


```r
res <- plot_KM(exprset = expr, marker = marker, clin = clin,optimal_cut = F)
```

![plot of chunk unnamed-chunk-12](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-12-174794620.png)

最佳截点效果还是很明的！

批量生存分析也是默认支持最佳截点的哦。
