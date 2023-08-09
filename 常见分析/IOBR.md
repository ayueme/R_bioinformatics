免疫浸润一直是生信数据挖掘的重点，免疫浸润的各种算法也是大家学习生信数据挖掘必学的知识点。

前几天看到还有单独介绍各种免疫浸润方法实现的文章，其实现在不用那么麻烦了，已经有人帮我们写好了整合的R包，它就是：`IOBR`。

这个包**整合了常见的10种免疫浸润方法，只需要1行代码即可实现，而且输出的格式统一，方便进行可视化和数据整合操作！**

完全不需要你自己分别进行操作，极大地简化了进行免疫浸润分析的步骤。

下面就给大家演示如何使用1行代码实现8种免疫浸润方法！

首先是安装R包。

## 安装

先安装依赖包，我发现还有很多人被困在R包安装这一步，实在是不应该，我专门做了视频版教程帮大家解决R包安装的问题，建议R包安装还有问题的赶紧去看：[可能是最好用的R包安装教程](https://www.bilibili.com/video/BV11g411o7be/?vd_source=2a81c5384000daae61949f58079f1cfd)

最近（2023.05.09）R更新到了4.3.0版本，bioconductor也已经更新到了3.17版本，版本不匹配就会安装失败哈，初学者要注意！


```r
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}
```

再安装`IOBR`包：


```r
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")
```

这个R包不仅可以进行各种免疫浸润分析，还有很多使用的分析，这些我们放到以后再说，今天就主要介绍它的免疫浸润分析功能。


```r
library(IOBR)
```

## 下载数据

我们就用`easyTCGA`包下载`TCGA-COAD`的数据进行演示，这个包专门为小白解决TCGA数据下载和整理而写，可以参考：[easyTCGA：让初学者也能享受“征服”TCGA的喜悦](https://mp.weixin.qq.com/s/kvGYVCOSBgKqVaeQU01JcA)


```r
library(easyTCGA)
getmrnaexpr("TCGA-COAD")
```

下载数据超级简单，如果你的网络通畅，只需要几分钟就可以得到6种表达矩阵和临床信息，且同时保存为`rdata`和`csv`两种格式。

我们使用`tpm`数据进行演示，首先进行`log2`转换，然后去掉一些低质量的基因。


```r
load(file = "output_mRNA_lncRNA_expr/TCGA-COAD_mrna_expr_tpm.rdata")
expr_coad <- log2(mrna_expr_tpm+0.1)
expr_coad <- expr_coad[apply(expr_coad,1,sd)>0.5,]

expr_coad[1:4,1:4]
##        TCGA-D5-6540-01A-11R-1723-07 TCGA-AA-3525-11A-01R-A32Z-07
## MT-CO2                     14.47089                     14.06391
## MT-CO3                     14.42028                     13.94628
## MT-ND4                     14.35429                     13.24382
## MT-CO1                     14.63867                     13.99546
##        TCGA-AA-3525-01A-02R-0826-07 TCGA-AA-3815-01A-01R-1022-07
## MT-CO2                     13.98376                     14.84026
## MT-CO3                     14.17448                     14.66966
## MT-ND4                     13.45308                     14.51663
## MT-CO1                     13.60652                     14.53353
dim(expr_coad)
## [1] 17170   524
```

这样我们的表达矩阵就整理好了，下面就开始演示1行代码实现8种免疫浸润分析。

在此之前，先看看支持哪8种方法：


```r
tme_deconvolution_methods
##         MCPcounter               EPIC              xCell          CIBERSORT 
##       "mcpcounter"             "epic"            "xcell"        "cibersort" 
## CIBERSORT Absolute                IPS           ESTIMATE                SVR 
##    "cibersort_abs"              "ips"         "estimate"              "svr" 
##               lsei              TIMER          quanTIseq 
##             "lsei"            "timer"        "quantiseq"
```

大家常见的方法全都囊括了，再也不用费功夫去找各种代码了！

并且每种方法都给出了参考文献和开源证书：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230508194953468.png)

大家可能发现怎么没有`ssGSEA`啊？这么优秀的算法，肯定有的，下面会给大家介绍的！

## 1行代码实现8种免疫浸润分析

所有方法只要通过`deconvo_tme()`函数即可完成，只要在`method`中选择不同的方法即可！也不需要自己再去寻找各种方法对应的免疫细胞集，都是内置好的！


```r
# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr_coad,
                            method = "mcpcounter"
                            )
## 
## >>> Running MCP-counter

# EPIC
im_epic <- deconvo_tme(eset = expr_coad,
                       method = "epic",
                       arrays = F
                       )
## 
## >>> Running EPIC

# xCell
im_xcell <- deconvo_tme(eset = expr_coad,
                        method = "xcell",
                        arrays = F
                        )
## 
## >>> Running xCell
## [1] "Num. of genes: 9854"
## Estimating ssGSEA scores for 489 gene sets.
## [1] "Calculating ranks..."
## [1] "Calculating absolute values from ranks..."
  |                                                                            
  |======================================================================| 100%

# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr_coad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
                            )
## 
## >>> Running CIBERSORT

# IPS
im_ips <- deconvo_tme(eset = expr_coad,
                      method = "ips",
                      plot = F
                      )
## 
## >>> Running Immunophenoscore

# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr_coad,
                            method = "quantiseq",
                            scale_mrna = T
                            )
## 
## Running quanTIseq deconvolution module
## Gene expression normalization and re-annotation (arrays: FALSE)
## Removing 17 noisy genes
## Removing 15 genes with high expression in tumors
## Signature genes found in data set: 134/138 (97.1%)
## Mixture deconvolution (method: lsei)
## Deconvolution sucessful!
# ESTIMATE
im_estimate <- deconvo_tme(eset = expr_coad,
                           method = "estimate"
                           )
## 
## >>> Running ESTIMATE
## [1] "Merged dataset includes 9232 genes (1180 mismatched)."
## [1] "1 gene set: StromalSignature  overlap= 136"
## [1] "2 gene set: ImmuneSignature  overlap= 139"

# TIMER
im_timer <- deconvo_tme(eset = expr_coad
                        ,method = "timer"
                        ,group_list = rep("coad",dim(expr_coad)[2])
                        )
## ## Enter batch mode
## ## Loading immune gene expression
## [1] "Outlier genes: ACTB ACTG1 ACTG2 CLCA1 COL1A1 COL3A1 CXCL8 DEFA5 DEFA6 FTL GAPDH H4C3 HLA-C HLA-DRA JCHAIN LCN2 LGALS1 MT-ATP6 MT-ATP8 MT-CO1 MT-CO2 MT-CO3 MT-CYB MT-ND1 MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-ND5 MT-ND6 OLFM4 PAEP PI3 PIGR PLA2G2A PPBP PRSS2 REG1A REG1B REG3A RPL8 RPS11 RPS12 RPS18 RPS2 RPS21 RPS6 S100A6 S100P SLC25A6 SPINK4 TFF1 TFF3 TMSB10 TMSB4X"
## ## Removing the batch effect of C:\Users\liyue\AppData\Local\Temp\Rtmpw79DY6\file108c19af34b0
## Found2batches
## Adjusting for0covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data

# 需要提供reference，暂不演示！
#im_svr <- deconvo_tme(eset = expr_coad
#                      ,method = "svr"
#                      ,arrays = F
#                      ,reference = 
#                      )
#im_lsei <- deconvo_tme(eset = expr_coad
#                       ,method = "lsei"
#                       ,arrays = F
#                       )
```

这样8种免疫浸润方法就做好了！是不是太easy了！

而且作者对结果进行了整理，都是统一的格式，方便你直接进行合并操作！随便给大家放几个看看：


```r
dim(im_cibersort)
## [1] 524  26
im_cibersort[1:4,1:4]
## # A tibble: 4 × 4
##   ID                           B_cells_naive_CIBERSORT B_cells_memory_…¹ Plasm…²
##   <chr>                                          <dbl>             <dbl>   <dbl>
## 1 TCGA-D5-6540-01A-11R-1723-07                  0.0504                 0   0    
## 2 TCGA-AA-3525-11A-01R-A32Z-07                  0.116                  0   0.183
## 3 TCGA-AA-3525-01A-02R-0826-07                  0.0655                 0   0.149
## 4 TCGA-AA-3815-01A-01R-1022-07                  0.0390                 0   0    
## # … with abbreviated variable names ¹​B_cells_memory_CIBERSORT,
## #   ²​Plasma_cells_CIBERSORT

dim(im_xcell)
## [1] 524  68
im_xcell[1:4,1:4]
## # A tibble: 4 × 4
##   ID                           aDC_xCell Adipocytes_xCell Astrocytes_xCell
##   <chr>                            <dbl>            <dbl>            <dbl>
## 1 TCGA-D5-6540-01A-11R-1723-07    0.119          1.32e-18         2.03e- 2
## 2 TCGA-AA-3525-11A-01R-A32Z-07    0.451          0                1.13e-17
## 3 TCGA-AA-3525-01A-02R-0826-07    0.0792         6.63e-19         2.08e-19
## 4 TCGA-AA-3815-01A-01R-1022-07    0.418          4.95e-19         0
```

看到了吗？所有的结果都是：**行是样本，列是细胞，行名是样本名，列名是细胞名，并且每种细胞名后面都有相应的方法名字的后缀**，方便标识！

而且提供了方便的可视化函数，比如可视化`cibersort`的结果：


```r
library(tidyr)
# 取前12个样本做演示
res<-cell_bar_plot(input = im_cibersort[1:12,], title = "CIBERSORT Cell Fraction")
## There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap
## >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
```

![plot of chunk unnamed-chunk-9](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-9-174361934.png)

关于免疫浸润的可视化，我们以后再专门介绍！

有了这样的数据方便大家直接进行合并操作！比如：


```r
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID")

tme_combine[1:4,1:4]
## # A tibble: 4 × 4
##   ID                           T_cells_MCPcounter CD8_T_cells_MCPcounter Cytot…¹
##   <chr>                                     <dbl>                  <dbl>   <dbl>
## 1 TCGA-D5-6540-01A-11R-1723-07              0.876                   1.47   0.323
## 2 TCGA-AA-3525-11A-01R-A32Z-07              2.37                    1.90   0.655
## 3 TCGA-AA-3525-01A-02R-0826-07              0.519                  -1.00  -0.269
## 4 TCGA-AA-3815-01A-01R-1022-07              2.78                    1.98   0.851
## # … with abbreviated variable name ¹​Cytotoxic_lymphocytes_MCPcounter
dim(tme_combine)
## [1] 524 138
```

看看这统一的语法和格式，真是令人赏心悦目啊！

那`ssGSEA`实现免疫浸润该怎么做呢？其实这个是该包的另一个重要的函数`calculate_sig_score`实现的。

因为`ssGSEA`是计算**富集分数**的算法，所以作者把它放到了这里，大家常见的这种方法对应的可能是28个免疫细胞的基因集，但是这个算法的牛逼之处就在于，你给它不同的基因集，它就可以计算不同的富集分数，并不是局限于28种细胞的那个基因集（上面8种方法都是局限于1种基因集的）。

除`ssGSEA`之外，`calculate_sig_score`还提供了`pca/zscore/intergration`一共4种方法计算分数，真的是非常强！


```r
## 基于ssGSEA计算免疫浸润分数
load(file = "../000files/ssGSEA28.Rdata")
im_ssgsea <- calculate_sig_score(eset = expr_coad
                                 , signature = cellMarker # 这个28种细胞的文件需要自己准备
                                 , method = "ssgsea" # 选这个就好了
                                 )
## 
## >>> Calculating signature score using ssGSEA method
## >>> log2 transformation is not necessary
## Estimating ssGSEA scores for 28 gene sets.
## [1] "Calculating ranks..."
## [1] "Calculating absolute values from ranks..."
## 
  |                                                                            
  |======================================================================| 100%
## 
## [1] "Normalizing..."
im_ssgsea[1:4,1:4]
## # A tibble: 4 × 4
##   ID                           `Activated B cell` `Activated CD4 T cell` Activ…¹
##   <chr>                                     <dbl>                  <dbl>   <dbl>
## 1 TCGA-3L-AA1B-01A-11R-A37K-07             -0.192                0.0120   0.170 
## 2 TCGA-4N-A93T-01A-11R-A37K-07             -0.249                0.00395  0.121 
## 3 TCGA-4T-AA8H-01A-11R-A41B-07             -0.358                0.0552   0.0838
## 4 TCGA-5M-AAT4-01A-11R-A41B-07             -0.399                0.0909   0.107 
## # … with abbreviated variable name ¹​`Activated CD8 T cell`
```

这个28种细胞的文件需要自己准备，如果你需要，直接在本号后台回复`ssgsea28`即可得到这个文件。

我已经在作者的`github`提了issue，不知道会不会通过...

这个`calculate_sig_score`函数得到的结果也是和上面一模一样的格式！也可以直接进行合并操作！


```r
tme_combine <- tme_combine %>% 
  inner_join(im_ssgsea, by = "ID")
```

真的是太好用了！

前面我们说过，`calculate_sig_score`可以根据不同的基因集计算不同的分数，作者为了方便大家，直接**内置了已经公开发表的255种基因集，包括肿瘤微环境相关的、代谢相关的、m6A、外泌体、铁死亡和错配修复等**。

可通过以下函数查看，大家可以自己探索一下：


```r
# 总的
signature_collection
# 代谢相关
signature_metabolism
# 微环境相关
signature_tme
# 肿瘤相关
signature_tme
```

每一个`signature`都是list结构：


```r
signature_metabolism[1:4]
## $Cardiolipin_Metabolism
## [1] "CKMT1A"         "CKMT1B"         "CYCS"           "NME4"          
## [5] "TAZ"            "TMEM256-PLSCR3"
## 
## $Cardiolipin_Biosynthesis
## [1] "CRLS1"  "PGS1"   "PTPMT1" "TAMM41"
## 
## $Cholesterol_Biosynthesis
##  [1] "ACAT2"   "CYP51A1" "DHCR24"  "DHCR7"   "EBP"     "FDFT1"   "FDPS"   
##  [8] "GGPS1"   "HMGCR"   "HMGCS1"  "HSD17B7" "IDI1"    "LBR"     "LIPA"   
## [15] "LSS"     "MSMO1"   "MVD"     "MVK"     "NSDHL"   "PMVK"    "SC5D"   
## [22] "SOAT1"   "SQLE"    "TM7SF2" 
## 
## $Citric_Acid_Cycle
##  [1] "ACLY"   "ACO1"   "ACO2"   "CS"     "DLAT"   "DLD"    "DLST"   "FH"    
##  [9] "IDH1"   "IDH2"   "IDH3A"  "IDH3B"  "IDH3G"  "MDH1"   "MDH2"   "MPC1"  
## [17] "OGDH"   "OGDHL"  "PC"     "PCK1"   "PCK2"   "PDHA1"  "PDHA2"  "PDHB"  
## [25] "SDHA"   "SDHB"   "SDHC"   "SDHD"   "SUCLA2" "SUCLG1" "SUCLG2" "PDHX"
signature_metabolism[[1]]
## [1] "CKMT1A"         "CKMT1B"         "CYCS"           "NME4"          
## [5] "TAZ"            "TMEM256-PLSCR3"
```

每种`signature`还给出了参考文献，方便大家引用，绝对好用！


```r
signature_collection_citation[1:2,]
## # A tibble: 2 × 6
##   Signatures      `Published year` Journal Title                     PMID  DOI  
##   <chr>                      <dbl> <chr>   <chr>                     <chr> <chr>
## 1 CD_8_T_effector             2018 Nature  TGFβ attenuates tumour r… 2944… 10.1…
## 2 DDR                         2018 Nature  TGFβ attenuates tumour r… 2944… 10.1…
```

这么优秀的R包，大家快去用起来，使用时别忘记引用：

`Zeng D, Ye Z, Shen R, Yu G, Wu J, Xiong Y,…, Liao W (2021) IOBR: Multi-Omics Immuno-Oncology Biological Research to Decode Tumor Microenvironment and Signatures. Frontiers in Immunology. 12:687975. doi: 10.3389/fimmu.2021.687975`

## 免疫浸润结果可视化

在之前的推文中我们介绍了2行代码实现9种免疫浸润方法，今天给大家介绍下常见的免疫浸润结果的可视化。

就以大家最常见的`cibersort`为例进行介绍。

首先大家要对每种免疫浸润方法的结果有一个大体的认知，比如`cibersort`的结果是各种免疫细胞在样本中的比例，所以一个样本中所有的免疫细胞比例加起来总和是1！

但是`ssGSEA`就不是这样了。

只有理解了结果是什么样的，你才能选择合适的可视化方法。**数就是图，图就是数**


```r
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)

# 首先变为长数据
cibersort_long <- im_cibersort %>% 
  select(`P-value_CIBERSORT`,Correlation_CIBERSORT, RMSE_CIBERSORT,ID,everything()) %>% 
  pivot_longer(- c(1:4),names_to = "cell_type",values_to = "fraction") %>% 
  dplyr::mutate(cell_type = gsub("_CIBERSORT","",cell_type),
                cell_type = gsub("_"," ",cell_type))

head(cibersort_long[,4:6])
## # A tibble: 6 × 3
##   ID                           cell_type                  fraction
##   <chr>                        <chr>                         <dbl>
## 1 TCGA-D5-6540-01A-11R-1723-07 B cells naive                0.0504
## 2 TCGA-D5-6540-01A-11R-1723-07 B cells memory               0     
## 3 TCGA-D5-6540-01A-11R-1723-07 Plasma cells                 0     
## 4 TCGA-D5-6540-01A-11R-1723-07 T cells CD8                  0.119 
## 5 TCGA-D5-6540-01A-11R-1723-07 T cells CD4 naive            0     
## 6 TCGA-D5-6540-01A-11R-1723-07 T cells CD4 memory resting   0.0951
```

如果你是初学者，也可以直接用`IOBR`为大家提供的`cell_bar_plot()`函数，可以直接帮你转换为长数据，并且还可以把图画出来。


```r
res_cibersort <- cell_bar_plot(im_cibersort)
## There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap
## >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
```

![plot of chunk unnamed-chunk-17](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-17-174361934.png)

但是这个图你可能不能接受，没关系，转换好的长数据已经在`res_cibersort`这个对象里了。


```r
head(res_cibersort$data)
##                             ID     cell_type   fraction
## 1 TCGA-D5-6540-01A-11R-1723-07 B cells naive 0.05044811
## 2 TCGA-AA-3525-11A-01R-A32Z-07 B cells naive 0.11573408
## 3 TCGA-AA-3525-01A-02R-0826-07 B cells naive 0.06545120
## 4 TCGA-AA-3815-01A-01R-1022-07 B cells naive 0.03903166
## 5 TCGA-D5-6923-01A-11R-A32Z-07 B cells naive 0.02560987
## 6 TCGA-G4-6322-01A-11R-1723-07 B cells naive 0.08727238
```

和我们自己转换的是差不多的，这个数据就可以直接使用`ggplot2`自己画图了：


```r
p1 <- cibersort_long %>% 
  ggplot(aes(ID,fraction))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
        )
p1
```

![plot of chunk unnamed-chunk-19](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-19-174361934.png)

除了这个最常见的堆叠条形图，还可以画箱线图，热图。


```r
# 有顺序的箱线图
library(forcats)

p2 <- ggplot(cibersort_long,aes(fct_reorder(cell_type, fraction),fraction,fill = cell_type)) + 
  geom_boxplot() + 
  #geom_jitter(width = 0.2,aes(color=cell_type))+
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = palette4)
p2
```

![plot of chunk unnamed-chunk-20](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-20-174361934.png)

如果你的样本有不同分组，就可以根据画出分组箱线图。比如我这里就根据`tumor/normal`把样本分组，然后再组间进行非参数检验，并添加P值。

这些都是R语言基础操作，本号的**可视化**合集中介绍了太多这些基本绘图知识了。


```r
library(ggpubr)
library(stringr)

# 分组
cibersort_long$Group = ifelse(as.numeric(str_sub(cibersort_long$ID,14,15))<10,"tumor","normal")

p3 <- ggplot(cibersort_long,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 0.4)
p3
```

![plot of chunk unnamed-chunk-21](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-21-174361934.png)

热图也是一样的easy，而且有了`tidyHeatmap`的加持，直接使用长数据即可，不用在变为宽数据了！

可以参考文章：[tidyHeatmap完美使用长数据的热图可视化](https://mp.weixin.qq.com/s/-PkppnmepYmlm-RJlxykeA)


```r
library(tidyHeatmap)

p4 <- heatmap(.data = cibersort_long
        ,.row = cell_type
        ,.column = ID
        ,.value = fraction
        ,scale = "column"
        ,palette_value = circlize::colorRamp2(
            seq(-2, 2, length.out = 11), 
            RColorBrewer::brewer.pal(11, "RdBu")
        )
        ,show_column_names=F
        ,row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 7),
        column_title_gp = gpar(fontsize = 7),
        row_title_gp = gpar(fontsize = 7)
        ) %>% 
  add_tile(Group) # 新版本已经改了，注意
p4
```

![plot of chunk unnamed-chunk-22](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-22-174361934.png)

`ssGSEA`的可视化方法也是类似的，不过就没有堆叠条形图了，因为加和不是1，堆叠起来就会参差不齐，毫无美感。

常用的还是箱线图这种样式的。


```r
ssgsea_long <- im_ssgsea %>% 
  pivot_longer(- ID,names_to = "cell_type",values_to = "Score")
head(ssgsea_long)
## # A tibble: 6 × 3
##   ID                           cell_type                         Score
##   <chr>                        <chr>                             <dbl>
## 1 TCGA-3L-AA1B-01A-11R-A37K-07 Activated B cell               -0.192  
## 2 TCGA-3L-AA1B-01A-11R-A37K-07 Activated CD4 T cell            0.0120 
## 3 TCGA-3L-AA1B-01A-11R-A37K-07 Activated CD8 T cell            0.170  
## 4 TCGA-3L-AA1B-01A-11R-A37K-07 Activated dendritic cell        0.00545
## 5 TCGA-3L-AA1B-01A-11R-A37K-07 CD56bright natural killer cell  0.199  
## 6 TCGA-3L-AA1B-01A-11R-A37K-07 CD56dim natural killer cell     0.351
```


```r
ggplot(ssgsea_long, aes(cell_type, Score))+
  geom_violin(width=2.0,aes(color=cell_type))+
  geom_boxplot(width=0.2,fill="black") + 
  theme_bw() + 
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  #scale_fill_manual(values = palette4)+
  scale_color_manual(values = palette4,name=NULL)
```

![plot of chunk unnamed-chunk-24](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-24-174361934.png)

除此之外也是可以添加分组，画热图等，其他的免疫浸润结果也是同样的可视化思路，这里就不再重复了，大家自己尝试下即可。

## 和分子联系起来

如果和某个分子联系起来，又可以画出各种花里胡哨的图，比如棒棒糖图，热图，散点图等。

我这里是以`ssGSEA`的结果为例进行演示的，其他的都是一样的。

我们就以`CTLA4` `EGFR` `PDL1`  `BRAF`  "VEGFB" "VEGFA" "VEGFC" "VEGFD" "NTRK2" "NTRK1" "NTRK3"这几个分子为例吧。


```r
genes <- c("CTLA4","EGFR","PDL1","BRAF","KRAS","VEGFA","VEGFB","VEGFC","VEGFD","NTRK1","NTRK2","NTRK3")

genes_expr <- as.data.frame(t(expr_coad[rownames(expr_coad) %in% genes,]))
genes_expr <- genes_expr[match(im_ssgsea$ID,rownames(genes_expr)),]
identical(im_ssgsea$ID,rownames(genes_expr))
## [1] TRUE
```

接下来就是批量计算每一个基因和28种细胞之间的相关系数和P值，这个需求你可以写循环实现，或者apply系列，purrr系列等，但是我试过，都太慢了，尤其是数据量很大的时候。

所以我这里给大家介绍一种更快的方法，借助`linkET`包实现，这个包在之前也介绍过了：[mantel-test可视化，别再只知道ggcor！](https://mp.weixin.qq.com/s/-2HPfesOhiUfHZDB25uDhg)


```r
library(linkET)
## 
## Attaching package: 'linkET'
## The following object is masked from 'package:purrr':
## 
##     simplify
## The following object is masked from 'package:ComplexHeatmap':
## 
##     anno_link

cor_res <- correlate(genes_expr, im_ssgsea[,-1],method = "spearman")
  
qcorrplot(cor_res) +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
```

![plot of chunk unnamed-chunk-26](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-26-174361934.png)

但是这个图并没有明显的表示出P值，所以我知道大家想自己画的更加花里胡哨一点，在很久之前我就介绍过了这个方法了：[R语言ggplot2画相关性热图](https://mp.weixin.qq.com/s/E5ev-Ltu-4oVnzk9gRF-NQ)

画图前先准备下数据，把P值数据和相关系数数据整合到一起，所以借助`linkET`包也是有缺点的，如果自己写函数肯定是直接弄好的。


```r
# 先整理下数据
df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "correlation")

df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`

head(df_cor)
## # A tibble: 6 × 5
##   gene  cell_type                      correlation   pvalue stars
##   <chr> <chr>                                <dbl>    <dbl> <fct>
## 1 VEGFB Activated B cell                    0.0837 5.55e- 2 " "  
## 2 VEGFB Activated CD4 T cell               -0.338  2.35e-15 "***"
## 3 VEGFB Activated CD8 T cell                0.0417 3.41e- 1 " "  
## 4 VEGFB Activated dendritic cell           -0.0181 6.79e- 1 " "  
## 5 VEGFB CD56bright natural killer cell      0.287  2.48e-11 "***"
## 6 VEGFB CD56dim natural killer cell         0.289  1.82e-11 "***"
```

数据都有了，画图就行了。`ggplot2`搞定一切，求求大家赶紧学学`ggplot2`吧，别再天天问图怎么画了。两本说明书，随便买一本认真看看就搞定了：《R数据可视化手册》或者《ggplot2：数据分析与图形艺术》


```r
library(ggplot2)

ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                      limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/df_cor.png)

```r
ggsave("df_cor.png",height = 4,width = 9) # 保存时不断调整长宽比即可得到满意的图形
```

棒棒糖图也是一样的简单，我们之前也介绍过了：[你还不会画棒棒糖图？](https://mp.weixin.qq.com/s/LP_C8OR62jvY4UCC_aryDw)

不过这种展示的是1个基因和其他细胞的关系，就是1对多的关系展示，上面的热图是多对多的关系展示。


```r
# 以EGFR为例
df_egfr <- df_cor %>% 
  filter(gene=="EGFR")

text_color <- c(rep("red",4),rep("#FDD819",1),rep("#67B26F",1),rep("#C4E0E5",12),rep("#67B26F",2),rep("#FDD819",3),rep("red",5))

ggplot(df_egfr, aes(correlation, fct_reorder(cell_type,correlation)))+
  geom_segment(aes(xend = 0,yend = cell_type),color="grey70",size=1.2)+ 
  geom_point(aes(size = abs(correlation), color=factor(stars)))+ 
  scale_color_manual(values = c("red","#FDD819","#67B26F","#C4E0E5"),name="P-value")+
  scale_size_continuous(range = c(3,8),name="Correlation")+ 
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(color="black",size = 14),
        axis.text.y = element_text(color = text_color,size = 14)#给标签上色没想到特别好的方法
        )
```

![plot of chunk unnamed-chunk-29](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-29-174361934.png)

```r
ggsave("df_egfr.png",width = 8,height = 8)
```

除此之外，你还可以可视化1个基因和1个细胞之间的关系，用散点图或者散点图矩阵的形式。

我们可以直接使用`ggplot2`里面的分面，画一张图。


```r
# 还是以EGFR为例
df_egfr_scatter <- im_ssgsea %>% 
  mutate(EGFR = genes_expr[,"EGFR"],.before = 1) %>% 
  pivot_longer(-c(1,2),names_to = "cell_type",values_to = "score")
```

画图即可：


```r
ggplot(df_egfr_scatter, aes(EGFR,score))+
  geom_point()+
  geom_smooth(method = "lm",color="blue")+
  stat_cor(method = "spearman",color="red")+
  facet_wrap(~cell_type,scales = "free_y",ncol = 5)
## `geom_smooth()` using formula = 'y ~ x'
```

![plot of chunk unnamed-chunk-31](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-31-174361934.png)

```r
ggsave("df_facet.png",width = 14,height = 8)
## `geom_smooth()` using formula = 'y ~ x'
```

>注意：到目前为止我们用的都是所有样本，tumor和normal都有！

我们也可以按照tumor和normal分个组，再画图：


```r
df_egfr_scatter$sample_type <- ifelse(as.numeric(substr(df_egfr_scatter$ID,14,15))<10,"tumor","normal")

ggplot(df_egfr_scatter, aes(EGFR,score,color=sample_type))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  scale_color_manual(values = c('#FDD819','#028EA1'))+
  facet_wrap(~cell_type,scales = "free_y",ncol = 5)
## `geom_smooth()` using formula = 'y ~ x'
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/df_facet_two.png)

```r
ggsave("df_facet_two.png",width = 14,height = 8)
## `geom_smooth()` using formula = 'y ~ x'
```

如果你需要单独画图，并保存，也是可以的，都没有问题，只要你基础够扎实，想做什么都可以，你R语言基础不行，啥都做不出来。


```r
# 还是以EGFR为例
df_egfr_scatter %>% 
  filter(cell_type == "Activated B cell") %>% 
  ggplot(aes(EGFR, score,color=sample_type))+
  geom_point()+
  geom_rug(aes(color=sample_type))+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  scale_color_manual(values = c('#FDD819','#028EA1'))+
  theme_bw()
## `geom_smooth()` using formula = 'y ~ x'
```

![plot of chunk unnamed-chunk-33](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-33-174361934.png)

是不是很简单呢？

然后你可以循环出图并保存到本地，不过我并没有使用上面这种花里胡哨的图，你可以自己修改：


```r
library(purrr)

plot_list <- df_egfr_scatter %>% 
  split(.$cell_type) %>% # group_split没有名字
  map(~ ggplot(., aes(EGFR,score,color=sample_type))+
        geom_point()+
        geom_smooth(method = "lm")+
        stat_cor(method = "spearman")+
        scale_color_manual(values = c('#FDD819','#028EA1'))
      )
paths <- paste0(names(plot_list),".png")

pwalk(list(paths, plot_list),ggsave,width=6,height=3)
```

28张图片已保存到本地：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230509182427698.png)

每一张都长这样：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230509182545366.png)

OK，就先介绍到这里，关于结果的可视化，我这里介绍的只是最常见的，冰山一角而已，毕竟可视化方法太多了，不可能全都介绍到。

大家如果有喜欢的图形，可通过评论区，粉丝QQ群等方式发给我~



