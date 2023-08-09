2022年都快结束了！我终于要学习单细胞分析了！

曾老师的单细胞视频看好几遍了，但是关于单细胞的代码一句也没写过。。。主要原因是一直没有这方面的需求。。

眼看着2022年就快结束了，再不学就来不及了！



如果你没有常规转录组分析的基础，有一些概念可能会难以理解，这时候我推荐你去B站看一看生信技能树的相关视频，看看不吃亏，反正是免费的。

关于单细胞的技术原理以及主要应用方向之类的，大家可以自行学习，这里就不赘述了，我们直接进入单细胞的下游分析标准流程！

![单细胞分析常见流程](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/Snipaste_2022-08-11_20-38-25.png)

下面是我的学习笔记。

目录：

[toc]

## 安装Seurat

Seurat可以说是目前单细胞下游分析最常用的R包之一了，单细胞必备R包。


```r
# 竟然不是bioconductor！
install.packages('Seurat')
library(Seurat)
```

## 建立Seurat对象

读取表达数据，这里使用PBMC，先去下载数据：https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz，然后解压。

这个数据是人类外周血单核细胞的10X单细胞数据，使用NextSeq 500测了2700个细胞。

类似于普通转录组的表达矩阵，一个样本可以测多个细胞，每个细胞里有多个基因。


```r
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(Seurat)
## Attaching SeuratObject
## Attaching sp
library(patchwork)
library(ggplot2)

pbmc.data <- Read10X("../000files/filtered_gene_bc_matrices/hg19/") 
dim(pbmc.data)
## [1] 32738  2700
```

这就是一个稀疏矩阵，也是一个表达矩阵，行是基因，列是细胞。

可以随便取子集：


```r
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
## 3 x 30 sparse Matrix of class "dgCMatrix"
##    [[ suppressing 30 column names 'AAACATACAACCAC-1', 'AAACATTGAGCTAC-1', 'AAACATTGATCAGC-1' ... ]]
##                                                                    
## CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
## TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
## MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .
```

可以看到有很多`.`，代表0，啥都没测到，单细胞测序中有很多0，使用稀疏矩阵可以节省内存、提高速度。和普通的转录组表达矩阵不太一样。

接下来是单细胞分析的最基本步骤：构建`seurat`对象，很重要！下面所有分析都是基于这个对象的。

众所周知，学会单细胞，你将拥有很多**对象**！


```r
pbmc <- CreateSeuratObject(counts = pbmc.data,
                           project = "pbmc3k",
                           min.cells = 3, # 基因至少在多少个细胞中有表达
                           min.features = 200 # 至少表达多少基因
                           )
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
```

这个对象里面包含13714行（基因），2700列（细胞）。

## 标准预处理流程

质量控制，标准化和归一化，挑选高变基因，降维，聚类，分群。

### 质量控制，挑选细胞

主要就是根据`nfeature\nCount\线粒体基因比例`挑选合适的细胞。

- 每个细胞中测到的基因太少或太多的（nfeature），
- 每个细胞内检测到的分子总数nCount(UMI总数)太多或太少的，
- 线粒体基因过多的，

这些细胞都不符合要求，需要去除。


```r
# 计算线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 查看前6个细胞的QC指标
head(pbmc@meta.data, 6)
##                  orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
## AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
## AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363
## AAACCGTGCTTCCG-1     pbmc3k       2639          960  1.7430845
## AAACCGTGTATGCG-1     pbmc3k        980          521  1.2244898
## AAACGCACTGGTAC-1     pbmc3k       2163          781  1.6643551
```

- nCount_RNA：每个细胞的UMI数
- nFeature_RNA：每个细胞的基因数
- 一个基因可以有多个UMI。

用图形的方式查看QC结果：


```r
# 看看结果
VlnPlot(pbmc,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"), #纵坐标
        ncol = 3
        )
```

![plot of chunk unnamed-chunk-6](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-6-150864622.png)

可以看出大部分细胞的`nFeature_RNA`都在2000以下，`nCount_RNA`大部分都在6000以下，大部分细胞的线粒体基因比例在5%以下。


```r
# 任选两个指标，看看两者之间的关系
plot1 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2
```

![plot of chunk unnamed-chunk-7](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-7-150864622.png)

`nCount_RNA`和`nFeature_RNA`非常一致哦~，相关性0.95，这些图其实可以自己用`ggplot2`画。


```r
# 根据上面的QC结果挑选合适的细胞
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc
## An object of class Seurat 
## 13714 features across 2638 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
```

还剩13714个基因，2638个细胞。

### 标准化


```r
pbmc <- NormalizeData(pbmc)
```


### 挑选高变基因

挑选在不同细胞中高度变化的基因，比如某些基因在一些细胞中高表达，在另一些基因中低表达，把这些基因找出来继续进行下一步分析。


```r
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures = 2000 # 选择变化最大的前2000个基因
                             )

# 看看前10个
top10 <- head(VariableFeatures(pbmc),10)
top10
##  [1] "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"  
##  [9] "GNG11"  "S100A8"
```

图形方式查看结果：


```r
# 均值-方差图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1,points = top10,repel = T) # 给前10个加名字
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1+plot2
```

![plot of chunk unnamed-chunk-11](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-11-150864622.png)

### 缩放

主要有以下作用： 
- 改变基因表达量，使细胞间表达量为0
- 改变基因表达量，使细胞间方差是1
  - 这一步可以让高表达基因不会占太多主导地位


```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) # 对所有基因进行scale，也可以只选择高变基因
## Centering and scaling data matrix
```

### PCA降维

进行PCA降维：


```r
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
## PC_ 1 
## Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
## 	   FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
## 	   PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
## Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
## 	   CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
## 	   MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
## PC_ 2 
## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
## 	   HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
## 	   BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
## Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
## 	   CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
## 	   TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
## PC_ 3 
## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
## 	   HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
## 	   PLAC8, BLNK, MALAT1, SMIM14, PLD4, LAT2, IGLL5, P2RX5, SWAP70, FCGR2B 
## Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
## 	   HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
## 	   NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
## PC_ 4 
## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HLA-DPB1, HIST1H2AC, PF4, TCL1A 
## 	   SDPR, HLA-DPA1, HLA-DRB1, HLA-DQA2, HLA-DRA, PPBP, LINC00926, GNG11, HLA-DRB5, SPARC 
## 	   GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, GZMB, CLU, TUBB1 
## Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
## 	   AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
## 	   LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
## PC_ 5 
## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
## 	   GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, CCL5 
## 	   RBP7, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
## Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
## 	   PTGES3, LILRB2, MAL, CD27, HN1, CD2, GDI2, ANXA5, CORO1B, TUBA1B 
## 	   FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB
```


```r
# 查看PCA结果
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  CST3, TYROBP, LST1, AIF1, FTL 
## Negative:  MALAT1, LTB, IL32, IL7R, CD2 
## PC_ 2 
## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
## Negative:  NKG7, PRF1, CST7, GZMB, GZMA 
## PC_ 3 
## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1 
## Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
## PC_ 4 
## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
## Negative:  VIM, IL7R, S100A6, IL32, S100A8 
## PC_ 5 
## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY 
## Negative:  LTB, IL7R, CKB, VIM, MS4A7
```

下面是用各种图形的方式展示PCA结果：


```r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```

![plot of chunk unnamed-chunk-15](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-15-150864622.png)


```r
DimPlot(pbmc, reduction = "pca")
```

![plot of chunk unnamed-chunk-16](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-16-150864622.png)


```r
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```

![plot of chunk unnamed-chunk-17](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-17-150864622.png)


```r
# 展示多个维度，可以看到越往后越不好
DimHeatmap(pbmc,dims = 1:15,cells = 500,balanced = T)
```

![plot of chunk unnamed-chunk-18](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-18-150864622.png)

### 确定挑选哪几个维度


```r
# 挑选合适的维度进行下一步分析，这一步很浪费时间
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
```

用图形的方式查看结果，显示每个成分的P值分布，一般越主要的维度P值越小。


```r
JackStrawPlot(pbmc,dims = 1:15) # 展示前15个维度
```

![plot of chunk unnamed-chunk-20](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-20-150864622.png)

另一种方法：


```r
ElbowPlot(pbmc) # 更节省时间的方法
```

![plot of chunk unnamed-chunk-21](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-21-150864622.png)

可以看到在10个主成分以后基本就是平的了，没什么变化了。

>官方建议使用多个不同的维度进行下游分析试试看，选择维度不要太小！

### 聚类

- resolution参数控制聚类分析的粒度，越大细胞类群越多。
- 3000左右的细胞量，设置resolution为0.4-1.2是比较合适的。
- 细胞数据集越大，需要的resolution越大, 也会获得更多的细胞聚类。


```r
pbmc <- FindNeighbors(pbmc, dims = 1:10) # 前10个维度
## Computing nearest neighbor graph
## Computing SNN
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95965
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8723
## Number of communities: 9
## Elapsed time: 0 seconds
```


```r
# 查看前5个细胞的类群ID
head(Idents(pbmc), 5)
## AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
##                2                3                2                1 
## AAACCGTGTATGCG-1 
##                6 
## Levels: 0 1 2 3 4 5 6 7 8
```

### 非线性降维方法 (UMAP/tSNE)分群


```r
pbmc <- RunUMAP(pbmc,dims = 1:10)
## 17:02:54 UMAP embedding parameters a = 0.9922 b = 1.112
## 17:02:54 Read 2638 rows and found 10 numeric columns
## 17:02:54 Using Annoy for neighbor search, n_neighbors = 30
## 17:02:54 Building Annoy index with metric = cosine, n_trees = 50
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## 17:02:54 Writing NN index file to temp file C:\Users\liyue\AppData\Local\Temp\RtmpoXaiZI\file3a9c15526d0a
## 17:02:54 Searching Annoy index using 1 thread, search_k = 3000
## 17:02:54 Annoy recall = 100%
## 17:02:54 Commencing smooth kNN distance calibration using 1 thread
## 17:02:55 Initializing from normalized Laplacian + noise
## 17:02:55 Commencing optimization for 500 epochs, with 105124 positive edges
## 17:03:01 Optimization finished
```

可视化结果，比PCA好多了！


```r
DimPlot(pbmc, reduction = "umap",label = T)
```

![plot of chunk unnamed-chunk-25](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-25-150864622.png)

把结果保存下，方便下次使用。


```r
saveRDS(pbmc, file = "../000files/pbmc_tutorial.rds")
```

基本上所有的单细胞分析都要经过上面的标准分析流程，得到现在的这个`seurat`对象才是接下来一切分析的基础。

所以，**把中间可视化的步骤去掉后，基本的流程就是下面这几步：**

```{r,eval=FALSE}
pbmc.data <- Read10X("../000files/filtered_gene_bc_matrices/hg19/") 
pbmc <- CreateSeuratObject(counts = pbmc.data)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc,dims = 1:10)
pbmc <- RunTSNE(pbmc,dims = 1:10)
```

下面就是差异分析、亚群注释、细胞轨迹分析、细胞通讯、转录因子分析、拷贝数变异等。

### 挑选差异基因

`FindMarkers()`默认会找单个类群和其他所有类群相比的marker基因。


```r
# 找类群2的marker基因（和其他所有类群比）
cluster2.markers <- FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25) 
head(cluster2.markers,5)
##             p_val avg_log2FC pct.1 pct.2    p_val_adj
## IL32 2.593535e-91  1.2154360 0.949 0.466 3.556774e-87
## LTB  7.994465e-87  1.2828597 0.981 0.644 1.096361e-82
## CD3D 3.922451e-70  0.9359210 0.922 0.433 5.379250e-66
## IL7R 1.130870e-66  1.1776027 0.748 0.327 1.550876e-62
## LDHB 4.082189e-65  0.8837324 0.953 0.614 5.598314e-61
```


```r
# 找类群5的marker基因（和类群0及类群3比）
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
## FCGR3A        2.150929e-209   4.267579 0.975 0.039 2.949784e-205
## IFITM3        6.103366e-199   3.877105 0.975 0.048 8.370156e-195
## CFD           8.891428e-198   3.411039 0.938 0.037 1.219370e-193
## CD68          2.374425e-194   3.014535 0.926 0.035 3.256286e-190
## RP11-290F20.3 9.308287e-191   2.722684 0.840 0.016 1.276538e-186
```


```r
# 一步到位，找到所有类群的marker基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, # 保留表达上调的marker
                               min.pct = 0.25, # 在两个细胞群中任何一个被检测到的百分比，默认0.1
                               test.use = "wilcox", # wlicox,bimod,roc,Student`s t test,poisson,negbinom,LR,MAST,DEseq2
                               logfc.threshold = 0.25)
## Calculating cluster 0
## Calculating cluster 1
## Calculating cluster 2
## Calculating cluster 3
## Calculating cluster 4
## Calculating cluster 5
## Calculating cluster 6
## Calculating cluster 7
## Calculating cluster 8

# 查看每个类群中前2个marker基因
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
## # A tibble: 18 × 7
## # Groups:   cluster [9]
##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
##  1 1.17e- 83       1.33 0.435 0.108 1.60e- 79 0       CCR7    
##  2 1.74e-109       1.07 0.897 0.593 2.39e-105 0       LDHB    
##  3 0               5.57 0.996 0.215 0         1       S100A9  
##  4 0               5.48 0.975 0.121 0         1       S100A8  
##  5 7.99e- 87       1.28 0.981 0.644 1.10e- 82 2       LTB     
##  6 2.61e- 59       1.24 0.424 0.111 3.58e- 55 2       AQP3    
##  7 0               4.31 0.936 0.041 0         3       CD79A   
##  8 9.48e-271       3.59 0.622 0.022 1.30e-266 3       TCL1A   
##  9 4.93e-169       3.01 0.595 0.056 6.76e-165 4       GZMK    
## 10 1.17e-178       2.97 0.957 0.241 1.60e-174 4       CCL5    
## 11 3.51e-184       3.31 0.975 0.134 4.82e-180 5       FCGR3A  
## 12 2.03e-125       3.09 1     0.315 2.78e-121 5       LST1    
## 13 6.82e-175       4.92 0.958 0.135 9.36e-171 6       GNLY    
## 14 1.05e-265       4.89 0.986 0.071 1.44e-261 6       GZMB    
## 15 1.48e-220       3.87 0.812 0.011 2.03e-216 7       FCER1A  
## 16 1.67e- 21       2.87 1     0.513 2.28e- 17 7       HLA-DPB1
## 17 3.68e-110       8.58 1     0.024 5.05e-106 8       PPBP    
## 18 7.73e-200       7.24 1     0.01  1.06e-195 8       PF4
```



可视化结果：


```r
# 所选基因在不同类群中的表达量
VlnPlot(pbmc,features = c("MS4A1","CD79A"))
```

![plot of chunk unnamed-chunk-32](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-32-150864622.png)



纵坐标做变为log counts：


```r
VlnPlot(pbmc, features = c("NKG7", "PF4"),
        slot = "counts", log = TRUE)
```

![plot of chunk unnamed-chunk-33](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-33-150864622.png)


```r
# 可以看到每个基因都是在某个亚群中高表达
FeaturePlot(pbmc, 
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
```

![plot of chunk unnamed-chunk-34](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-34-150864622.png)



热图的形式展示每个亚群中表达量前10的基因：


```r
top10 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

# 展示每个类群前10个基因
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

![plot of chunk unnamed-chunk-35](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-35-150864622.png)



以气泡图的形式展示每个亚群中表达量较高的基因：


```r
# 这里用20个基因展示
DotPlot(pbmc, features = unique(top10$gene)[1:20])+RotatedAxis()
```

![image-20220815125150676](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220815125150676.png)

### 细胞类群注释

不同的细胞类型会特异性表达一些基因，于是就可以根据某些特定基因的表达来确定这一类群是属于哪类细胞。

可以用R包，也可以用自己找的，比如cellmarker网站上面的数据。

现在亚群的名字还是数字，我们改一下：


```r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet") # 自定义类群名字

names(new.cluster.ids) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.ids) # 重命名类群

# 重新画图
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![plot of chunk unnamed-chunk-38](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-38-150864622.png)

保存下数据。


```r
saveRDS(pbmc, file = "../000files/pbmc3k_final.rds")
```


## 参考资料

- Seurat官网
- 生信技能树、单细胞天地
- 菲沙基因单细胞培训课程
