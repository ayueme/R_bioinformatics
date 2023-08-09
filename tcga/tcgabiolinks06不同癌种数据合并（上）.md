
很多文章对于TCGA中的一些癌症都是联合分析的，比如TCGA-COAD和TCGA-READ，首先是它们的疾病特点和治疗方式存在很多相似之处，同时这样做也可以增大样本量。

如果你是使用`TCGAbiolinks`包下载的数据，那么它们的合并超级简单，直接`cbind()`即可！

## 加载数据和R包

数据都是之前下载好的，可以参考之前的推文：

[1.新版TCGA数据库学习：批量下载数据](https://mp.weixin.qq.com/s/m8w1L4N2aXAIers_ZJvp_g)

[2.新版TCGA数据库学习：表达矩阵提取（mRNA/lncRNA/counts/tpm/fpkm）](https://mp.weixin.qq.com/s/wI0_GyVl5LiKAjX5C3f-NQ)

[3.手动下载的TCGA数据也是可以用TCGAbiolinks包整理的](https://mp.weixin.qq.com/s/DHj9wp6hkae2Zrl61sU1fQ)

我们直接加载TCGA-COAD和TCGA-READ的数据。


```r
#library(TCGAbiolinks)

# COAD
load(file = "./TCGA-mRNA/TCGA-COAD_mRNA.Rdata")
coad <- data

# READ
load(file = "./TCGA-mRNA/TCGA-READ_mRNA.Rdata")
read <- data
```

## 合并数据

现在`coad`和`read`都是`SummarizedExperiment`对象，并且具有相同的行和行名：


```r
coad
## Loading required package: SummarizedExperiment
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## Loading required package: GenomicRanges
## Loading required package: stats4
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following object is masked from 'package:grDevices':
## 
##     windows
## Loading required package: GenomeInfoDb
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
## class: RangedSummarizedExperiment 
## dim: 60660 521 
## metadata(1): data_release
## assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(521): TCGA-A6-5664-01A-21R-1839-07
##   TCGA-D5-6530-01A-11R-1723-07 ... TCGA-A6-2683-01A-01R-0821-07
##   TCGA-A6-2683-11A-01R-A32Z-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status

read
## class: RangedSummarizedExperiment 
## dim: 60660 177 
## metadata(1): data_release
## assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(177): TCGA-AG-3580-01A-01R-0821-07
##   TCGA-AF-2692-11A-01R-A32Z-07 ... TCGA-AG-3894-01A-01R-1119-07
##   TCGA-AG-3574-01A-01R-0821-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

对于这样的数据我们直接合并即可，我认为这是目前合并两个癌种最方便的方法了！


```r
# 直接cbind
colrectal <- cbind(coad,read)

colrectal
## class: RangedSummarizedExperiment 
## dim: 60660 698 
## metadata(2): data_release data_release
## assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
## rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
##   ENSG00000288674.1 ENSG00000288675.1
## rowData names(10): source type ... hgnc_id havana_gene
## colnames(698): TCGA-A6-5664-01A-21R-1839-07
##   TCGA-D5-6530-01A-11R-1723-07 ... TCGA-AG-3894-01A-01R-1119-07
##   TCGA-AG-3574-01A-01R-0821-07
## colData names(107): barcode patient ... paper_vascular_invasion_present
##   paper_vital_status
```

得到的结果也是一个`SummarizedExperiment`对象。并且这个对象中各种信息也是保存好的，想用什么直接提取即可，非常方便。

但是这样合并可能涉及批次效应的问题，大家在实际使用时可根据自己的情况选择要不要去除批次效应！

## 提取信息

比如提取样本的临床信息，非常简单，甚至不需要重新下载：


```r
clin <- as.data.frame(colData(colrectal))

clin[1:10,1:10]
##                                                   barcode      patient
## TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664
## TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530
## TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556
## TCGA-AA-3660-11A-01R-1723-07 TCGA-AA-3660-11A-01R-1723-07 TCGA-AA-3660
## TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818
## TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660
## TCGA-DM-A28G-01A-11R-A16W-07 TCGA-DM-A28G-01A-11R-A16W-07 TCGA-DM-A28G
## TCGA-AA-3976-01A-01R-1022-07 TCGA-AA-3976-01A-01R-1022-07 TCGA-AA-3976
## TCGA-G4-6307-01A-11R-1723-07 TCGA-G4-6307-01A-11R-1723-07 TCGA-G4-6307
## TCGA-AA-3522-11A-01R-A32Z-07 TCGA-AA-3522-11A-01R-A32Z-07 TCGA-AA-3522
##                                        sample shortLetterCode
## TCGA-A6-5664-01A-21R-1839-07 TCGA-A6-5664-01A              TP
## TCGA-D5-6530-01A-11R-1723-07 TCGA-D5-6530-01A              TP
## TCGA-AA-3556-01A-01R-0821-07 TCGA-AA-3556-01A              TP
## TCGA-AA-3660-11A-01R-1723-07 TCGA-AA-3660-11A              NT
## TCGA-AA-3818-01A-01R-0905-07 TCGA-AA-3818-01A              TP
## TCGA-AA-3660-01A-01R-1723-07 TCGA-AA-3660-01A              TP
## TCGA-DM-A28G-01A-11R-A16W-07 TCGA-DM-A28G-01A              TP
## TCGA-AA-3976-01A-01R-1022-07 TCGA-AA-3976-01A              TP
## TCGA-G4-6307-01A-11R-1723-07 TCGA-G4-6307-01A              TP
## TCGA-AA-3522-11A-01R-A32Z-07 TCGA-AA-3522-11A              NT
##                                       definition sample_submitter_id
## TCGA-A6-5664-01A-21R-1839-07 Primary solid Tumor    TCGA-A6-5664-01A
## TCGA-D5-6530-01A-11R-1723-07 Primary solid Tumor    TCGA-D5-6530-01A
## TCGA-AA-3556-01A-01R-0821-07 Primary solid Tumor    TCGA-AA-3556-01A
## TCGA-AA-3660-11A-01R-1723-07 Solid Tissue Normal    TCGA-AA-3660-11A
## TCGA-AA-3818-01A-01R-0905-07 Primary solid Tumor    TCGA-AA-3818-01A
## TCGA-AA-3660-01A-01R-1723-07 Primary solid Tumor    TCGA-AA-3660-01A
## TCGA-DM-A28G-01A-11R-A16W-07 Primary solid Tumor    TCGA-DM-A28G-01A
## TCGA-AA-3976-01A-01R-1022-07 Primary solid Tumor    TCGA-AA-3976-01A
## TCGA-G4-6307-01A-11R-1723-07 Primary solid Tumor    TCGA-G4-6307-01A
## TCGA-AA-3522-11A-01R-A32Z-07 Solid Tissue Normal    TCGA-AA-3522-11A
##                              sample_type_id
## TCGA-A6-5664-01A-21R-1839-07             01
## TCGA-D5-6530-01A-11R-1723-07             01
## TCGA-AA-3556-01A-01R-0821-07             01
## TCGA-AA-3660-11A-01R-1723-07             11
## TCGA-AA-3818-01A-01R-0905-07             01
## TCGA-AA-3660-01A-01R-1723-07             01
## TCGA-DM-A28G-01A-11R-A16W-07             01
## TCGA-AA-3976-01A-01R-1022-07             01
## TCGA-G4-6307-01A-11R-1723-07             01
## TCGA-AA-3522-11A-01R-A32Z-07             11
##                                                         sample_id
## TCGA-A6-5664-01A-21R-1839-07 3048539a-b914-4e43-b1cc-43ea707e3b3d
## TCGA-D5-6530-01A-11R-1723-07 50560725-c72d-4bab-b602-5e50e6bececd
## TCGA-AA-3556-01A-01R-0821-07 4794413c-ed92-451c-a3ce-f411fed5ca82
## TCGA-AA-3660-11A-01R-1723-07 a0832917-75b9-45c8-9273-009c3737a43a
## TCGA-AA-3818-01A-01R-0905-07 0cf35153-2c04-4bdd-91e0-d63cb98da5bf
## TCGA-AA-3660-01A-01R-1723-07 87cf1a20-2dc5-4c06-b0c4-16103be40ef0
## TCGA-DM-A28G-01A-11R-A16W-07 f5acd8b8-32c4-4f3c-aacd-259f8e1fdfee
## TCGA-AA-3976-01A-01R-1022-07 8d529023-abca-4ddc-a265-d5bd1fd48708
## TCGA-G4-6307-01A-11R-1723-07 801b8d05-2d29-4f8c-8d8d-634b2e21b867
## TCGA-AA-3522-11A-01R-A32Z-07 218bbd07-5fa3-4946-a2c1-0ece13466441
##                                      sample_type days_to_collection
## TCGA-A6-5664-01A-21R-1839-07       Primary Tumor                 NA
## TCGA-D5-6530-01A-11R-1723-07       Primary Tumor                 NA
## TCGA-AA-3556-01A-01R-0821-07       Primary Tumor                 NA
## TCGA-AA-3660-11A-01R-1723-07 Solid Tissue Normal                 NA
## TCGA-AA-3818-01A-01R-0905-07       Primary Tumor                 NA
## TCGA-AA-3660-01A-01R-1723-07       Primary Tumor                 NA
## TCGA-DM-A28G-01A-11R-A16W-07       Primary Tumor               3419
## TCGA-AA-3976-01A-01R-1022-07       Primary Tumor                 NA
## TCGA-G4-6307-01A-11R-1723-07       Primary Tumor                 NA
## TCGA-AA-3522-11A-01R-A32Z-07 Solid Tissue Normal                 NA

dim(clin)
## [1] 698 107

colnames(clin)[10:30]
##  [1] "days_to_collection"        "state"                    
##  [3] "initial_weight"            "intermediate_dimension"   
##  [5] "pathology_report_uuid"     "submitter_id"             
##  [7] "shortest_dimension"        "oct_embedded"             
##  [9] "longest_dimension"         "is_ffpe"                  
## [11] "tissue_type"               "synchronous_malignancy"   
## [13] "ajcc_pathologic_stage"     "days_to_diagnosis"        
## [15] "treatments"                "last_known_disease_status"
## [17] "tissue_or_organ_of_origin" "days_to_last_follow_up"   
## [19] "age_at_diagnosis"          "primary_diagnosis"        
## [21] "prior_malignancy"
```

现在一共有698行，107列临床信息，**你想要的生存时间、生存状态、样本类型、分期等信息都在里面，都不需要自己手动划分，想要什么直接取子集就好了。**


比如大家最喜欢的生存信息：


```r
clin_subset <- clin[,c("days_to_last_follow_up","vital_status")]

head(clin_subset)
##                              days_to_last_follow_up vital_status
## TCGA-A6-5664-01A-21R-1839-07                    672        Alive
## TCGA-D5-6530-01A-11R-1723-07                    621        Alive
## TCGA-AA-3556-01A-01R-0821-07                    700        Alive
## TCGA-AA-3660-11A-01R-1723-07                   2375        Alive
## TCGA-AA-3818-01A-01R-0905-07                     NA         Dead
## TCGA-AA-3660-01A-01R-1723-07                   2375        Alive
```

## 合并miRNA

也是一样的操作。


```r
rm(list = ls())

load(file = "./TCGA-mirna/TCGA-COAD_miRNA.Rdata")
coad <- data

load(file = "./TCGA-mirna/TCGA-READ_miRNA.Rdata")
read <- data
```

可以看到两个表达矩阵的第一列（miRNA的名字），完全一样：


```r
identical(coad$miRNA_ID,read$miRNA_ID)
## [1] TRUE
```

所以我们直接合并即可：


```r
# 第一列都是
colrectal_mi <- cbind(coad,read[,-1])

colrectal_mi[1:5,1:4]
##       miRNA_ID read_count_TCGA-A6-5664-01A-21H-1838-13
## 1 hsa-let-7a-1                                    6959
## 2 hsa-let-7a-2                                    6941
## 3 hsa-let-7a-3                                    7120
## 4   hsa-let-7b                                   31616
## 5   hsa-let-7c                                    4211
##   reads_per_million_miRNA_mapped_TCGA-A6-5664-01A-21H-1838-13
## 1                                                    8143.201
## 2                                                    8122.137
## 3                                                    8331.598
## 4                                                   36996.038
## 5                                                    4927.578
##   cross-mapped_TCGA-A6-5664-01A-21H-1838-13
## 1                                         N
## 2                                         N
## 3                                         N
## 4                                         N
## 5                                         N
```

但是miRNA的表达矩阵现在还有点问题，它包含3种信息：count/rpm/cross-mapped，而我们只需要count，所以还是要处理一下。


```r
dim(colrectal_mi)
## [1] 1881 1891

# 只要count
colrec_mi <- colrectal_mi[,c(1,seq(2,1891,by=3))]
dim(colrec_mi)
## [1] 1881  631

# 改下列名
colnames(colrec_mi)[-1] <- substr(colnames(colrec_mi)[-1],12,39)

colrec_mi[1:5,1:5]
##       miRNA_ID TCGA-A6-5664-01A-21H-1838-13 TCGA-A6-2683-01A-01T-0822-13
## 1 hsa-let-7a-1                         6959                        50288
## 2 hsa-let-7a-2                         6941                        50537
## 3 hsa-let-7a-3                         7120                        51098
## 4   hsa-let-7b                        31616                       143822
## 5   hsa-let-7c                         4211                         3943
##   TCGA-D5-6530-01A-11H-1722-13 TCGA-DM-A28G-01A-11H-A16S-13
## 1                        35778                        11788
## 2                        35334                        11588
## 3                        35980                        11885
## 4                        68674                        12086
## 5                          605                         1171
```

简单！

## 合并CNV


```r
rm(list = ls())
load("G:/tcga/TCGA-CNV/TCGA-COAD_CNV.Rdata")
coad <- data

load("G:/tcga/TCGA-CNV/TCGA-READ_CNV.Rdata")
read <- data

colrec_cnv <- rbind(coad,read)

head(colrec_cnv)
##                            GDC_Aliquot Chromosome    Start       End Num_Probes
## 1 741d4882-3a2c-4862-8402-636e6aebfdc6          1  3301765 247650984     129758
## 2 741d4882-3a2c-4862-8402-636e6aebfdc6          2   480597 241537572     132218
## 3 741d4882-3a2c-4862-8402-636e6aebfdc6          3  2170634  25586863      14093
## 4 741d4882-3a2c-4862-8402-636e6aebfdc6          3 25587626  25587698          3
## 5 741d4882-3a2c-4862-8402-636e6aebfdc6          3 25588064 197812401      93106
## 6 741d4882-3a2c-4862-8402-636e6aebfdc6          4  1059384 124538250      69561
##   Segment_Mean                       Sample
## 1      -0.0019 TCGA-AA-3556-10A-01D-0819-01
## 2      -0.0007 TCGA-AA-3556-10A-01D-0819-01
## 3      -0.0009 TCGA-AA-3556-10A-01D-0819-01
## 4      -1.9166 TCGA-AA-3556-10A-01D-0819-01
## 5       0.0023 TCGA-AA-3556-10A-01D-0819-01
## 6       0.0025 TCGA-AA-3556-10A-01D-0819-01
```


这个文件稍加整理就可以拿去给gistic用了。



## 合并SNP


```r
rm(list = ls())

load("G:/tcga/TCGA-SNP/TCGA-READ_SNP.Rdata")
read <- data

load("G:/tcga/TCGA-SNP/TCGA-COAD_SNP.Rdata")
coad <- data

colrec_snp <- rbind(coad,read)
```




这样以后再分析就可以用合并后的数据了！
