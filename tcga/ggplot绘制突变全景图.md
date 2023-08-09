这篇是生信技能树的一个学徒作业：小队列的肿瘤外显子临床预后意义

![image-20220903185149650](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185149650.png)

主要学习的图是这几个：

![突变全景图](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185248626.png)

![fig2a](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185329843.png)

![fig2c](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185352286.png)

## 读取数据

附件下载地址：https://ehoonline.biomedcentral.com/articles/10.1186/s40164-021-00200-x


```r
s2 <- data.table::fread("./s240164_2021_200_MOESM2_ESM.csv",skip = 2,header = T)
```

加载R包：


```r
suppressMessages(library(tidyverse))
```

```
## Warning: package 'tidyverse' was built under R version 4.2.1
```

改个名字方便使用：


```r
names(s2)[3] <- "gene_symbol"
names(s2)
```

```
##  [1] "SampleID"    "PatientID"   "gene_symbol" "cHGVS"       "pHGVS"      
##  [6] "Function"    "Transcript"  "ExIn_ID"     "Cosmic ID"   "Vary Type"  
## [11] "caseAF"
```

## Fig2a

Fig2a其实就是突变全景图的右边条形图部分，但是作者给分开展示了。只要计算某个基因在多少个样本中突变了，再除以53即可得到纵坐标mutation percentage！


```r
# 把这几个基因挑出来
mutationper <- s2 |> 
  dplyr::filter(gene_symbol %in% c("TET2","RHOA","PCLO",
                                   "DNMT3A","IDH2","PIEZO1",
                                   "TP53","RELN","FAT3",
                                   "CHD3")) |> 
  group_by(PatientID,gene_symbol) |> 
  summarise(n=n()) |> 
  mutate(n = case_when(n > 3 ~ "3",
                       T ~ as.character(n)
                       )) |> 
  ungroup() |> 
  group_by(gene_symbol,n) |> 
  summarise(nsub=n(),per=nsub/53) |> 
  mutate(gene_symbol = factor(gene_symbol,
                              levels = c("TET2","RHOA","PCLO",
                                   "DNMT3A","IDH2","PIEZO1",
                                   "TP53","RELN","FAT3",
                                   "CHD3")),
         n = factor(n,levels = c("3","2","1"))
         )
```

```
## `summarise()` has grouped output by 'PatientID'. You can override
## using the `.groups` argument.
## `summarise()` has grouped output by 'gene_symbol'. You can override
## using the `.groups` argument.
```

```r
mutationper
```

```
## # A tibble: 20 × 4
## # Groups:   gene_symbol [10]
##    gene_symbol n      nsub    per
##    <fct>       <fct> <int>  <dbl>
##  1 CHD3        1         6 0.113 
##  2 CHD3        2         1 0.0189
##  3 CHD3        3         1 0.0189
##  4 DNMT3A      1         9 0.170 
##  5 DNMT3A      2         1 0.0189
##  6 FAT3        1         6 0.113 
##  7 FAT3        2         2 0.0377
##  8 IDH2        1         8 0.151 
##  9 IDH2        2         1 0.0189
## 10 PCLO        1         9 0.170 
## 11 PCLO        2         3 0.0566
## 12 PIEZO1      1         8 0.151 
## 13 PIEZO1      2         1 0.0189
## 14 RELN        1         7 0.132 
## 15 RELN        2         1 0.0189
## 16 RHOA        1        23 0.434 
## 17 TET2        1         7 0.132 
## 18 TET2        2        17 0.321 
## 19 TET2        3        10 0.189 
## 20 TP53        1         8 0.151
```

画图即可：


```r
mutationper |> 
  ggplot(aes(gene_symbol,per))+
  geom_bar(stat = "identity",aes(fill=n))+
  scale_x_discrete(name=NULL)+
  scale_y_continuous(name="Mutation Percentage %",expand = c(0,0),
                     limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8),
                     labels = c(0,20,40,60,80)
                     )+
  scale_fill_discrete(name=NULL,labels=c("≥ 3 mutations",
                                         "2 mutations","1 mutation"
                                         ))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,color = "black",size = 12,hjust = 1),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.y = element_text(color = "black",size = 14),
        legend.position = c(0.8,0.6),
        axis.line = element_line(color = "black",size = 1.1),
        axis.ticks = element_line(color = "black",size = 1.1)
        
        )
```

![fig2a](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185604214.png)


## Fig1

最开始想用`complexheatmap`画，但是发现是长数据，可以直接用`ggplot2`画。


```r
# 预处理数据

heat_df <- s2

# 挑选展示的基因
genes <- c("CHD3","APC","TP53","PALB2","FANCA","TET2","DNMT3A","IDH2","ARID1A","ARID1B","MLL3","TYK2","STAT3","LRRK2","MAP2K1","PCLO","PIEZO1","FAT3","CSMD1","NSD1","MKI67","WDR90","MGA","CPS1","SPEN","ATP10B","ANKRD11","RELN","PLCG1","ALK","FLT4","RHOA","NOTCH1","NOTCH4")

# 挑选数据
aa <- heat_df |> 
  select(gene_symbol,PatientID,Function) |> 
  filter(gene_symbol %in% genes)

# 变成因子方便排序
aa$gene_symbol <- factor(aa$gene_symbol,levels = genes)

# 类型修改
aa[aa=="splice-3"] <- "splicing"
aa[aa=="splice-5"] <- "splicing"
```

画热图部分：


```r
p1 <- ggplot(aa, aes(factor(PatientID),fct_rev(gene_symbol)))+
  geom_tile(aes(fill=Function))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()
        )
p1
```

![热图部分](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185635959.png)

画条形图部分：


```r
up.df <- aa |> count(PatientID,Function)

p2 <- ggplot(up.df, aes(factor(PatientID),n))+
  geom_bar(stat = "identity", aes(fill=Function))+
  scale_y_continuous(name = NULL,expand = c(0,0))+
  scale_x_discrete(name=NULL)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none"
        )
p2
```

![条形图部分](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185656325.png)

拼图：


```r
library(aplot)

p3 <- p2 |> insert_bottom(p1,height = 5)
p3
```

![拼图](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185717409.png)

由于纵坐标变成了因子，没有突变的样本会直接移除，所以和原图还是有差别的。


## fig2c

可以用`trackviewer`画，但是我还不会。。只能偷个懒用`maftools`画。

这个附件也不像常见的几个软件的结果，应该是作者修改过的，所以也没办法直接用技能树的方法转换为maf。


```r
s2 <- data.table::fread("./s240164_2021_200_MOESM2_ESM.csv",skip = 2,header = T)
```


```r
names(s2)
```

```
##  [1] "SampleID"    "PatientID"   "Gene Symbol" "cHGVS"       "pHGVS"      
##  [6] "Function"    "Transcript"  "ExIn_ID"     "Cosmic ID"   "Vary Type"  
## [11] "caseAF"
```

改名字，变成`maftools`需要的列名：


```r
names(s2)[1] <- "Tumor_Sample_Barcode"
names(s2)[3] <- "Hugo_Symbol"
names(s2)[6] <- "Variant_Classification"
names(s2)[10] <- "Variant_Type"
```

增加几列`maftools`需要的列名，内容随便填即可：


```r
s2$Chromosome <- c(rep(paste0("chr",1:22),38),paste0("chr",1:20))
s2$Start_Position <- 1
s2$End_Position <- 2
s2$Reference_Allele <- 3
s2$Tumor_Seq_Allele2 <- 4
```

还要修改`Variant_Classification`的内容，不然`maftools`会报错。


```r
table(s2$Variant_Classification)
```

```
## 
##       cds-del       cds-ins    frameshift      missense         ncRNA 
##            24             5            43           712             2 
##      nonsense          span      splice-3      splice-5     stop-loss 
##            45             4             6            13             1 
## stop-retained 
##             1
```


这几个类型转换可能有问题，网络上没找到合适的信息.


```r
s3 <- s2 |> 
  mutate(Variant_Type = ifelse(Variant_Type == "SNV","SNP","DEL"),
         Variant_Classification = 
           case_when(Variant_Classification == "cds-del" ~ "In_Frame_Del",
                     Variant_Classification == "cds-ins" ~ "In_Frame_Ins",
                     Variant_Classification == "frameshift" ~ "Frame_Shift_Ins",
                     Variant_Classification == "missense" ~ "Missense_Mutation",
                     Variant_Classification == "ncRNA" ~ "RNA",
                     Variant_Classification == "nonsense" ~ "Nonsense_Mutation",
                     Variant_Classification == "span" ~ "Intron",
                     Variant_Classification %in% c("splice-3","splice-5") ~ "Splice_Site",
                     Variant_Classification %in% c("stop-loss","stop-retained") ~ "3'Flank"
                     )
         )
```




```r
library(maftools)
```

```r
ptclMaf <- read.maf(s3)
```

```
## -Validating
## --Removed 108 duplicated variants
## -Silent variants: 8 
## -Summarizing
## -Processing clinical data
## --Missing clinical data
## -Finished in 0.050s elapsed (0.030s cpu)
```

由于类型转换问题，比例差异很大。


```r
oncoplot(ptclMaf,top = 30)
```

![oncoprint](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185824186.png)

```r
lollipopPlot(ptclMaf, gene = "TET2",
             AACol = "pHGVS"
             
             )
```

```
## Removed 4 mutations for which AA position was not available
```

![image-20220903185900483](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185900483.png)


```r
lollipopPlot(ptclMaf, gene = "TP53", AACol = "pHGVS")
```

```
## 8 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.
```

```
##    HGNC    refseq.ID   protein.ID aa.length
## 1: TP53    NM_000546    NP_000537       393
## 2: TP53 NM_001126112 NP_001119584       393
## 3: TP53 NM_001126118 NP_001119590       354
## 4: TP53 NM_001126115 NP_001119587       261
## 5: TP53 NM_001126113 NP_001119585       346
## 6: TP53 NM_001126117 NP_001119589       214
## 7: TP53 NM_001126114 NP_001119586       341
## 8: TP53 NM_001126116 NP_001119588       209
```

```
## Using longer transcript NM_000546 for now.
```

```
## Removed 1 mutations for which AA position was not available
```

![image-20220903185923147](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220903185923147.png)







