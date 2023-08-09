生信技能树之前推过一篇推文：多次差异分析难道就需要多个火山图吗？

里面介绍过一种使用条形图展示多组差异分析结果。

今天学习下这个图的画法，不过不是多组差异分析的结果，而是多个癌症的差异分析结果，就是用TCGA的33种癌症。

## 差异分析

首先读取33个癌症的raw counts表达矩阵，这个是从xena网站下载的，简单处理下即可。生息技能树也有教程的，搜索即可获得！

### 读取数据

```R
# 清空变量
rm(list = ls())
# 加载R包
library(tinyarray)

# 获取路径
fs <- list.files('Rdata/',pattern = 'htseq_counts')
fs

[1] "TCGA-ACC.htseq_counts.Rdata" 
 [2] "TCGA-BLCA.htseq_counts.Rdata"
 [3] "TCGA-BRCA.htseq_counts.Rdata"
 [4] "TCGA-CESC.htseq_counts.Rdata"
 [5] "TCGA-CHOL.htseq_counts.Rdata"
 [6] "TCGA-COAD.htseq_counts.Rdata"
 [7] "TCGA-DLBC.htseq_counts.Rdata"
 [8] "TCGA-ESCA.htseq_counts.Rdata"
 [9] "TCGA-GBM.htseq_counts.Rdata" 
[10] "TCGA-HNSC.htseq_counts.Rdata"
[11] "TCGA-KICH.htseq_counts.Rdata"
[12] "TCGA-KIRC.htseq_counts.Rdata"
[13] "TCGA-KIRP.htseq_counts.Rdata"
[14] "TCGA-LAML.htseq_counts.Rdata"
[15] "TCGA-LGG.htseq_counts.Rdata" 
[16] "TCGA-LIHC.htseq_counts.Rdata"
[17] "TCGA-LUAD.htseq_counts.Rdata"
[18] "TCGA-LUSC.htseq_counts.Rdata"
[19] "TCGA-MESO.htseq_counts.Rdata"
[20] "TCGA-OV.htseq_counts.Rdata"  
[21] "TCGA-PAAD.htseq_counts.Rdata"
[22] "TCGA-PCPG.htseq_counts.Rdata"
[23] "TCGA-PRAD.htseq_counts.Rdata"
[24] "TCGA-READ.htseq_counts.Rdata"
[25] "TCGA-SARC.htseq_counts.Rdata"
[26] "TCGA-SKCM.htseq_counts.Rdata"
[27] "TCGA-STAD.htseq_counts.Rdata"
[28] "TCGA-TGCT.htseq_counts.Rdata"
[29] "TCGA-THCA.htseq_counts.Rdata"
[30] "TCGA-THYM.htseq_counts.Rdata"
[31] "TCGA-UCEC.htseq_counts.Rdata"
[32] "TCGA-UCS.htseq_counts.Rdata" 
[33] "TCGA-UVM.htseq_counts.Rdata" 
```

```r
# 批量读取
dat_list <- lapply(fs, function(x){
  pro <- gsub('.htseq_counts.Rdata','',x)
  load(file = file.path('Rdata/',x)) 
  dat <- pd_mat
  return(dat)
})

# 获取名称
types <- gsub(".htseq_counts.Rdata","",fs)
types

[1] "TCGA-ACC"  "TCGA-BLCA" "TCGA-BRCA" "TCGA-CESC"
 [5] "TCGA-CHOL" "TCGA-COAD" "TCGA-DLBC" "TCGA-ESCA"
 [9] "TCGA-GBM"  "TCGA-HNSC" "TCGA-KICH" "TCGA-KIRC"
[13] "TCGA-KIRP" "TCGA-LAML" "TCGA-LGG"  "TCGA-LIHC"
[17] "TCGA-LUAD" "TCGA-LUSC" "TCGA-MESO" "TCGA-OV"  
[21] "TCGA-PAAD" "TCGA-PCPG" "TCGA-PRAD" "TCGA-READ"
[25] "TCGA-SARC" "TCGA-SKCM" "TCGA-STAD" "TCGA-TGCT"
[29] "TCGA-THCA" "TCGA-THYM" "TCGA-UCEC" "TCGA-UCS" 
[33] "TCGA-UVM" 
```

### 批量对33种癌症进行差异分析

使用`DESeq2`进行差异分析，真的是3个差异分析R包中最简单的一个了，直接一句代码即可！

```R
library(DESeq2)

for(i in 1:length(dat_list)){
  group <- make_tcga_group(dat_list[[i]])
  if(length(unique(group)) <2) next # 有的癌症没有正常样本，不要
    
  coldata <- data.frame(sample = colnames(dat_list[[i]]),
                        group = group)
    
  dds <- DESeqDataSetFromMatrix(
    countData = dat_list[[i]],
    colData = coldata,
    design = ~ group
  )
    
  dds <- DESeq(dds) # deseq2用起来是真的简单
  res <- results(dds,tidy = T)
  res2 <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,] # 提取差异基因
  save(res2, file = paste0(types[i],"_diff_genes.Rdata"))
}
```

最终获得了24个癌症的差异基因。

![image-20220311210128705](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220311210128705.png)

把数据都保存好，接下来就可以画图了。



## 多个差异分析结果的条形图

### 构造画图数据

```R
# 读取数据
rm(list = ls())

files <- list.files("./", pattern = ".Rdata")

alldiff <- lapply(files, function(x){
  load(file = x)
  dat <- res2
  return(dat)
})
  
names(alldiff) <- gsub("-diff.Rdata","",files) 

alldiff <- lapply(alldiff,na.omit) # 去掉NA
```

**数就是图，图就是数**，所以最重要的就是构造数据。

```R
df_dim <- lapply(alldiff, function(x){
  dime <- cut(x$log2FoldChange,
               breaks = c(-Inf,-4,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,4,Inf)) # 切割值
  return(data.frame(table(dime)))
})

library(dplyr)

df <- bind_rows(df_dim) # 列表变成数据框
df$type <- rep(names(alldiff),each=13) # 添加癌症分类

str(df)

'data.frame':	312 obs. of  3 variables:
 $ dime: Factor w/ 13 levels "(-Inf,-4]","(-4,-3]",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ Freq: int  80 240 211 311 454 742 0 960 576 401 ...
 $ type: chr  "TCGA-BLCA" "TCGA-BLCA" "TCGA-BLCA" "TCGA-BLCA" ...
```

`df`就长这样，整洁的长数据，适合画图：

![image-20220311211035909](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220311211035909.png)



### 画图

```R
library(ggplot2)


ggplot(df, aes(x=type, y=Freq, fill=dime))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Freq),position = position_stack(vjust = 0.5),size=2)
```

第一个图长这样，不太好看的样子...下面美化一下。

![image-20220311211201211](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220311211201211.png)



美化一下：

```R
df_sum <- df %>% group_by(type) %>% summarise(n=sum(Freq)) # 计算每个癌症总的差异基因数目

# 构造色盘
library(RColorBrewer)
my_col <- c(rev(brewer.pal(6,"Reds")),"white",brewer.pal(6,"Blues"))

# 画图
ggplot()+
  geom_bar(data = df,mapping = aes(x=type, y=Freq, fill=dime),stat = "identity")+
  geom_text(data = df_sum,mapping = aes(x=type,y=n,label=n),vjust=-1)+
  scale_fill_manual(values = my_col,name = NULL)+
  scale_y_continuous(expand = c(0,0),name = "Counts",limits = c(0,8000))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title.x = element_blank())
```

![tcga_batch](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/tcga_batch.png)

这样看起来还是不错的！

但是这个图有什么用呢？