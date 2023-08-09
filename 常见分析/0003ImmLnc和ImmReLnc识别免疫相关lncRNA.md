和免疫相关的各种分析感觉都挺火的......免疫相关的lncRNA也很火。

## ImmReg获取免疫相关lncRNA

我在查找免疫相关lncRNA时发现了一个专门的网站：[ImmReg](http://bio-bigdata.hrbmu.edu.cn/ImmReg/index.jsp): http://bio-bigdata.hrbmu.edu.cn/ImmReg/index.jsp

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521110220342.png)

这个网站可以下载**免疫相关的转录因子、RNA结合蛋白、miRNA、lncRNA**，下载后就可以直接用在你的文章里！

然后顺便在网站下方的参考文献里找到了这篇发表在NC上的文章：*Pan-cancer characterization of immune-related lncRNAs identifies potential oncogenic biomarkers*（DOI：https://doi.org/10.1038/s41467-020-14802-2）

里面提到了`ImmLnc`这个R包可以计算免疫相关的lncRNA，但是非常不幸运的是文章中提到的可以下载这个R包的网址已经更新了，变成了上面的那个网站，而且已经不再提供该R包的下载了，只提供了2段代码......

但是不要紧，所有数据都可以在网站上直接下载，毕竟我们只是要一个结果而已。

在网站的`download`界面可以非常方便的下载免疫相关lncRNA：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521110850489.png)

认真读过上面提到的文献，就知道免疫相关的lncRNA只要下载红框里的`sig.txt`即可！

下载完成后我们读取到R里面。


```r
# 读取文件
lnc_pathway_sig <- data.table::fread("E:/projects/ImmLnc/lncRNA_Pathway_sig.txt")
```

这个免疫相关lncRNA的寻找过程在文献里写的还蛮清楚的，主要是2个步骤：

1. 计算mrna和lncRNA的相关性
2. 基于免疫相关通路进行GSEA富集分析

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521111814774.png)

具体细节我们就不去管了，下面我们探索下这个下载的结果。


```r
names(lnc_pathway_sig)
##  [1] "Cancer"                "lncRNA Id"             "lncRNA Symbol"        
##  [4] "Immune Pathway"        "P Value"               "P Adjust"             
##  [7] "ES"                    "Score"                 "Number of Marker Gene"
## [10] "Marker Gene"

lnc_pathway_sig[1:4,1:4]
##    Cancer       lncRNA Id lncRNA Symbol                      Immune Pathway
## 1:    ACC ENSG00000082929     LINC01587                           Cytokines
## 2:    ACC ENSG00000117242      PINK1-AS Antigen Processing and Presentation
## 3:    ACC ENSG00000130600           H19               TCR signaling Pathway
## 4:    ACC ENSG00000145063    AC062028.1    Natural Killer Cell Cytotoxicity
```

非常详细，给出了每个癌种的lncRNA和免疫相关的通路信息，还给出了P值和富集分数。

根据网站和文献所说，`FDR<0.05 & |score| > 0.995`可被认为是免疫相关lncRNA。

看下这个标准下和`TCGA-COAD/TCGA-READ`相关的免疫相关lncRNA有多少：


```r
suppressMessages(library(tidyverse))

tt <- lnc_pathway_sig %>% 
  filter(Cancer %in% c("COAD","READ")) %>% 
  filter(abs(Score)>0.995, `P Adjust` < 0.05)

length(unique(tt$`lncRNA Symbol`))
## [1] 2828
```

2828个，如果你觉得多，可以调低筛选标准。

下面我们画个图，看看在每一条免疫通路中有多少lncRNA。这个免疫相关通路也是从其他网站下载的，大家可以读原文，写的很清楚。


```r
library(tidyverse)
source("tools_plot.R")#我自己常用的一些画图相关函数

plot_df <- lnc_pathway_sig %>% 
  filter(Cancer %in% c("COAD","READ")) %>% 
  filter(abs(Score)>0.995, `P Adjust` < 0.01) %>% 
  count(`Immune Pathway`) %>% 
  mutate(n=log2(n)) %>% 
  arrange(n)

ggplot(plot_df, aes(n, fct_reorder(`Immune Pathway`,n)))+
  geom_bar(stat = "identity", aes(fill=`Immune Pathway`))+
  scale_fill_manual(values = sample(col_vector, 16))+
  labs(x="log2(number of immune-related lncrnas)", y=NULL)+
  theme_classic()
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/unnamed-chunk-4-17644417376444244.png)

这个方法简单实用哈，下载就能用，NC出品，作为参考文献也不寒碜。

## ImmReLnc识别免疫相关lncRNA

顺着这个方法我继续搜索，竟然又发现了一个类似的算法：`ImmReLnc`。

文章发表在`frontiers in genetics`，题目：*ImReLnc: Identifying Immune-Related LncRNA Characteristics in Human Cancers Based on Heuristic Correlation Optimization* (DOI: doi: 10.3389/fgene.2021.792541)

作者提出了自己的免疫相关lncRNA鉴定算法：`ImmReLnc`

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521113618633.png)

和上面那篇非常相似，而且也给出了相关的代码，在github：https://github.com/meihonggao/ImReLnc

我已经下载下来跑过了，能跑通，可以得到和文章中一样的结果。

而且作者还和`ImmLnc`的结果做了比较：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521114106887.png)

这种方法得到的lncRNA数量很少哈，比如BRCA，只有几个，不太够用。

## TilSig识别免疫相关lncRNA

继续找资料，很快就发现了还有一个识别免疫相关lncRNA的方法！

文章发表在`Journal For Immunotherapy of Cancer`，题目：*Identification of tumor immune infiltration- associated lncRNAs for improving prognosis and immunotherapy response of patients with non- small cell lung cancer* （DOI:10.1136/jitc-2019-000110）

作者提到了一种自己开发的`TILSig`的方法：
![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230521114902960.png)

然而这篇文章并没有提供相关的代码和数据，大家感兴趣的可以根据文章中提到的方法进行复现试一试。

## 其他方法

比较常见的是WGCNA识别免疫相关lncRNA，放到下一次再讲。

其实还有很多简单的方法，比如直接粗暴的计算lncRNA和免疫相关mRNA的相关性，然后再做差异分析，取个交集就是了。

但是这种方法很显然不如上面介绍的其他几种方法炫酷啦。

本次涉及的3篇文献可在后台回复**免疫相关lncRNA**获取。

