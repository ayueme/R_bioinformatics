在进行差异分析、生存分析等下游分析时，有很多粉丝朋友对到底使用哪种类型的数据非常纠结，所以我们今天比较一下counts、tpm、fpkm、vst、cpm的表达量差异，让大家对这些数据类型有有一个直观的感受。

以TCGA-CHOL为例。

首先获取counts、tpm、fpkm表达矩阵，这个过程建议使用1行代码系列，一步到位：

- [1行代码提取6种TCGA表达矩阵和临床信息](https://mp.weixin.qq.com/s/1OBGjUKnGyiALmLafYNPUQ)
- [1行代码提取6种TCGA表达矩阵2.0版](https://mp.weixin.qq.com/s/QFGCtrIeaAIichovw6OBVw)
- [1行代码提取TCGA的6种表达矩阵是有视频教程的](https://mp.weixin.qq.com/s/u6VkBcYqakZkaNXjzNTZcw)

```R
rm(list = ls())
load(file = "G:/tcga/TCGA-mRNA/TCGA-CHOL_mRNA.Rdata")

library(tidyverse)
library(SummarizedExperiment)
```

然后我们再准备下vst格式的表达矩阵：

```R
library(DESeq2)

mrna_expr_vst <- vst(as.matrix(mrna_expr_counts))
```

再准备下cpm格式的表达矩阵：

```R
library(edgeR)

mrna_expr_cpm <- cpm(mrna_expr_counts)
```

简单看下数据情况，都是19938行，44列。

```R
dim(mrna_expr_counts)
dim(mrna_expr_fpkm)
dim(mrna_expr_tpm)
dim(mrna_expr_vst)
dim(mrna_expr_cpm)

[1] 19938    44
[1] 19938    44
[1] 19938    44
[1] 19938    44
[1] 19938    44
```

然后简单画个箱线图看看表达量分布情况：

```R
opar <- par(mfrow=c(3,2))
boxplot(mrna_expr_counts)
boxplot(mrna_expr_fpkm)
boxplot(mrna_expr_tpm)
boxplot(mrna_expr_vst)
boxplot(mrna_expr_cpm)
```

![image-20230109192100488](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230109192100488.png)

结果很清晰了吧？这里面只有vst是另类，这也是为什么vst不需要再log的原因，其他4种类型的表达量都是很大且很分散的。

接下来我们再看看其他几个数据log之后的情况。

```R
opar <- par(mfrow=c(3,2))
boxplot(log2(mrna_expr_counts+1))
boxplot(log2(mrna_expr_fpkm+1))
boxplot(log2(mrna_expr_tpm+1))
boxplot(mrna_expr_vst)
boxplot(log2(mrna_expr_cpm+1))
```

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20230109192605821.png)

这样看是不是很接近了呢？

所以大家不要纠结了！对于TCGA这种转录组数据，差异分析就用counts，使用`DESeq2`包，后续的各种分析都用vst，没啥问题。你看这篇cell的文章用的就是vst后的数据：
![doi:10.1016/j.cell.2018.03.052.](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221229100729382.png)

当然log2之后的tpm也可以用于后续的各种分析，你去pubmed搜一下就知道，大把文章用的都是log2(tpm+1)这种，当然你用log后的tpm做差异分析(limma包)也是可以的(不推荐)，可以多看看文献~