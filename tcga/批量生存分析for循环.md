批量生存分析是生信数据挖掘的基础操作，我们在之前的推文中介绍了批量logrank检验和批量单因素cox分析：[批量生存分析(logrank和单因素COX)](https://mp.weixin.qq.com/s/o-gCc_1B9SQmNFrG-I6yAQ)

TCGA的数据在进行批量生存分析前需要进行一些预处理，具体请看上面这篇推文，但是即使是这样，还是偶尔会遇到`there is only 1 group`的报错，主要有以下原因：

- 有些基因表达量太低了，
- 或者在normal和tumor样本中的表达量都一样或者很接近，
- 或者根据中位数分组后，其中一个组别只有1到2个样本

我碰到的最多就是以上几个原因，如果还有，欢迎评论区留言~

避免这个报错最简单的方法就是直接用`for`循环，可以定位到有问题的基因，并且可以跳过有问题的基因！

## 准备数据

首先用我的1行代码整理表达矩阵和临床信息，会自动得到6个表达矩阵（mRNA和lncRNA的counts/fpkm/tpm）和临床信息，这样得到的表达矩阵的样本顺序和临床信息的顺序自动就是一致的，不需要再次整理。

- [1行代码提取6种TCGA表达矩阵和临床信息](https://mp.weixin.qq.com/s/1OBGjUKnGyiALmLafYNPUQ)            
- [1行代码提取6种TCGA表达矩阵2.0版](https://mp.weixin.qq.com/s/QFGCtrIeaAIichovw6OBVw)
- [1行代码提取TCGA的6种表达矩阵是有视频教程的](https://mp.weixin.qq.com/s/u6VkBcYqakZkaNXjzNTZcw)

```R
rm(list = ls())

# 加载整理好的表达矩阵和临床信息
load(file = "G:/tcga/output_expr/TCGA-COAD_clinical.rdata")
load(file = "G:/tcga/output_expr/TCGA-COAD_mrna_expr_counts.rdata")
```

根据样本名字进行分组，然后去掉没有生存信息的样本。

我们这里用的是vst后的数据，所以先用`DESeq2`包获取vst后的数据。

```R
mrna_expr_counts <- DESeq2::vst(as.matrix(mrna_expr_counts))

# 根据列名的第14,15位数判断，小于10 就是tumor
keep_samples <- as.numeric(substr(colnames(mrna_expr_counts),14,15))<10
exprset <- mrna_expr_counts[,keep_samples]

# 临床信息只要两列就够了：生存时间和生存状态
clin <- clin_info[keep_samples, c("days_to_last_follow_up","vital_status")] 
names(clin) <- c("time","event")

# 2个样本没有结局，去掉
drop <- is.na(clin$event)
clin <- clin[!drop,]
exprset <- exprset[,!drop]

# 最后把生存结局用0,1表示
clin$event <- ifelse(clin$event == "Dead",1,0)

dim(clin)
## 478   2
dim(exprset)
## 19938   478
```

这样简单的操作后还剩下19938个基因和478个样本。

>通常到这一步后还需要过滤一下低表达的基因，但是使用`for`循环的话可以省略这一步，因为可以在计算的时候跳过这样的基因，非常方便的解决`there is only 1 group`的报错！

接下来就是大家最常见的一种操作，把这两个数据框合并到一起，然后用for循环做批量生存分析。

```R
expr_clin <- cbind(clin,t(exprset))
dim(expr_clin)
expr_clin[1:4,1:4]

##                              time event   MT-CO1   MT-ND4
## TCGA-AA-A03F-01A-11R-A16W-07    0     1 19.95212 20.11403
## TCGA-G4-6314-01A-11R-1723-07 1093     0 17.86275 17.61481
## TCGA-A6-3809-01A-01R-A278-07  996     0 15.04427 12.25016
## TCGA-AZ-6605-01A-11R-1839-07   NA     1 18.52669 18.21352
```

## 批量logrank

下面就是批量进行logrank检验：

```R
library(survival)

gene <- colnames(expr_clin)[-c(1:2)]

logrank.result <- list()

for(i in 1:length(gene)){
  print(i)
  group <- ifelse(expr_clin[,gene[i]]>median(expr_clin[,gene[i]]),"high","low")
  if(length(table(group)) == 1) next
  if(length(grep("high",group)) < 3) next
  surv <- as.formula(paste('Surv(time, event)~', "group"))
  tmp <- cbind(expr_clin[,1:2],group)
  x <- survdiff(surv, data = tmp)
  pValue <- 1-pchisq(x$chisq,df=1) 
  logrank.result[[i]] <- c(gene[i],pValue)
}
## 下面会把每个基因的序号打印出来！
## 1
## 2
## 3
## 4
## 5
## 6
## 7
## 8
## 9
## 10
## 11
## 12
## 13
## 14
## 15
## 16
# 省略。。。。。
## 38
# 省略。。。。。

res.logrank <- data.frame(do.call(rbind,logrank.result))
names(res.logrank) <- c("gene","p.value")
```

这样结果就好了，可以看到我们一开始并没有进行过滤低表达基因，过程中也没有任何报错哦！

如果你非常不幸运的还是遇到了各种报错，你可以根据打印出来的序号定位到具体的基因，单独把有问题的基因拿出来，看看到底是哪里的问题。

简单看下结果：

```R
library(dplyr)

res.logrank %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) %>% 
  head()

##    gene  p.value
## 	TAOK2	0.000103602075000153		
## 	LMAN2L	0.000105203887372118		
## 	C2orf50	0.00010953355328791		
##	FOXD4	0.000110947806679973		
##	INTS3	0.000111323387814499		
##	LBX2	0.00011999292506526
```

## 批量cox

批量cox回归也可以使用同样的思路进行分析。

```{r,echo=FALSE}
library(survival)

gene <- colnames(expr_clin)[-c(1:2)]

cox.result <- list()

for(i in 1:length(gene)){
  print(i)
  group <- ifelse(expr_clin[,gene[i]]>median(expr_clin[,gene[i]]),"high","low")
 
  if(length(table(group)) == 1) next
  if(length(grep("high",group)) < 3) next
  
  surv <- as.formula(paste('Surv(time, event)~', "group"))
  tmp <- cbind(expr_clin[,1:2],group)
  x <- coxph(surv, data = tmp)
  tmp1 <- broom::tidy(x,exponentiate = T, conf.int = T)
  cox.result[[i]] <- c(gene[i],tmp1)
}
res.cox <- data.frame(do.call(rbind,cox.result))
names(res.cox)[1] <- "gene"
```

随便找几个看看结果：

```{r}
head(res.cox)
# 省略
```

关于每一列的意义，我在之前的推文中说过超多次了，大家自己翻一翻：

- [R语言生存分析之Cox回归](https://mp.weixin.qq.com/s/eeQ8PJQZunYzQmssGTimNA)

```{r}
table(res.cox$p.value<0.05)
## FALSE  TRUE 
## 18003  1387 
```

批量cox分析也做好了，是不是很简单呢？

取交集的操作就太简单了，这里就不演示了。

后期会把这些函数写进1个R包里，实现1行代码即可完成两种批量生存分析！