前几天有小伙伴问怎么能批量进行wilcoxon检验，我立马就想到了`rstatix`包。然后才是for循环。

接下来就演示下怎么批量进行检验。使用`tidyverse`系列和`base R` 两种方法。

## 加载R包和数据

既然要优雅，就必须少不了`tidyverse`系列！

这个`rstatix`之前介绍过很多次了，你如果想要优雅的做医学统计学，它很重要！


```r
library(rstatix)
## 
## 载入程辑包：'rstatix'
## The following object is masked from 'package:stats':
## 
##     filter
library(tidyverse)
## -- Attaching packages ----------------------------- tidyverse 1.3.1 --
## v ggplot2 3.3.5     v purrr   0.3.4
## v tibble  3.1.6     v dplyr   1.0.8
## v tidyr   1.2.0     v stringr 1.4.0
## v readr   2.1.1     v forcats 0.5.1
## -- Conflicts -------------------------------- tidyverse_conflicts() --
## x dplyr::filter() masks rstatix::filter(), stats::filter()
## x dplyr::lag()    masks stats::lag()

# 加载数据
expr <- readRDS(file = "../000files/20220409.rds")
expr[1:3,1:3]
##      GSM1026687 GSM1026688 GSM1026689
## PAX8    6.69199    7.10748    7.47710
## THRA    4.28230    5.71306    5.23010
## CCL5    8.47402    6.37751    6.02629
```


优雅的第一步，把你的数据变成整洁的长数据！毕竟，优雅的数据才能配得上优雅的操作！


```r
expr_long <- expr %>% 
  rownames_to_column(var = 'genes') %>% 
  pivot_longer(cols = - genes, names_to = 'samples',values_to = 'values') %>% # 变长变整洁
  mutate(groups = rep(c(rep('group1',3),rep('group2',3)),9025)) # 加组别

head(expr_long)
## # A tibble: 6 x 4
##   genes samples    values groups
##   <chr> <chr>       <dbl> <chr> 
## 1 PAX8  GSM1026687   6.69 group1
## 2 PAX8  GSM1026688   7.11 group1
## 3 PAX8  GSM1026689   7.48 group1
## 4 PAX8  GSM1026693   7.69 group2
## 5 PAX8  GSM1026694   6.90 group2
## 6 PAX8  GSM1026696   7.49 group2
```

接下来就是批量进行wilcoxon检验：


```r
df <- expr_long %>% 
  group_by(genes) %>% 
  wilcox_test(values ~ groups, detailed = T)

df[1:5,1:13] # 查看结果
## # A tibble: 5 x 13
##   genes    estimate .y.    group1 group2    n1    n2 statistic     p conf.low
##   <chr>       <dbl> <chr>  <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>
## 1 A1BG      -0.277  values group1 group2     3     3         3   0.7   -0.681
## 2 A1BG-AS1  -0.0408 values group1 group2     3     3         3   0.7   -0.312
## 3 A2M-AS1   -0.887  values group1 group2     3     3         1   0.2   -1.14 
## 4 A4GALT    -0.384  values group1 group2     3     3         2   0.4   -1.04 
## 5 AAAS      -0.469  values group1 group2     3     3         1   0.2   -1.59 
## # ... with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>
```


是不是非常**优雅**？确实，非常的`tidy`，耗时大概2分钟，取决于你的电脑配置，如果你有10万个基因呢？

优雅的代价就是速度慢。当然我们也有更快的方法，可能就没有那么的优雅了！

下面简单说下不进行数据转换进行批量wilcoxon检验的思路。

以我的这个数据为例，前3列是一组，后3列是一组，这样我们也是可以进行检验的。

如果你的数据和我的不一样，也可以，这里只是提供一种思路，你可以采取各种方法解决你的问题。



```r
# 用第一行试一下
fit <- wilcox.test(as.numeric(expr[1,1:3]),as.numeric(expr[1,4:6]))
fit
## 
## 	Wilcoxon rank sum exact test
## 
## data:  as.numeric(expr[1, 1:3]) and as.numeric(expr[1, 4:6])
## W = 2, p-value = 0.4
## alternative hypothesis: true location shift is not equal to 0
```

结果很完美，就是我们想要的东西。

接下来就用base R，完成9025个基因的wilcoxon检验。


```r
# 定义一个空列表
res <- list()

for( i in 1:9025){
  
  fit <- wilcox.test(as.numeric(expr[i,1:3]),as.numeric(expr[i,4:6]))
  res[[i]] <- data.frame(pvalue = fit$p.value,stats = fit$statistic)
}

res_df <- do.call(rbind, res)
res_df$gene <- rownames(expr)
```


不错，虽然是for循环，但是速度比上面那个方法快很多！


后台回复*20220409*即可获得今日数据！配合代码复制粘贴即可运行！注意路径！
