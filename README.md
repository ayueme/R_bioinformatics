# R_bioinformatics
> 一个临床医生的生信学习之路。

这里是我在公众号**医学和生信笔记**中发表的和**生信数据挖掘**有关的推文合集，并且进行了简单的分类。

生信分析需要的数据量巨大，且多数都是可以直接下载的，所以这里并没有数据，个别的外部数据我都在推文中注明了获取方式。

## 上游分析

1. 常规转录组上游分析流程
2. 参考基因组注释文件下载
3. GISTIC2.0分析拷贝数变异数据

## 基础知识

- [生信初学者基础知识资源推荐]()
- [学生信必知的镜像设置（加快你的下载速度）]()

## 常见下游分析

### 生存分析

**批量生存分析已通过easyTCGA实现，1行代码即可完成2种批量生存分析，支持最佳截点**，详情请参考：[GitHub - ayueme/easyTCGA: Speed up your TCGA analysis!](https://github.com/ayueme/easyTCGA)

- [count、tpm、fpkm等表达量差异](https://mp.weixin.qq.com/s/aff-AX9aA2tSDa2zbB8ZRQ)
- [批量生存分析(logrank和单因素COX)](https://mp.weixin.qq.com/s/o-gCc_1B9SQmNFrG-I6yAQ)
- [解决TCGA批量生存分析报错：there is only 1 group](https://mp.weixin.qq.com/s/8AYsnJ9KiEigVPKI5NunCg)
- 在医学统计合集中还有其他关于生存分析的内容

### 差异分析

**差异分析已通过easyTCGA实现，1行代码即可完成3种差异分析分析**，详情请参考：[GitHub - ayueme/easyTCGA: Speed up your TCGA analysis!](https://github.com/ayueme/easyTCGA)

- [count、tpm、fpkm等表达量差异](https://mp.weixin.qq.com/s/aff-AX9aA2tSDa2zbB8ZRQ)
- [DESeq2差异分析及VST变换的探索](https://mp.weixin.qq.com/s/CBznByKNGwPEKIKM5U0Oyw)
- [limma差异分析，谁和谁比很重要吗？](https://mp.weixin.qq.com/s/vdkDcBzuoqCASts61efjBw)
- [单基因差异分析并绘制火山图和热图](https://mp.weixin.qq.com/s/IZ7_46zJjXba7I73-Im9hw)

### 富集分析

- [富集分析常见类型](https://mp.weixin.qq.com/s/RtF7DPXYaObiDauIQTnkFg)
- [enrichplot富集分析可视化](https://mp.weixin.qq.com/s/1mpoaZqdgymhSsMGFrCP_A)
- [GSEA富集分析可视化](https://mp.weixin.qq.com/s/cusiasAAVPBq3uIHP0EKZw)
- [Goplot富集分析可视化](https://mp.weixin.qq.com/s/DckdtQcPv48DDLyA6oZQew)
- [GseaVis富集分析可视化](https://mp.weixin.qq.com/s/hdGkcemBdRuayA2ySMH3hw)
- [simplifyEnrichment的使用示例](https://mp.weixin.qq.com/s/BmROSJCTEzHRj9yiM8rcmA)
- [单基因富集分析](https://mp.weixin.qq.com/s/q6nkujgTYlbOQpkENjyyxA)
- [GSVA和ssGSEA](https://mp.weixin.qq.com/s/aUEP6XnejtHohaPeeEOMOQ)

### 免疫浸润

- [1行代码完成8种免疫浸润分析](https://mp.weixin.qq.com/s/JqO7rVBMGGmOXRA8w8nDSg)
- [免疫浸润结果可视化](https://mp.weixin.qq.com/s/YcUVElp0BEj5TxEqfSEkIQ)
- [免疫浸润结果分子分型（一致性聚类ConsensusClusterPlus）](https://mp.weixin.qq.com/s/96s_hfBH0HjLvvTfNgTIlQ)

### 批次效应

- [批次效应去除之combat和removebatcheffect](https://mp.weixin.qq.com/s/yRUmVTimI9f9itoHWxyYrA)

### WGCNA

- [WGCNA实战：识别免疫相关lncRNA](https://mp.weixin.qq.com/s/Pr33WscVtNQQaoryxTiJ-Q)

### 分子分型

- [免疫浸润结果分子分型（一致性聚类ConsensusClusterPlus）](https://mp.weixin.qq.com/s/96s_hfBH0HjLvvTfNgTIlQ)

## TCGA

- [新版TCGA数据库学习：批量下载数据](https://mp.weixin.qq.com/s/m8w1L4N2aXAIers_ZJvp_g)；注意：该方法已被easyTCGA替代
- [新版TCGA数据库学习：表达矩阵提取（mRNA/lncRNA/counts/tpm/fpkm）](https://mp.weixin.qq.com/s/wI0_GyVl5LiKAjX5C3f-NQ)；注意：该方法已被easyTCGA替代
- [手动下载的TCGA数据也是可以用TCGAbiolinks包整理的](https://mp.weixin.qq.com/s/DHj9wp6hkae2Zrl61sU1fQ)；注意：该方法已被easyTCGA替代
- [可能是最适合初学者的TCGA官网下载和表达矩阵整理教程](https://mp.weixin.qq.com/s/rbnWvstRsfhbi9il-qSYpQ)；注意：该方法可以配合easyTCGA一起用
- [TCGA官网下载的文件数量竟然和TCGAbiolinks不一致！](https://mp.weixin.qq.com/s/EuEaaBFjK6n6rxmUo27UMw)
- [1行代码提取6种TCGA表达矩阵和临床信息](https://mp.weixin.qq.com/s/1OBGjUKnGyiALmLafYNPUQ)；注意：该方法已被easyTCGA替代
- [1行代码提取6种TCGA表达矩阵2.0版](https://mp.weixin.qq.com/s/QFGCtrIeaAIichovw6OBVw)；注意：该方法已被easyTCGA替代
- [1行代码提取TCGA的6种表达矩阵是有视频教程的](https://mp.weixin.qq.com/s/u6VkBcYqakZkaNXjzNTZcw)；注意：该方法已被easyTCGA替代，b站有教程
- [新版TCGA数据库不同癌种的组学数据合并](https://mp.weixin.qq.com/s/0hcQ1m_9l1TtvXgEG20F5Q)；注意：该方法已被easyTCGA替代
- [TCGA不同癌种数据合并问题](https://mp.weixin.qq.com/s/2umNgpBSI19gqeqcamOv4A)；注意：该方法已被easyTCGA替代
- [TCGA的maf突变文件不能下载了？直接用TCGAbiolinks包搞定！](https://mp.weixin.qq.com/s/GpXovlWS_MAKdoRv3OAjCw)；注意：该方法已被easyTCGA替代
- [maftools需要的文件如何自己整理](https://mp.weixin.qq.com/s/1cR3Cnfd5Co9U3jIoIWJBA)
- [新版TCGAbiolinks包学习：差异分析](https://mp.weixin.qq.com/s/0SLQOZRkZ4hOQY1ETnQRUA)；注意：该方法已被easyTCGA替代
- [新版TCGAbiolinks包学习：富集分析和生存分析](https://mp.weixin.qq.com/s/z4Pl7D8tA24bHJL6eyTMlw)；注意：生存分析已被easyTCGA替代
- [新版TCGAbiolinks包学习：可视化](https://mp.weixin.qq.com/s/j0f1MDwlNmViqUeXU_Ikow)
- [1行代码计算肿瘤突变负荷TMB](https://mp.weixin.qq.com/s/TPURe613FXKi1tMHzAcJFA)
- [解决TCGA批量生存分析报错：there is only 1 group](https://mp.weixin.qq.com/s/8AYsnJ9KiEigVPKI5NunCg)；注意：该方法已被easyTCGA替代

### easyTCGA

详情请参考：[GitHub - ayueme/easyTCGA: Speed up your TCGA analysis!](https://github.com/ayueme/easyTCGA)

- [easyTCGA：让初学者也能享受“征服”TCGA的喜悦](https://mp.weixin.qq.com/s/kvGYVCOSBgKqVaeQU01JcA)
- [easyTCGA：1行代码搞定TCGA的6种表达矩阵和临床信息](https://mp.weixin.qq.com/s/z1fgyXLZXwmoaI39f2ftYw)
- [easyTCGA：1行代码搞定TCGA突变maf文件下载和整理](https://mp.weixin.qq.com/s/GBkB8Hv45l06BVnyFNFzzw)
- [easyTCGA生存分析支持最佳截点，任意基因在不同组中的表达量箱线图](https://mp.weixin.qq.com/s/Qc9m6hX-qKVJt5GzrXY9bA)

- [TCGA、GTEx的泛癌数据也是1行代码整理](https://mp.weixin.qq.com/s/SzGB1wVH_DNBbXxvkBe5NA)

- [任意基因在泛癌中的表达量展示](https://mp.weixin.qq.com/s/MIDRG57oRSMTyX6Gm99-3w)

- [b站视频教程](https://space.bilibili.com/42460432)；视频教程略旧，很多新特性没介绍

## 甲基化

- [TCGAbiolinks的甲基化数据分析](https://mp.weixin.qq.com/s/xbgQvGr0Q5DzBUqg8b__Zg)
- [ChAMP分析甲基化数据：样本信息csv的制作和IDAT读取](https://mp.weixin.qq.com/s/O_W-P_HpziXtNMZXZm8b4w)
- [ChAMP分析甲基化数据：标准流程](https://mp.weixin.qq.com/s/1xpT1E4BaWG-ulrCzylwrA)
- [ChAMP分析甲基化数据：从β值矩阵开始的流程](https://mp.weixin.qq.com/s/5x4oeJ6E0BPqtTjmEFPMcg)
- [ChAMP分析TCGA结直肠癌的甲基化数据！](https://mp.weixin.qq.com/s/TB3LTaq55yqL-Z95wY-rQA)
- [minfi包处理甲基化数据](https://mp.weixin.qq.com/s/E8j6KhEigcALgXA8fZIs9Q)

## 泛癌

- [TCGA、GTEx的泛癌数据也是1行代码整理](https://mp.weixin.qq.com/s/SzGB1wVH_DNBbXxvkBe5NA)

- [任意基因在泛癌中的表达量展示](https://mp.weixin.qq.com/s/MIDRG57oRSMTyX6Gm99-3w)

## 单细胞

- [单细胞入门之Seurat标准流程](https://mp.weixin.qq.com/s/ymdhvgcqyek2wGsDgKChfg)
- [单细胞入门之多样本整合](https://mp.weixin.qq.com/s/3w_-rYSdA31xxH83qaUy2Q)
- [单细胞入门之细胞类型鉴定](https://mp.weixin.qq.com/s/Sdx9oLC9LII7iyYl0VLKlg)

## SIC图表学习

### [Exp Hematol Oncol Fig1 Fig2a/b](https://ehoonline.biomedcentral.com/articles/10.1186/s40164-021-00200-x)

- [ggplot2绘制突变全景图](https://mp.weixin.qq.com/s/IOk1Lbi3sVIDjwMk5Jz-iA)

### [Nat.Commun Fig1](https://www.nature.com/articles/s41467-022-28421-6)

- [批次效应去除之combat和removebatcheffect](https://mp.weixin.qq.com/s/yRUmVTimI9f9itoHWxyYrA)
- [免疫浸润结果分子分型（一致性聚类ConsensusClusterPlus）](https://mp.weixin.qq.com/s/96s_hfBH0HjLvvTfNgTIlQ)
- [免疫相关lncRNA的识别](https://mp.weixin.qq.com/s/jrgZ6brGyrh1cAnW6Ddp3w)
- [WGCNA实战：识别免疫相关lncRNA](https://mp.weixin.qq.com/s/Pr33WscVtNQQaoryxTiJ-Q)

### [Molecular Cancer Fig1E Fig2A](https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-020-01170-0)

- [3D-PCA图](https://mp.weixin.qq.com/s/LTQIWYW86QCOEu7fctF8xQ)
- [R语言生信图表学习之网络图](https://mp.weixin.qq.com/s/t8UrYMO5fDkFjB2GI8WuXQ)

## 其他

- [免疫相关lncRNA的识别](https://mp.weixin.qq.com/s/jrgZ6brGyrh1cAnW6Ddp3w)

