`TCGA`的数据由于正常样本较少，满足不了大家的需求，GTEx刚好有很多正常组织的测序数据，二者一起用非常合适。

但是由于两个来源的数据上游分析流程不一样，不能直接使用。那问题就来了，怎么一起用这两个数据呢？

如果你是一个热爱学习的人，肯定不会有这样的疑问，去`pubmed`一搜便知，有无数人已经用这两个数据发了不知道多少篇SCI了，你随便找一篇看看方法部分，就可以模仿了！我下面随便找了几篇，方法五花八门，有的甚至根本没说怎么处理的...各种都有：

![doi: 10.7717/peerj.8961](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114174703634.png)

![doi: 10.3389/fimmu.2021.688215](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114174918909.png)

![doi: 10.3389/fonc.2020.605097](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114175456712.png)



已发表的文章用法各有不同，对错我不去评判，这里**我给大家推荐两个比较靠谱的数据来源。**

- 一个是来自Xena
- 另一个是来自于scientific data上的一篇文章

这两种都是下载了原始数据，从头处理得到的数据，并且矫正了批次效应，也都有详细的方法说明，并且也有相关的参考文献，因此不用担心审稿人的质疑，参考文献甩他脸上...

今天先给大家介绍下scientific data文章里面的做法及数据。

文献标题：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114191135925.png)



这篇文章详细介绍了他们的数据处理流程：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114191437322.png)

原始数据下载、比对工具、流程参数等细节都写的很清楚，大家可以去阅读原文，我这里之放了一部分：

![image-20221114191550482](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114191550482.png)

并且他们的数据处理pipeline的**代码也放在了github：https://github.com/mskcc/RNAseqDB**

除此之外，**作者把处理好的数据也放在了figshare中，需要的直接去下载即可**，我已经都下载好了，肯定是没有问题的。唯一不足是他们**只处理了编码基因的数据，非编码rna和其他组学数据都没有。**

**数据处理为3种，分别放在3个仓库中，根据癌症类型都给你直接整理成表达矩阵了，并且自带gene symbol和entrez id：**

- *Data record 1*: 下载地址：https://figshare.com/articles/dataset/Data_record_1/5330539     The maximum likelihood gene expression levels computed using RSEM, i.e., the expected_count in RSEM’s output, are in Data Citation 1. This dataset includes 52 data files, each being a sample-gene matrix of a certain tissue type (see Table 1 for the tissues we processed). This dataset can be provided to programs such as edgeR for identifying differentially expressed genes.
- *Data record 2*: 下载地址：https://figshare.com/articles/dataset/Data_record_2/5330575    The gene expression levels calculated from the FPKM (Fragments Per Kilobase of transcript per Million) in RSEM’s output are in Data Citation 2. This dataset (of data files) was quantile normalized, but not corrected for batch effects.
- *Data record 3*: 下载地址：https://figshare.com/articles/dataset/Data_record_3/5330593    The normalized gene expression levels (FPKM) are in Data Citation 3. This dataset was not only quantile normalized, but was corrected for batch effects (using ComBat).

以上3种，大家选择自己喜欢的就好，别纠结！

我把expect_count的数据给大家看看：

分为3种不同的表达矩阵，稍微有点TCGA数据库的常识都能看出来，3个表达矩阵分别是TCGA的tumor、TCGA的normal(癌旁)、GTEx的表达矩阵。

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114195018665.png)

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114195041420.png)

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114195122962.png)



有了这个数据，大家就可以非常方便的去做一些下游分析，不用再担心正常样本不够的问题！

我这里演示的是**RSEM得到的结果，这里的expected_count不是整数**，你可以用三大差异分析R包中的任何一个进行差异分析。如果是用`DEseq2`，那可以直接四舍五入把小数去掉，或者可以用`tximport`包处理一下再使用`DESeqDataSetFromTximport`函数进行后续的分析。









一共需要下载3个文件：表达矩阵，样本信息文件，基因信息文件



首先到这个网址下载GTEx的数据：[UCSC Xena (xenabrowser.net)](https://xenabrowser.net/datapages/?cohort=GTEX&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

![image-20221114121526050](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20221114121526050.png)

下载id转换需要的文件：https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap



整合流程：[GitHub - mskcc/RNAseqDB](https://github.com/mskcc/RNAseqDB)

xena的建议和处理：[How do I compare tumor vs normal expression? - User Help Pages (gitbook.io)](https://ucsc-xena.gitbook.io/project/how-do-i/tumor-vs-normal)

xena处理好的数据：[UCSC Xena (xenabrowser.net)](https://xenabrowser.net/datapages/?cohort=TCGA TARGET GTEx)

nature的文章：[Unifying cancer and normal RNA sequencing data from different sources | Scientific Data (nature.com)](https://www.nature.com/articles/sdata201861)