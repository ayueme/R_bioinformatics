这篇推文适合初学者看，大佬酌情阅读！

从打开网址开始教你一步一步的下载TCGA的数据，图文并茂，真的是详细的不能再详细了！

如果你看完了这篇还不会下载TCGA的数据，那不是你疯了就是我疯了！

**非常适合初学者，因为使用这个方法下载TCGA数据后，只要2行代码即可提取表达矩阵，包含count/FPKM/tpm，自带gene symbol，并且附带和表达矩阵对应的详细临床信息，无需再次下载！**

在下载TCGA数据之前，你可能需要一些背景知识，比如TCGA的33癌症简称和英文名，拷贝数变异、单核苷酸多态性、甲基化等的英文，建议自己百度下哦~

首先你要到这个网址：https://portal.gdc.cancer.gov/，进入下面这个界面，如果你打不开这个页面，那你的下载大概率也会有问题的，因为这个对网络有要求！

![image-20220912163652737](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912163652737.png)



打开这个页面后，你需要选择你想要下载的东西，这个数据库下载东西逻辑是很清晰的，比如你想要下载**TCGA的直肠癌的常规转录组的mRNA数据**，首先你要点击`Repository`，下面箭头指的两个地方，任意点一个就行，都是一样的：

![image-20220912164329640](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912164329640.png)



点完了之后会进入到这个界面:

![image-20220912165913901](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912165913901.png)

这里你只要关注左侧这一栏的东西就好了，你现在的需求是下载**TCGA的直肠癌的常规转录组的mRNA数据**，所以你需要先找到**TCGA的直肠癌**，点击**Cases**。

**重点来了！！！这里是决定你能不能用2行代码整理表达矩阵的关键！！**

有的教程会让你在**Primary Site**中找到**直肠癌**，勾选它，像下图这样，但**我建议你直接跳过这一步！！！**

<img src="https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912170857082.png" alt="image-20220912170857082" style="zoom:33%;" />

<img src="https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912171149012.png" alt="image-20220912171149012" style="zoom:33%;" />



但是**我建议你直接跳过上面这一步！！直接在`Project`里面选中`TCGA-READ`即可，不要在`Primary Site`中勾选任何东西！！**

![image-20220912214851137](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912214851137.png)

这样**Cases**这边的选项就选好了，其他东西就不用选了，比如**Disease Type/Gender/Age at Diagnosis等等**。

下一步就是到**Files**里面选择数据类型，你想要的**转录组的的mRNA数据**，所以先点击**Files**，然后在下面的**Data Category**里面选择**transcriptome profiling**，在**Data Type**里面选择**Gene Expression Quantification**：

<img src="https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912215053313.png" alt="image-20220912215053313" style="zoom:33%;" />

OK，到这里你就选择好你想要的数据了：**TCGA的直肠癌的常规转录组的mRNA数据**，其他的都不用选了，你可以看到一共177个文件！

下一步，把你的所有数据添加到购物车，也就是点击右侧**Add All Files to Cart**，点完之后你的右上角购物车会出现数字，就像下面这样：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912215417751.png)



加入购物车之后，点击**Cart**，进入下面的下载界面，这个界面内容很多，比如167个cases，177个files，747.58M，还有各种下载选项，都给你标出来了。

![image-20220912215730994](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912215730994.png)

此时，你点击**Download**，下面会出现**Manifest和Cart**的选项。

你如果选择**Cart**，会下载一个压缩包，里面就是你选择好的177个文件。点击**Cart**会下载下图中标号**1**的文件，解压后得到标号**2**的文件，把**2**继续解压，就得到标号为**3**的文件夹：![image-20220912182021322](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912182021322.png)

打开这个文件夹，里面就是你的178个文件，因为多了一个**Manifest**文件。

<img src="https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912220006828.png" style="zoom:50%;" />

如果你的网络不错，直接下载**Cart**其实是非常简单的方法，比如我这里177个文件，700多M，不到10分钟就下载完了！:smile:下载完成后你如果需要整理成表达矩阵，那你还需要点击**Metadata**，下载一个metadata文件！

![image-20220912184217970](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912184217970.png)

这个metadata是json格式的，里面记录着文件名和样本名的对应关系，整理成表达矩阵需要这个文件。



如果你选择了**Manifest**，那么会下载一个manifest文件（这个文件内容和上面通过cart方式下载得到的MANIFEST文件内容完全一样）：

![image-20220912183518323](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912183518323.png)

这个文件里面包含了你选择好的177个文件的信息，有了它，你就可以用GDC官方推荐的**gdc client**工具下载了，后面我们会说怎么用这个文件进行下载。如果你要整理成表达矩阵，那你还需要点击**Metadata**，下载一个metadata文件！

除此之外，你还可以在这个界面**下载临床信息**，点击**Clinical**，下载**TSV**格式的临床数据。

![image-20220912181338874](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912181338874.png)

其他的数据用处不大，可以不用下载。

下面说说如何用GDC官方推荐的**gdc client**工具进行下载。

首先，需要到`gdc client`的下载地址：https://gdc.cancer.gov/access-data/gdc-data-transfer-tool，下载这个软件，往下拉即可看到各个平台的版本：

![image-20220912162758228](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912162758228.png)

左边的是命令行版本，右边是图形化界面，可以用鼠标点点点的那种！

我先给大家演示下点点点的版本，适合不会写代码的人！根据你的系统，下载合适的版本，比如我是Windows，我就下载了下面这个：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912184930655.png)

然后解压它，安装它，打开它，就会出现下面这个界面：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912190252435.png)

点击那个**Select Manifest File**，会让你选择你的**manifest**文件，也就是上一步下载的那个，选好之后会出现下面的界面：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912211050696.png)



稍等片刻，如果你的网络没问题就会变成下面这个界面，然后你点击右下角的**Download**就可以下载了（下载前你可以先设置下，见下一张图），下载过程会告诉你一共多少几个，下载中几个，失败几个，停止几个，完成几个等，非常清晰明了：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912220530992.png)

下载前可以先进行一些设置：比如更改你的保存路径，每次下载的大小，自动重连的次数等等，改好之后记得点击**Save Settings**：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912211611086.png)

都设置好之后，点击**Download**即可愉快的下载了！

如果有失败的，会在下面显示，直接选中继续下载即可：

![](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912213741905.png)



下载完之后也是177个文件：

![image-20220912190839901](https://aliyun-bucket0324.oss-cn-shanghai.aliyuncs.com/img/image-20220912190839901.png)

这个方法本质上和命令行版本的gdc client没有区别！

这个点点点的图形化界面真的太香了，这不比命令行版本的gdc client香多了么！还要命令行那么复杂的东西干啥？？

**根据这个教程下载后，可以无缝衔接另一篇教程：只要2行代码即可整理成表达矩阵！**

[3.手动下载的TCGA数据也是可以用TCGAbiolinks包整理的](https://mp.weixin.qq.com/s/DHj9wp6hkae2Zrl61sU1fQ)

TCGA傻瓜版下载教程未完待续，下一次说说怎么用gdc client的命令行进行下载，后面还会介绍如何整理成表达矩阵！

万里长城第一步，这才开始！