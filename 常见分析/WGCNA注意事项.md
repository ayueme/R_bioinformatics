> We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples. In a typical high-throughput setting, correlations on fewer than 15 samples will simply be too noisy for the network to be biologically meaningful. If at all possible, one should have at least 20 samples; as with any analysis methods, more samples usually lead to more robust and refined results.



> **We do not recommend filtering genes by differential expression.** WGCNA is designed to be an unsupervised analysis method that clusters genes based on their expression profiles. Filtering genes by differential expression will lead to a set of correlated genes that will essentially form a single (or a few highly correlated) modules. It also completely invalidates the scale-free topology assumption, so choosing soft thresholding power by scale-free topology fit will fail.



> **Signed networks.** The choice of signed vs. unsigned networks is complex, but in general we prefer signed (or "signed hybrid") networks to unsigned networks. 



>**Can WGCNA be used to analyze RNA-Seq data?**



Yes. As far as WGCNA is concerned, working with (properly normalized) RNA-seq data isn't really any different from working with (properly normalized) microarray data.

We suggest removing features whose counts are consistently low (for example, removing all features that have a count of less than say 10 in more than 90% of the samples) because such low-expressed features tend to reflect noise and correlations based on counts that are mostly zero aren't really meaningful. The actual thresholds should be based on experimental design, sequencing depth and sample counts.

We then recommend a variance-stabilizing transformation. For example, package DESeq2 implements the function `varianceStabilizingTransformation` which we have found useful, but one could also start with normalized counts (or RPKM/FPKM data) and log-transform them using `log2(x+1)`. For highly expressed features, the differences between full variance stabilization and a simple log transformation are small.

Whether one uses RPKM, FPKM, or simply normalized counts doesn't make a whole lot of difference for WGCNA analysis as long as all samples were processed **the same way**. These normalization methods make a big difference if one wants to compare expression of gene A to expression of gene B; but WGCNA calculates correlations for which gene-wise scaling factors make no difference. (Sample-wise scaling factors of course do, so samples do need to be normalized.)

If data come from different batches, we recommend to check for batch effects and, if needed, adjust for them. We use ComBat for batch effect removal but other methods should also work.

Finally, we usually check quantile scatterplots to make sure there are no systematic shifts between samples; if sample quantiles show correlations (which they usually do), quantile normalization can be used to remove this effect.









- If the reader has access to a large workstation with more than 4 GB of memory, the parameter maxBlockSize can be increased. A 16GB workstation should handle up to 20000 probes; a 32GB workstation should handle perhaps 30000. A 4GB standard desktop or a laptop may handle up to 8000-10000 probes, depending on operating system and ihow much memory is in use by other running programs. 

- If a computer with large-enough memory is not available, the reader should follow Section 2.c, Dealing with large datasets, and adapt the code presented there for their needs. In general it is preferable to analyze a data set in one block if possible, although in Section 2.c we present a comparison of block-wise and single-block analysis that indicates that the results are very similar.





一些常见的注意事项：

- 15个样本以上
- counts、tpm、fpkm，芯片数据都可以，count建议使用`DESeq2`进行表转化，其他建议进行`log2`转换
- 选signed