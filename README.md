
# SCISSOR (Shape Changes In Selecting Sample Outliers in RNA-seq)

Abstract

In cancer research, a substantial proportion of human genes differ in function through shape changes in read coverage of RNA-seq. For example, tumor-suppressor genes lose their function through various changes in expression such as aberrant splicing, frameshift indel, large deletions, or overexpression of noncoding RNAs. Although previous studies have examined mutational effect on such shape changes, it has been observed that a large fraction of shape aberration occur in the absence of mutations. We have developed a statistical method to systematically detect abnormal RNA-seq samples using read coverage data alone independently of mutational analysis. To tackle the high dimensionality of read level RNA-seq data, we model the underlying mechanism possibly generating aberrant shapes by multiple unknown mixture distributions, recasting the problem as a high-dimensional latent variables framework. Based on this underlying mechanism, we normalize the read coverage to adjust for different library sizes, extract the latent information in a robust way, and determine the cases that are strongly involved in abnormality. This approach allows us to detect not only local changes in expression such as alternative splicing events and frameshift indels but also landscape changes such as fusion and a wide range of deletions. 


