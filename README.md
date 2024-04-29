# CNV_MCD
CNV_MCD: Detection of Copy Number Variations based on Minimum Covariance Determinant using Next-generation Sequencing Data
# Required Dependencies
```
platform: linux-64
pyod=1.0.9
rpy2=3.4.5
scikit-learn=1.2.2
pysam==0.21.0
numba==0.56.4
numpy=1.23.5
cnvpytor=1.3.1
```
# Usage
### 1、Open the file main.py and preprocess.py and modify the variables bamFilePath and refPath inside;
```
if __name__ == '__main__':
    chr_num = 21;
    train_bam_path = '/realdata/1000Genome/NA19239.chrom21.SLX.maq.SRP000032.2009_11.bam'
    output_path = "/realdata/1000Genome/CNV_MCD/NA19239_result.txt"
    detect_cnv(train_bam_path, output_path, chr_num)
```
```
 for i in range(chrNum):
        refList = read_ref_file("../hg19/chr" + str(chr_num) + ".fa", refList, i)
```
### 2、Run main.py (for testing a single BAM file) or run run.py (for testing multiple BAM files).
# Real Datasets
### The real datasets can be obtained in the following way.
-Clink this link：(https://pan.baidu.com/s/1bHGFtC4_EA5Hc9bURXEsJw) 
-extraction code：o122 
