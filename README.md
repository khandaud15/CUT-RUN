### CUT-RUN Analysis
Full pipeline to process the CUT&RUN Dataset on the uchicago Midway3/Midway3 cluster

### Install all the tool required for running the workflow by creating anaconda environment 

```
module load gcc
module load python
conda create -n cut_run
conda activate cut_run
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda bowtie2
conda install -c bioconda samtools
conda install -c bioconda picard
conda install -c bioconda seacr
conda install -c bioconda bedtools
pip install deeptools
```
