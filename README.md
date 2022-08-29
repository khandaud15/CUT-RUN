# CUT-RUN
full pipeline to process the cut&amp;run Data on slurm

`module load gcc
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
`
