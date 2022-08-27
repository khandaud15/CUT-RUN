#!/bin/bash

#SBATCH --job-name=CUT_RUN
#SBATCH --output=CUT_RUN.out
#SBATCH --error=CUT_RUN.err
#SBATCH --time=34:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=56000
#SBATCH --account=pi-piunti

module load gcc
module load python
conda activate  cut_tag

#### Note : the script is supposed to run from the directory where the fastq files are sitting
export PATH="$PATH:/home/khandaud15/.conda/envs/cut_tag"
Fasta=/scratch/midway3/khandaud15/Midway3/Genome/hg38/hg38.fa
Genome=/scratch/midway3/khandaud15/Midway3/Genome/hg38/hg38
Hg38gene=/scratch/midway3/khandaud15/Midway3/Genome/hg38/Homo_sapiens.GRCh38.93.gene.bed
Hg38chrSize=/scratch/midway3/khandaud15/Midway3/Genome/hg38/hg38.chrom.sizes
set current_dir=%cd%




### lets do the normal QC on raw data using fastqc
workdir="${PWD}" 
FastqDir=$workdir/FastQ  ### This is where your raw files are kept
FastQcDir=$workdir/Fastqc_Reports
[[ -d ${FastQcDir} ]] || mkdir -p ${FastQcDir}  


cd ${FastqDir}

for file in *.gz; do
	fastqc --threads 10 --outdir ${FastQcDir}  $file; 
done



# create bowtie2 index database (database name: hg38)
bowtie2-build $Fasta  $Genome


#### Trim the adapters using Trimmomatic
### please keep TruSeq3-PE.fa  in the directory where you ruuning the script
### To find the file  "find /home/khandaud15/ -name "TruSeq3-PE-2.fa"", This will give the apth where the file is there 
TrimmedDir=$workdir/TrimmedFiles
Adapters=$workdir/TruSeq3-PE.fa
[[ -d ${TrimmedDir} ]] || mkdir -p ${TrimmedDir}  

cd ${FastqDir}

for infile in *_R1_001.fastq.gz
do
  base=$(basename ${infile} _R1_001.fastq.gz)
   trimmomatic PE ${infile} ${base}_R2_001.fastq.gz \
                ${TrimmedDir}/${base}_1.trim.fastq.gz ${TrimmedDir}/${base}_1un.trim.fastq.gz \
                ${TrimmedDir}/${base}_2.trim.fastq.gz ${TrimmedDir}/${base}_2un.trim.fastq.gz \
                ILLUMINACLIP:${Adapters}:2:15:4:1:true MINLEN:20 -threads 10
 done


### Bowtie2 allignement 

BamDir=$workdir/AllignedBam
[[ -d ${BamDir} ]] || mkdir -p ${BamDir}  

cd ${TrimmedDir}
rm  -rf *un.trim.fastq.gz
for infile in *_1.trim.fastq.gz
do
	base=$(basename ${infile} _1.trim.fastq.gz)	
bowtie2  --dovetail  -p 16 -x $Genome -1 ${infile} -2 ${base}_2.trim.fastq.gz | samtools sort -o ${BamDir}/${base}.bam -
done



### picard mark duplicates  and remove duplicates using samtools
PicardDir=$workdir/DupBam
SamDir=$workdir/RmBam
[[ -d ${PicardDir} ]] || mkdir -p ${PicardDir}    
[[ -d ${SamDir} ]] || mkdir -p ${SamDir}  


cd ${BamDir}
for infile in *.bam
	do
	base=$(basename ${infile} .bam)	
	picard MarkDuplicates \
	I=$infile \
	O=${PicardDir}/${base}.md.bam \
	M=${PicardDir}/${base}.out.metrics.txt \
	ASSUME_SORT_ORDER=coordinate \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
	REMOVE_DUPLICATES=false  | samtools view -F 1024 -f 2 -b ${PicardDir}/${base}.md.bam -o ${SamDir}/${base}.rm.bam 

done



### alignments are converted into bedpe format using bedtools
BedDir=$workdir/BedFiles
[[ -d ${BedDir} ]] || mkdir -p ${BedDir}

cd ${SamDir}
for infile in *.rm.bam
	do
	base=$(basename ${infile} .rm.bam)
	bedtools bamtobed -i $infile -bedpe > ${BedDir}/${base}.bedpe
done

### Finally, the files pass through a final set of cleaning and sorting recommended for peak calling with SEACR 

BedGDir=$workdir/Bedgraph
[[ -d ${BedGDir} ]] || mkdir -p ${BedGDir}  


cd ${BedDir}
for infile in *.bedpe
	do
	base=$(basename ${infile} .bedpe)
  awk '$1==$4 && $6-$2 < 1000 {print $0}' ${infile} > ${BedGDir}/${base}.clean.bed
  cut -f 1,2,6 ${BedGDir}/${base}.clean.bed| sort -k1,1 -k2,2n -k3,3n > ${BedGDir}/${base}.fragments.bed
  bedtools genomecov -bg -i ${BedGDir}/${base}.fragments.bed -g $workdir/hg38.chrom.sizes> ${BedGDir}/${base}.fragments.bedgraph
done

#### lets call the peak on the bedgraph file using IgG control

IgGDir=$workdir/IgGCTRL
PeakDir=$workdir/SEACR_Peaks
[[ -d ${IgGDir} ]] || mkdir -p ${IgGDir} 
[[ -d ${PeakDir} ]] || mkdir -p ${PeakDir} 

find ${BedGDir} -name '*IgG*' -exec mv -t ${IgGDir} {} +  ### find and move to new directory

cd ${BedGDir}

for infile in *.fragments.bedgraph
	do
	base=$(basename ${infile} .fragments.bedgraph)
       bash  SEACR_1.3.sh $infile ${IgGDir}/*.fragments.bedgraph  norm stringent ${PeakDir}/${base}_seacr_control.peaks
       bash  SEACR_1.3.sh $infile  0.01 norm stringent ${PeakDir}/${base}_seacr_top0.01.peaks
done

#### create an bigwig file using bedgraphtobigwig
export PATH=$PATH:/scratch/midway3/khandaud15/sorftware/
BigwigDir=$workdir/Bigwigs
[[ -d ${BigwigDir} ]] || mkdir -p ${BigwigDir}  


cd ${BedGDir}
for infile in *.fragments.bedgraph
	do
	base=$(basename ${infile} .fragments.bedgraph)
	bedGraphToBigWig $infile $Hg38chrSize  ${BigwigDir}/${base}.bw
done


#### Heatmap over Transcription Unit using deeptools
HeatMapTU=$workdir/HeatMapTU
[[ -d ${HeatMapTU} ]] || mkdir -p ${HeatMapTU}  


cd ${BigwigDir}
##== linux command ==##
cores=8
computeMatrix scale-regions -S ${BigwigDir}/*.bw \
				-R $Hg38gene \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 3000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o ${HeatMapTU}/matrix_gene.mat.gz -p $cores

plotHeatmap -m ${HeatMapTU}/matrix_gene.mat.gz -out ${HeatMapTU}/Histone_gene.png --sortUsing sum


##### heatmap for peaks using midpoint
HeatMapCR=$workdir/HeatMapCR
[[ -d ${HeatMapCR} ]] || mkdir -p ${HeatMapCR}  


cd ${PeakDir} 
for infile in *_seacr_control.peaks.stringent.bed
	do
	base=$(basename ${infile} _seacr_control.peaks.stringent.bed)
        awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $infile > ${PeakDir}/${base}_seacr_control.peaks.summitRegion.bed
  	computeMatrix reference-point -S ${BigwigDir}/${base}.bw \
              -R ${PeakDir}/${base}_seacr_control.peaks.summitRegion.bed \
              --skipZeros -o ${HeatMapCR}/${base}_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

	plotHeatmap -m ${HeatMapCR}/${base}_SEACR.mat.gz -out ${HeatMapCR}/${base}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${base}"

done

echo "# "
echo "#        Congrats! The CUT & RUN  analysis is complete!"
echo "# "
echo "--------------------------------------------------------------------------"

