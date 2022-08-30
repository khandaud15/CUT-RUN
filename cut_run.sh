#!/bin/bash

#SBATCH --job-name=CUT_RUN
#SBATCH --output=CUT_RUN.out
#SBATCH --error=CUT_RUN.err
#SBATCH --time=15:05:00
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
Hg38gene=/scratch/midway3/khandaud15/Midway3/Genome/hg38/hg38_gene.tsv
Hg38chrSize=/scratch/midway3/khandaud15/Midway3/Genome/hg38/hg38.chrom.sizes
GetSummits=/scratch/midway3/khandaud15/sorftware/get_summits_broadPeak.py
ChangeBdg=/scratch/midway3/khandaud15/sorftware/change.bdg.py
GetSEACRSummits=/scratch/midway3/khandaud15/sorftware/get_summits_seacr.py
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

echo "# Congrats! Fastqc  analysis is complete!"

# create bowtie2 index database (database name: hg38)
bowtie2-build $Fasta  $Genome

echo "#  Congrats! Bowtie2 Indexing completed!"
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

echo "#  Congrats! Adapter Trimming for all sample  is complete!"
### Bowtie2 allignement 

BamDir=$workdir/AllignedBam
logdir=$workdir/logdir
[[ -d ${BamDir} ]] || mkdir -p ${BamDir}  
[[ -d ${logdir} ]] || mkdir -p ${logdir}

cd ${TrimmedDir}
rm  -rf *un.trim.fastq.gz
for infile in *_1.trim.fastq.gz
do
	base=$(basename ${infile} _1.trim.fastq.gz)	
bowtie2  --dovetail --phred33   -p 16 -x $Genome -1 ${infile} -2 ${base}_2.trim.fastq.gz  2> $logdir/${base}.bowtie2 | samtools view -bS - > ${BamDir}/${base}.bam
echo "[info] Botwie2 allignemnet for ... ""$base" completed
done

echo "# Congrats! Bowtie Allignent for all sample  is complete!"



#### filtering the unmapped read from alligned bam 

cd ${BamDir}
for infile in *.bam
do
	base=$(basename ${infile} .bam)
echo "[info] Filtering unmapped fragments... ""$base".bam
samtools view -bh -f 3 -F 4 -F 8 $infile  > ${BamDir}/${base}_filtered.bam
echo "[info] Sorting BAM... ""$base"_filtered.bam
picard SortSam  INPUT=${BamDir}/${base}_filtered.bam  OUTPUT=${BamDir}/${base}_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
done

rm -rf ${BamDir}/${base}_filtered.bam


### picard mark duplicates  and remove duplicates using samtools
PicardDir=$workdir/DupBam
SamDir=$workdir/RmBam
[[ -d ${PicardDir} ]] || mkdir -p ${PicardDir}    
[[ -d ${SamDir} ]] || mkdir -p ${SamDir}  


cd ${BamDir}
for infile in *_sorted.bam
	do
	base=$(basename ${infile} _sorted.bam)
       echo "[info] Marking duplicates... ""$base"_sorted.bam
	picard MarkDuplicates \
	I=$infile \
	O=${PicardDir}/${base}.md.bam \
	VALIDATION_STRINGENCY=SILENT  \
	M=${PicardDir}/${base}.out.metrics.txt   2> $logdir/${base}.picard.out
samtools view -bh -F 1024 ${PicardDir}/${base}.md.bam  -o  ${SamDir}/${base}.rm.bam
done

conda deactivate
echo "#  Congrats! duplicat marking and removal is completed "
conda activate macs2

### picard mark duplicates  and remove duplicates using samtools
outdir=$workdir/peakcalling/macs2.narrow 
outdirbroad=$workdir/peakcalling/macs2.broad
outdirseac=$workdir/peakcalling/seacr
[[ -d ${outdir} ]] || mkdir -p ${outdir}
[[ -d ${outdirbroad} ]] || mkdir -p ${outdirbroad}
[[ -d ${outdirseac} ]] || mkdir -p ${outdirseac}

cd ${SamDir}
for infile in *.rm.bam
       do
       base=$(basename ${infile} .rm.bam)
       echo "[info] macs2 narrow peak calling on" "$base".rm.bam
	 macs3 callpeak -t $infile -g hs -f BAMPE -n ${base} --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/${base}.macs2.narrow
       echo "[info] macs2 broad peak calling on" "$base".rm.bam
       macs3 callpeak -t $infile  -g hs -f BAMPE -n ${base} --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> $logdir/${base}.macs2.broad
       echo "[info] Getting broad peak summits of"  "$base".rm.bam
       python $GetSummits $outdirbroad/"$base"_peaks.broadPeak | sort-bed - > $outdirbroad/"$base"_summits.bed
       echo "[info] SEACR stringent peak calling" "$base".rm.bam
       macs3 callpeak -t $infile  -g hs -f BAMPE -n ${base} --outdir $outdirseac -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base".seacr
       python $ChangeBdg $outdirseac/"$base"_treat_pileup.bdg > $outdirseac/"$base"_treat_integer.bdg
done

conda deactivate

####  call peaks using seacr

conda activate cut_tag

cd $outdirseac

for infile in *_treat_integer.bdg 
	do
	base=$(basename ${infile} _treat_integer.bdg)
	echo "[info] SEACR stringent peak calling" "$base"_treat_integer.bdg
bash  SEACR_1.3.sh $infile 0.01 non stringent $outdirseac/"$base"_treat
sort-bed $outdirseac/"$base"_treat.stringent.bed > $outdirseac/"$base"_treat.stringent.sort.bed
python  $GetSEACRSummits  $outdirseac/"$base"_treat.stringent.bed | sort-bed - > $outdirseac/"$base"_treat.stringent.sort.summits.bed
done


for infile in *_treat_integer.bdg
do
base=$(basename ${infile} _treat_integer.bdg) 
rm -rf "$base"_summits.bed "$base"_peaks.xls "$base"_peaks.narrowPeak "$base"_control_lambda.bdg "$base"_treat_pileup.bdg
done

#Generating the normalized signal file with BigWig format
Normbigwig=$workdir/Bigwig/Normalized
[[ -d ${Normbigwig} ]] || mkdir -p ${Normbigwig}

eGenomeSize=2913022398
cores=8

cd ${SamDir}
for infile in *.rm.bam
	do
	base=$(basename ${infile} .rm.bam)
	samtools index $infile  -@ 10
	echo "[info] normalizing " "$base".rm.bam
	bamCoverage --bam $infile -o $Normbigwig/"$base".cpm.norm.bw \
        --binSize 10 \
        --normalizeUsing CPM \
        --effectiveGenomeSize $eGenomeSize \
        --numberOfProcessors $cores 2> $logdir/"$base".gene.bw
done


#Generating the unnormalized signal file with BigWig format
export PATH=$PATH:/scratch/midway3/khandaud15/sorftware/
unNormbigwig=$workdir/Bigwig/unNormalized
[[ -d ${unNormbigwig} ]] || mkdir -p ${unNormbigwig}


cd ${outdir}
for infile in *_treat_pileup.bdg
       do
       base=$(basename ${infile} _treat_pileup.bdg)
		echo "[info] normalizing " "$base"_treat_pileup.bdg
		LC_ALL=C sort -k1,1 -k2,2n $infile > $outdir/"$base".sort.bdg
       bedGraphToBigWig $outdir/"$base".sort.bdg  $Hg38chrSize  $unNormbigwig/"$base".sorted.bw
done

#### Heatmap over Transcription Unit using deeptools
HeatMapTU=$workdir/HeatMapTU
[[ -d ${HeatMapTU} ]] || mkdir -p ${HeatMapTU}  


cd ${Normbigwig}
##== linux command ==##
cores=8
computeMatrix scale-regions -S ${Normbigwig}/*.bw \
				                        -R $Hg38gene \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o ${HeatMapTU}/matrix_gene.mat.gz -p $cores

plotHeatmap -m ${HeatMapTU}/matrix_gene.mat.gz -out ${HeatMapTU}/Histone_gene.png --sortUsing sum


##### heatmap for peaks using midpoint
HeatMapCR=$workdir/HeatMapCR
[[ -d ${HeatMapCR} ]] || mkdir -p ${HeatMapCR}  


cd ${outdirseac} 
for infile in *_treat.stringent.bed
	do
	base=$(basename ${infile} _treat.stringent.bed)
        awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $infile > ${outdirseac}/${base}_treat.stringent.peaks.summitRegion.bed
  	computeMatrix reference-point -S $Normbigwig/"$base".cpm.norm.bw \
              -R ${outdirseac}/${base}_treat.stringent.peaks.summitRegion.bed \
              --skipZeros -o ${HeatMapCR}/${base}_SEACR.cpm.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

	plotHeatmap -m ${HeatMapCR}/${base}_SEACR.cpm.mat.gz  -out ${HeatMapCR}/${base}_SEACR_heatmap_cpm.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${base}"

done

#### heatmap  using unnormalized  bigwig file

cd ${outdirseac}
for infile in *_treat.stringent.bed
        do
	base=$(basename ${infile} _treat.stringent.bed)
        awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $infile > ${outdirseac}/${base}_treat.stringent.peaks.summitRegion.bed
        computeMatrix reference-point -S ${unNormbigwig}/"$base".sorted.bw \
              -R ${outdirseac}/${base}_treat.stringent.peaks.summitRegion.bed \
       	      --skipZeros -o ${HeatMapCR}/${base}_SEACR.un.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

        plotHeatmap -m ${HeatMapCR}/${base}_SEACR.un.mat.gz  -out ${HeatMapCR}/${base}_SEACR_heatmap_un.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${base}"

done

echo "# "
echo "#        Congrats! The CUT & RUN  analysis is complete!"
echo "# "
echo "--------------------------------------------------------------------------"
