# RNAseq-pipeline

#!/bin/bash
#SBATCH -o /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2024/logs/%x_%j.out
#SBATCH -e /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2024/logs/%x_%j.err
#SBATCH --array=1-175
#SBATCH --account gniu
#SBATCH --ntasks=3
#SBATCH --mail-type ALL
#SBATCH --mem=10G
#SBATCH --time=0-10:00:00

module purge
module add hisat2
module add samtools
module add python-cbrg
module load fastqc/0.11.9
module load java/21.0.2
module load trimmomatic

out=/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2025

let j=1
let k=2
let m=3
let n=4

fileName=$(awk \
-v i="$SLURM_ARRAY_TASK_ID" -v col="$j" \
'NR==i {print $col}' \
/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2025/lookups/202503_samples.txt)

geno=$(awk \
-v i="$SLURM_ARRAY_TASK_ID" -v col="$k" \
'NR==i {print $col}' \
/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2024/lookups/2024.mbv.txt)

R1=$(awk \
-v i="$SLURM_ARRAY_TASK_ID" -v col="$k" \
'NR==i {print $col}' \
/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2025/lookups/202503_samples.txt)

R2=$(awk \
-v i="$SLURM_ARRAY_TASK_ID" -v col="$m" \
'NR==i {print $col}' \
/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2025/lookups/202503_samples.txt)

### QC
fastqc "$R1" "$R2" -o $out/fastqc

### Remove Adaptors
trimmomatic PE -phred33 \
    "$R1" "$R2" \
    $out/"$fileName".trimmed_R1.fastq.gz $out/"$fileName".unpaired_R1.fastq.gz \
    $out/"$fileName".trimmed_R2.fastq.gz $out/"$fileName".unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

### Alignment
#hisat2 -x /project/fairfaxlab/bfairfax/Leavers/nassisar/HISAT/GENOMEANNOT/GRCh38.dna -p 8 -1 $out/"$fileName".trimmed_R1.fastq.gz -2 $out/"$fileName".trimmed_R1.fastq.gz --rg-id $out/"$fileName" --rg SM:$out/"$fileName" -S $out/"$fileName".sam

### Convert
samtools view -bS $out/"$fileName".sam > $out/"$fileName".bam
samtools sort --threads 8 $out/"$fileName".bam > $out/"$fileName".sorted.bam

### Remove Duplication
MarkDuplicates I="$fileName".sorted.bam O="$fileName".nodup.bam M="$fileName".marked_dup_metrics.txt REMOVE_DUPLICATES=true
bamtools filter -tag NH:1 -in "$fileName".nodup.bam | bamtools filter -isProperPair true -out "$fileName".nodup_properPairs_NH.bam
BuildBamIndex I="$fileName".nodup_properPairs_NH.bam

rm "$fileName".sam

### Flagstat
samtools flagstat $out/"$fileName".sorted.bam > $out/"$fileName".sorted.flagstat.txt

### Gene Quantification
python -m HTSeq.scripts.count \
--format=bam \
--minaqual=0 \
--stranded=no \
--type=exon \
--mode=union \
"$fileName".nodup_properPairs_NH.bam \
/project/fairfaxlab/nassisar/HISAT/GENOMEANNOT/GRCh38.gtf > /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/2024.Jan/gene.counts1/"$fileName".txt

### STAR Mapping
vcfs=/project/fairfaxlab/gniu/bulk_map/vcf
path=/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.STAR.WASP/2024.Jan/star

module purge
module add STAR
module add samtools

STAR \
--genomeDir /project/fairfaxlab/gniu/ref_genome/ensembl/GenomeDir/ \
--readFilesIn "$R1" "$R2" \
--twopassMode Basic \
--outSAMstrandField intronMotif \
--readFilesCommand zcat \
--waspOutputMode SAMtag \
--runThreadN 5 \
--varVCFfile $vcfs/"$geno".vcf.recode.vcf \
--outSAMtype BAM Unsorted \
--outFileNamePrefix $path/"$fileName".

module unload STAR
module unload samtools


### Filtering BAMs

module add picard-tools
module add bcftools/1.9
module add samtools

#define output directories

vcfs=/project/fairfaxlab/gniu/bulk_map/vcf
path=/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.STAR.WASP/2024.Jan/star
out=/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/verifybamid/2024.Jan/filter
wasp=/project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.STAR.WASP/remappedData/working/wasp

#rather than add back header just re-give fasta 

samtools view $path/"$fileName".Aligned.out.bam | grep -v 'vW:i:2' \
| grep -v 'vW:i:3' \
| grep -v 'vW:i:4' \
| grep -v 'vW:i:5' \
| grep -v 'vW:i:6' \
| grep -v 'vW:i:7' \
| samtools view -bS -T /project/fairfaxlab/gniu/ref_genome/ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa \
> $out/"$fileName".bam

samtools \
sort \
--threads 8 \
$out/"$fileName".bam \
> $out/"$fileName".sorted.bam

#index bam
BuildBamIndex I=$out/"$fileName".sorted.bam


rm -f $out/"$fileName".bam

#filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read - using filter_remapped_reads.py

#non-differentially filter duplicate reads - using rmdup_pe.py script

module rm python
module add python-cbrg

python \
$wasp/rmdup_pe.py \
$out/"$fileName".sorted.bam \
$out/"$fileName".rm_dup.bam

rm -f $out/"$fileName".sorted.bam

samtools \
sort \
--threads 8 \
$out/"$fileName".rm_dup.bam \
> $out/"$fileName".sort.bam

#index bam
BuildBamIndex I=$out/"$fileName".sort.bam

rm -f $out/"$fileName".rm_dup.bam
rm -f $path/"$fileName".Aligned.out.bam

module unload picard-tools
module unload bcftools/1.9a
module unload samtools/1.9

vcf_all=/project/fairfaxlab/gniu/bulk_map/vcf/allchr.b38.for10X.vcf

### mbv

module add qtltools/1.2

QTLtools \
mbv \
--bam $out/"$fileName".sort.bam \
--vcf /project/fairfaxlab/xiasun/genotyping_analysis/Sep24_merging_92_402_94/final_vcf/indv_vcfs/"$geno".post_qc.recode.vcf \
-out /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/map.HISAT2/temp_directories/bulk.mapping/2025/mbv/single/"$fileName"

QTLtools \
mbv \
--bam $out/"$fileName".sort.bam \
--vcf $vcf_all \
--out /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/verifybamid/2024.Jan/mbv/"$fileName"

module unload qtltools/1.2

### verifybamid

module add samtools

#Verify Bam ID QC - sample contamination

/project/fairfaxlab/jagilchr/nk.project/may19_std_mapping/working/verifyBamfiles/verifyBamID.20120620 \
--vcf $vcf_all \
--bam $out/"$fileName".sort.bam \
--maxDepth 1000 \
--precise \
--self \
--ignoreRG \
--out /project/fairfaxlab/shared/transcriptomics/RNAseq/processed/verifybamid/2024.Jan/verifybamid/"$fileName" \
--verbose
