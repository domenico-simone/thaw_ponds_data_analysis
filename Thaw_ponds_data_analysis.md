# Thaw ponds data analysis

## Quality trimming

Sickle with default parameters.

```bash
module load sickle

for i in A-1_S21_L005 A-2_S22_L005 A-3_S23_L005 A-4_S24_L005 B-1_S25_L005 B-2_S26_L005 B-4_S27_L006 C-2_S28_L006 C-3_S29_L006 C-4_S30_L006 C-5_S31_L006 D-1_S32_L006 D-2_S33_L007 D-3_S34_L007 D-4_S35_L007 E-1_S36_L007 E-2_S37_L007 E-3_S38_L007 E-4_S39_L007 F-3_S40_L008 F-4_S41_L008 F-5_S42_L008 G-1_S43_L008 I-1_S44_L008 I-2_S45_L008 I-3_S46_L008; do
sickle pe -f ${i}_R1_001.fastq.gz \
-r ${i}_R2_001.fastq.gz \
-t sanger \
-o trimmed_${i}_R1_001.fastq \
-p trimmed_${i}_R2_001.fastq \
-s trimmed_singles_${i}.fastq
done
```

## Co-assembly of all samples

```bash
module load megahit

megahit -o megahit_assembly --out-prefix thawponds_assembly -1 trimmed_A-1_S21_L005_R1_001.fastq -2 trimmed_A-1_S21_L005_R2_001.fastq -1 trimmed_A-2_S22_L005_R1_001.fastq -2 trimmed_A-2_S22_L005_R2_001.fastq -1 trimmed_A-3_S23_L005_R1_001.fastq -2 trimmed_A-3_S23_L005_R2_001.fastq -1 trimmed_A-4_S24_L005_R1_001.fastq -2 trimmed_A-4_S24_L005_R2_001.fastq -1 trimmed_B-1_S25_L005_R1_001.fastq -2 trimmed_B-1_S25_L005_R2_001.fastq -1 trimmed_B-2_S26_L005_R1_001.fastq -2 trimmed_B-2_S26_L005_R2_001.fastq -1 trimmed_B-4_S27_L006_R1_001.fastq -2 trimmed_B-4_S27_L006_R2_001.fastq -1 trimmed_C-2_S28_L006_R1_001.fastq -2 trimmed_C-2_S28_L006_R2_001.fastq -1 trimmed_C-3_S29_L006_R1_001.fastq -2 trimmed_C-3_S29_L006_R2_001.fastq -1 trimmed_C-4_S30_L006_R1_001.fastq -2 trimmed_C-4_S30_L006_R2_001.fastq -1 trimmed_C-5_S31_L006_R1_001.fastq -2 trimmed_C-5_S31_L006_R2_001.fastq -1 trimmed_D-1_S32_L006_R1_001.fastq -2 trimmed_D-1_S32_L006_R2_001.fastq -1 trimmed_D-2_S33_L007_R1_001.fastq -2 trimmed_D-2_S33_L007_R2_001.fastq -1 trimmed_D-3_S34_L007_R1_001.fastq -2 trimmed_D-3_S34_L007_R2_001.fastq -1 trimmed_D-4_S35_L007_R1_001.fastq -2 trimmed_D-4_S35_L007_R2_001.fastq -1 trimmed_E-1_S36_L007_R1_001.fastq -2 trimmed_E-1_S36_L007_R2_001.fastq -1 trimmed_E-2_S37_L007_R1_001.fastq -2 trimmed_E-2_S37_L007_R2_001.fastq -1 trimmed_E-3_S38_L007_R1_001.fastq -2 trimmed_E-3_S38_L007_R2_001.fastq -1 trimmed_E-4_S39_L007_R1_001.fastq -2 trimmed_E-4_S39_L007_R2_001.fastq -1 trimmed_F-3_S40_L008_R1_001.fastq -2 trimmed_F-3_S40_L008_R2_001.fastq -1 trimmed_F-4_S41_L008_R1_001.fastq -2 trimmed_F-4_S41_L008_R2_001.fastq -1 trimmed_F-5_S42_L008_R1_001.fastq -2 trimmed_F-5_S42_L008_R2_001.fastq -1 trimmed_G-1_S43_L008_R1_001.fastq -2 trimmed_G-1_S43_L008_R2_001.fastq -1 trimmed_I-1_S44_L008_R1_001.fastq -2 trimmed_I-1_S44_L008_R2_001.fastq -1 trimmed_I-2_S45_L008_R1_001.fastq -2 trimmed_I-2_S45_L008_R2_001.fastq -1 trimmed_I-3_S46_L008_R1_001.fastq -2 trimmed_I-3_S46_L008_R2_001.fastq -r trimmed_singles_A-1_S21_L005.fastq -r trimmed_singles_A-2_S22_L005.fastq -r trimmed_singles_A-3_S23_L005.fastq -r trimmed_singles_A-4_S24_L005.fastq -r trimmed_singles_B-1_S25_L005.fastq -r trimmed_singles_B-2_S26_L005.fastq -r trimmed_singles_B-4_S27_L006.fastq -r trimmed_singles_C-2_S28_L006.fastq -r trimmed_singles_C-3_S29_L006.fastq -r trimmed_singles_C-4_S30_L006.fastq -r trimmed_singles_C-5_S31_L006.fastq -r trimmed_singles_D-1_S32_L006.fastq -r trimmed_singles_D-2_S33_L007.fastq -r trimmed_singles_D-3_S34_L007.fastq -r trimmed_singles_D-4_S35_L007.fastq -r trimmed_singles_E-1_S36_L007.fastq -r trimmed_singles_E-2_S37_L007.fastq -r trimmed_singles_E-3_S38_L007.fastq -r trimmed_singles_E-4_S39_L007.fastq -r trimmed_singles_F-3_S40_L008.fastq -r trimmed_singles_F-4_S41_L008.fastq -r trimmed_singles_F-5_S42_L008.fastq -r trimmed_singles_G-1_S43_L008.fastq -r trimmed_singles_I-1_S44_L008.fastq -r trimmed_singles_I-2_S45_L008.fastq -r trimmed_singles_I-3_S46_L008.fastq --cpu-only -m 250000000000 --min-contig-len 1000 --k-max 101
```

## Abundance of HMM families on thaw ponds contigs

### Generate bowtie2 dbs for co-assembly contigs

```bash
module load bioinfo-tools
module load bowtie2

bowtie2-build thawponds_assembly.fa.gz thawponds_assembly
```

### Map reads to contigs

```bash
cd /home/domeni/thaw_ponds_contig_mappings/new_contig_mappings

export homeDir=`pwd`
export wdir="/home/domeni/thaw_ponds"
export dataDir="$wdir/filtered_reads"
export sampleList="$wdir/sample_list"
export outdir="$wdir/alignments"

mkdir -p logs
mkdir -p $outdir

sbatch -p node -t 20:00:00 -A snic2018-3-22 \
-J tp_contig_mappings -o logs/tp_contig_mappings_%a.out -e logs/tp_contig_mappings_%a.err \
--array=1-$(wc -l < $sampleList) \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load bowtie2/2.2.9
#module load htseq

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleList)

cp ${dataDir}/${sample}*R* ${SNIC_TMP}
cp thawponds_assembly.*.bt2 ${SNIC_TMP}

cd ${SNIC_TMP}

mv ${sample}_R1_001.qc.fastq.gz ${sample}_R1_001.qc.fastq
mv ${sample}_R2_001.qc.fastq.gz ${sample}_R2_001.qc.fastq

bowtie2 -x thawponds_assembly \
-1 ${sample}_R1_001.qc.fastq \
-2 ${sample}_R2_001.qc.fastq \
--no-unal \
-S ${sample}.sam \
--threads 20

samtools view -F 4 -bS ${sample}.sam > ${sample}.bam
samtools sort -o ${sample}.sorted.bam ${sample}.bam

cp ${sample}.sorted.bam $outdir

EOF
```
