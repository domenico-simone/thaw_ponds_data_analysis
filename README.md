# Carbon degradation in thaw ponds: data analysis
<!-- TOC START min:2 max:3 link:true update:true -->
- [Parameters used in this protocol (edit this!)](#parameters-used-in-this-protocol-edit-this)
- [Download data](#download-data)
- [Install MetaDomain](#install-metadomain)
- [Quality trimming](#quality-trimming)
- [Co-assembly of all samples](#co-assembly-of-all-samples)
- [Abundance of HMM families in read datasets](#abundance-of-hmm-families-in-read-datasets)
    - [Run MetaDomain](#run-metadomain)
- [Abundance of HMM families on thaw ponds contigs](#abundance-of-hmm-families-on-thaw-ponds-contigs)
    - [Find CDS with Prodigal](#find-cds-with-prodigal)
    - [Annotate CDS against PFAM](#annotate-cds-against-pfam)
    - [Re-map hmmscan results on contig annotations](#re-map-hmmscan-results-on-contig-annotations)
    - [Generate bowtie2 dbs for co-assembly of contigs](#generate-bowtie2-dbs-for-co-assembly-of-contigs)
    - [Map reads to contigs](#map-reads-to-contigs)
    - [Get read counts](#get-read-counts)

<!-- TOC END -->

Paper results have been produced with the commands in this document, with the exception of differential abundance analyses which have been done with the commands in the Rmd file thaw_ponds_DA.Rmd.
 
## Parameters used in this protocol (edit this!)

**Run before anytime you run the protocol.**

```bash
export snicProj=
export wdir=
export PATH=${wdir}/scripts:$PATH
export PATH=${wdir}/ext_utils/Metadomain:$PATH
```
## Download data

```bash
cd $wdir
mkdir -p data && cd data
wget https://export.uppmax.uu.se/uppstore2018171/Pfam-A.hmm
wget https://export.uppmax.uu.se/uppstore2018171/Pfam-A.hmm.dat.lengths
wget https://export.uppmax.uu.se/uppstore2018171/carbohydrate_degradation_pfams_tveit.csv
wget https://export.uppmax.uu.se/uppstore2018171/thaw_ponds_datasets.tar.gz
wget https://export.uppmax.uu.se/uppstore2018171/sampleList
cd ..

mkdir -p scripts && cd scripts
wget https://export.uppmax.uu.se/uppstore2018171/splitSeqFile.py && chmod u+x splitSeqFile.py
cd ..
```

## Install MetaDomain

```bash
cd $wdir
mkdir -p ext_utils && cd ext_utils

curl -L https://sourceforge.net/projects/metadomain/files/MetaDomain.tar.gz/download > metadomain.tar.gz
tar -xvzf metadomain.tar.gz && cd MetaDomain && make
```

## Quality trimming

Sickle with default parameters.

```bash
cd $wdir
export sampleList="$wdir/data/sampleList"

mkdir -p $wdir/logs
mkdir -p $wdir/reads_filtered

sbatch -p node -t 20:00:00 -A $snicProj \
-J read_q_trimming -o logs/read_q_trimming.out -e logs/read_q_trimming.err \
--array=1-$(wc -l < $sampleList)<<'EOF'
#!/bin/bash

module load sickle

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleList)

sickle pe -f ${wdir}/data/${sample}_R1_001.fastq.gz \
-r ${wdir}/data/${sample}_R2_001.fastq.gz \
-t sanger \
-o ${wdir}/reads_filtered/trimmed_${sample}_R1_001.fastq \
-p ${wdir}/reads_filtered/trimmed_${sample}_R2_001.fastq \
-s ${wdir}/reads_filtered/trimmed_singles_${sample}.fastq

EOF
```

## Co-assembly of all samples

```bash
sbatch -p node -t 20:00:00 -A $snicProj \
-J tp_contig_mappings -o logs/tp_contig_mappings_%a.out -e logs/tp_contig_mappings_%a.err \
--array=1-$(wc -l < $sampleList)<<'EOF'
#!/bin/bash

module load megahit

megahit -o megahit_assembly --out-prefix thawponds_assembly -1 $wdir/reads_filtered/trimmed_A-1_S21_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_A-1_S21_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_A-2_S22_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_A-2_S22_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_A-3_S23_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_A-3_S23_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_A-4_S24_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_A-4_S24_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_B-1_S25_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_B-1_S25_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_B-2_S26_L005_R1_001.fastq -2 $wdir/reads_filtered/trimmed_B-2_S26_L005_R2_001.fastq -1 $wdir/reads_filtered/trimmed_B-4_S27_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_B-4_S27_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_C-2_S28_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_C-2_S28_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_C-3_S29_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_C-3_S29_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_C-4_S30_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_C-4_S30_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_C-5_S31_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_C-5_S31_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_D-1_S32_L006_R1_001.fastq -2 $wdir/reads_filtered/trimmed_D-1_S32_L006_R2_001.fastq -1 $wdir/reads_filtered/trimmed_D-2_S33_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_D-2_S33_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_D-3_S34_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_D-3_S34_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_D-4_S35_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_D-4_S35_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_E-1_S36_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_E-1_S36_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_E-2_S37_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_E-2_S37_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_E-3_S38_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_E-3_S38_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_E-4_S39_L007_R1_001.fastq -2 $wdir/reads_filtered/trimmed_E-4_S39_L007_R2_001.fastq -1 $wdir/reads_filtered/trimmed_F-3_S40_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_F-3_S40_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_F-4_S41_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_F-4_S41_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_F-5_S42_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_F-5_S42_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_G-1_S43_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_G-1_S43_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_I-1_S44_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_I-1_S44_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_I-2_S45_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_I-2_S45_L008_R2_001.fastq -1 $wdir/reads_filtered/trimmed_I-3_S46_L008_R1_001.fastq -2 $wdir/reads_filtered/trimmed_I-3_S46_L008_R2_001.fastq -r $wdir/reads_filtered/trimmed_singles_A-1_S21_L005.fastq -r $wdir/reads_filtered/trimmed_singles_A-2_S22_L005.fastq -r $wdir/reads_filtered/trimmed_singles_A-3_S23_L005.fastq -r $wdir/reads_filtered/trimmed_singles_A-4_S24_L005.fastq -r $wdir/reads_filtered/trimmed_singles_B-1_S25_L005.fastq -r $wdir/reads_filtered/trimmed_singles_B-2_S26_L005.fastq -r $wdir/reads_filtered/trimmed_singles_B-4_S27_L006.fastq -r $wdir/reads_filtered/trimmed_singles_C-2_S28_L006.fastq -r $wdir/reads_filtered/trimmed_singles_C-3_S29_L006.fastq -r $wdir/reads_filtered/trimmed_singles_C-4_S30_L006.fastq -r $wdir/reads_filtered/trimmed_singles_C-5_S31_L006.fastq -r $wdir/reads_filtered/trimmed_singles_D-1_S32_L006.fastq -r $wdir/reads_filtered/trimmed_singles_D-2_S33_L007.fastq -r $wdir/reads_filtered/trimmed_singles_D-3_S34_L007.fastq -r $wdir/reads_filtered/trimmed_singles_D-4_S35_L007.fastq -r $wdir/reads_filtered/trimmed_singles_E-1_S36_L007.fastq -r $wdir/reads_filtered/trimmed_singles_E-2_S37_L007.fastq -r $wdir/reads_filtered/trimmed_singles_E-3_S38_L007.fastq -r $wdir/reads_filtered/trimmed_singles_E-4_S39_L007.fastq -r $wdir/reads_filtered/trimmed_singles_F-3_S40_L008.fastq -r $wdir/reads_filtered/trimmed_singles_F-4_S41_L008.fastq -r $wdir/reads_filtered/trimmed_singles_F-5_S42_L008.fastq -r $wdir/reads_filtered/trimmed_singles_G-1_S43_L008.fastq -r $wdir/reads_filtered/trimmed_singles_I-1_S44_L008.fastq -r $wdir/reads_filtered/trimmed_singles_I-2_S45_L008.fastq -r $wdir/reads_filtered/trimmed_singles_I-3_S46_L008.fastq --cpu-only -m 250000000000 --min-contig-len 1000 --k-max 101

EOF
```

## Abundance of HMM families in read datasets

### Run MetaDomain

To optimize running time, read datasets are divided in chunks of 1 million reads each.

```bash
cd ${wdir}
# create tmp directory
export inDir=${wdir}/filtered_reads
export dataDir=${wdir}/data
export outDir=${wdir}/metadomain_1M
export logDir=${wdir}/logs
#mkdir -p ${inDir}
mkdir -p ${outDir}
mkdir -p ${logDir}
for sample in $(cat $dataDir/sampleList); do echo ${sample} ; python -c "
import sys, os
from itertools import izip

def runMetaDomain(sample_name, chunks, outDir='./pfam_reads_chunks_MetaDomain', logDir='./logs', dataDir='./pfam_db_2018', direction=1, inDir='./filtered_reads'):
    #outDir=
    os.system('sbatch -J pfam_reads_MetaDomain_%s_%d -o %s/pfam_reads_MetaDomain_%s_%d.out -e %s/pfam_reads_MetaDomain_%s_%d.err \
        scripts/run_MetaDomain.sh %s %s %s %s %s %s' % (sample_name, chunks, logDir, sample_name, chunks, logDir, sample_name, chunks, \
            outDir, sample_name, chunks, direction, dataDir, inDir))

sample_name = sys.argv[1]
inDir = sys.argv[2]
outDir = sys.argv[3]
dataDir = sys.argv[4]

for i in (1, 2):
    print sample_name, i
    S1 = open('%s/trimmed_%s_R%d_001.fastq' % (inDir, sample_name, i), 'r')
    #S1 = open('%s/%s_R%d_001.qc.fastq.gz' % (inDir, sample_name, i), 'r')
    chunks = 0
    print 'chunk number: %d' % chunks
    outhandle1 = open('%s/%s.%d.%d.fasta' % (inDir, sample_name, chunks, i), 'w')
    c = 1
    n_seqs = 1
    for l1 in S1:
        if c == 1:
            l1 = l1.split()
            outhandle1.write('>%s.1\n' % l1[0][1:])
            #outhandle2.write('>%s.2\n' % l2[0][1:])
        elif c == 2:
            outhandle1.write(l1)
        c += 1
        if c == 4:
            c = 0
            n_seqs += 1
            if n_seqs%500000 == 0:
                print n_seqs
                outhandle1.close()
                runMetaDomain(sample_name, chunks, outDir=outDir, logDir='./logs', dataDir=dataDir, inDir=inDir, direction=i)
                chunks += 1
                print 'chunk number: %d' % chunks
                outhandle1 = open('%s/%s.%d.%d.fasta' % (inDir, sample_name, chunks, i), 'w')
" ${sample} ${inDir} ${outDir} ${dataDir}> ${logDir}/${sample}_chunk_MetaDomain_submission.log; done
```

MetaDomain results will be processed in the Rmd file thaw_ponds_DA.Rmd.

## Abundance of HMM families on thaw ponds contigs

### Find CDS with Prodigal

```bash
export prodigalOutDir=prodigal
mkdir -p $prodigalOutDir

sbatch -p node -t 20:00:00 -A $snicProj \
-J prodigal -o logs/prodigal.out -e logs/prodigal.err<<'EOF'
#!/bin/bash

module load prodigal

prodigal -i megahit_assembly/thawponds_assembly.contigs.fa \
-p meta \
-a ${prodigalOutDir}/thawponds_assembly.cds.faa \
-d ${prodigalOutDir}/thawponds_assembly.cds.ffn \
-o ${prodigalOutDir}/thawponds_assembly.cds.out \
-m \
-f gff

EOF
```

Split prodigal result file (with protein sequences) in 500K sequence chunks. Outfile names will be `thawponds_assembly.cds.split00001.fa`, `thawponds_assembly.cds.split00002.fa` etc. and will we placed in the same directory as the input. Make a list of the output files so they will be used in the next step.

```bash
scripts/splitSeqFile.py $prodigalOutDir/thawponds_assembly.cds.faa \
fasta \
fasta \
500000

ls $prodigalOutDir/thawponds_assembly.cds.split* > $prodigalOutDir/thawponds_assembly.cds.faa.files
```

### Annotate CDS against PFAM

```bash
cd $wdir
mkdir -p hmmer

sbatch -p node -t 20:00:00 -A $snicProj \
-J hmmer -o logs/hmmer_%a.out -e logs/hmmer_%a.err \
--array=1-$(wc -l < $prodigalOutDir/thawponds_assembly.cds.faa.files)<<'EOF'
#!/bin/bash

module load hmmer

s=$(sed -n "$SLURM_ARRAY_TASK_ID"p $prodigalOutDir/thawponds_assembly.cds.faa.files)
outs=hmmer/$(basename ${s/.faa/})

hmmsearch \
-o ${outs}.hmmer_pfam.out \
--tblout ${outs}.hmmer_pfam.tblout \
--domtblout ${outs}.hmmer_pfam.domtblout \
--pfamtblout ${outs}.hmmer_pfam.pfamtblout \
--cpu 20 \
-E 1e-05 \
data/Pfam-A.hmm \
${s}

EOF
```

Then concatenate all results

```bash
cat hmmer/*split*domtblout > hmmer/thawponds_assembly.cds.split500000.all.hmmer_pfam.domtblout
```

### Re-map hmmscan results on contig annotations

```bash
Rscript --vanilla \
parseHMMscan.R \
-a prodigal/thawponds_assembly.cds.out \
-b hmmer/thawponds_assembly.cds.split500000.all.hmmer_pfam.domtblout \
-o prodigal/thawponds_assembly.cds.hmmer_pfam.gff
```

### Generate bowtie2 dbs for co-assembly of contigs

No actual need to run this as a batch job.

```bash
module load bioinfo-tools
module load bowtie2/2.2.9

bowtie2-build megahit_assembly/thawponds_assembly.contigs.fa thawponds_assembly
mv thawponds_assembly.*.bt2 megahit_assembly
```

### Map reads to contigs

Mapped with bowtie2.
Alignments are name sorted to avoid memory issues with htseq-count.

```bash
export dataDir="$wdir/reads_filtered"
export sampleList="$wdir/data/sampleList"
export outdir="$wdir/alignments"
export logdir="$wdir/logs"

mkdir -p $logdir
mkdir -p $outdir

sbatch -p node -t 20:00:00 -A $snicProj \
-J tp_contig_mappings -o logs/tp_contig_mappings_%a.out -e logs/tp_contig_mappings_%a.err \
--array=1-$(wc -l < $sampleList)<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load bowtie2/2.2.9
#module load htseq

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleList)

cp ${dataDir}/*${sample}*R* ${SNIC_TMP}
cp megahit_assembly/thawponds_assembly.*.bt2 ${SNIC_TMP}

cd ${SNIC_TMP}

bowtie2 -x thawponds_assembly \
-1 trimmed_${sample}_R1_001.fastq \
-2 trimmed_${sample}_R2_001.fastq \
--no-unal \
-S ${sample}.sam \
--threads 20

samtools view -F 4 -bS ${sample}.sam > ${sample}.bam
samtools sort -n -o ${sample}.sorted.name.bam ${sample}.bam

cp ${sample}.sorted.bam $outdir

EOF
```

### Get read counts

Counts for uniquely mapped reads are used for differential expression analyses, counts for all mapped reads are used for descriptive statistics.

#### Counts for uniquely mapped reads (used for differential expression analyses)

```bash
cd $wdir
mkdir -p counts

export sampleList="$wdir/data/sample_list"
export gffFile="thawponds_assembly.cds.hmmer_pfam.gff"

sbatch -A $snicProj -p core -t 6:00:00 \
-J ssort_htseq_thawponds_pfam -o logs/ssort_htseq_thawponds_pfam_%a.out -e logs/ssort_htseq_thawponds_pfam_%a.err \
--array=1-$(wc -l < sampleList)<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load pysam/0.13-python2.7.11
module load htseq
module load samtools

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p sampleList)

cp prodigal/${gffFile} ${SNIC_TMP}
cp alignments/${sample}.sorted.name.bam ${SNIC_TMP} && cd ${SNIC_TMP}

#samtools sort -o ${sample}.sorted.name.bam ${sample}.sorted.bam

#cp ${sample}.sorted.name.bam ${wdir}/alignments

time htseq-count \
-f bam \
-r name \
-s no \
-t Domain \
-i PFAM_ID \
-m union \
--nonunique all \
${sample}.sorted.name.bam \
${gffFile} > ${sample}.pfam.counts.all.out

time htseq-count \
-f bam \
-r name \
-s no \
-t Domain \
-i PFAM_ID \
-m union \
--nonunique none \
${sample}.sorted.name.bam \
${gffFile} > ${sample}.pfam.counts.unique.out

cp *counts*out ${wdir}/counts

EOF

```

