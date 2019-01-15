#!/bin/bash

#SBATCH -p core -t 5-00:00:00 -A snic2018-3-22
#SBATCH -J testMetaDomain_chunk1M -o logs/testMetaDomain_chunk1M.out -e logs/testMetaDomain_chunk1M.err<<'EOF'

export PATH=/home/domeni/glob/MetaDomain:$PATH

outDir=$1
sample_name=$2
chunks=$3
direction=$4
dataDir=$5
inDir=$6

cp ${dataDir}/Pfam-A.hmm ${SNIC_TMP}
cp ${inDir}/${sample_name}.${chunks}.${direction}.fasta ${SNIC_TMP}
mkdir -p ${outDir}

#cp $wdir/filtered_reads/${sample_name}.${chunk}.2.fasta ${SNIC_TMP}
#cp $wdir/filtered_reads/C-2_S28_L006.2.fasta ${SNIC_TMP}
#head -400000 $wdir/filtered_reads/C-2_S28_L006.2.fasta > ${SNIC_TMP}/C-2_S28_L006.2.fasta
grep -v "^#" ${dataDir}/carbohydrate_degradation_pfams_tveit.csv | awk 'BEGIN{FS=","}NR>1{print $3}' > ${SNIC_TMP}/carbo_hmm

export outfolder_tmp=metadomain_1M_${sample_name}_${direction}_${chunks}

cd ${SNIC_TMP}

MetaDomain.sh \
-m Pfam-A.hmm \
-l carbo_hmm \
-f ${sample_name}.${chunks}.${direction}.fasta \
-o ${outfolder_tmp}

cp -R ${outfolder_tmp} ${outDir}
