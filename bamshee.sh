#!/bin/bash
#
#SBATCH -J bamshee
#SBATCH --output /lustre/imgge/PharmGenHub/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 1-00:00:00

module load fastp
module load bwa
module load samtools
module load gatk

set -eux

SAMPLE=$1
WDIR="/lustre/imgge/PharmGenHub"
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
#INTERVALS="/lustre/imgge/db/reference/hg38/wgs_analysis_intervals"
FLOWCELL=$(zcat "${WDIR}"/fastq/"${SAMPLE}"*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 3)
PLATFORM="ILLUMINA"
LIBRARY=$(zcat ${WDIR}/fastq/${SAMPLE}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
RG="@RG\tID:${FLOWCELL}\tPL:${PLATFORM}\tLB:${LIBRARY}\tSM:${SAMPLE}"
THREADS=$SLURM_CPUS_PER_TASK
#TMPDIR=$WDIR/tmp
LANE=4
LANES=$(seq 1 $LANE)

for i in $LANES
do
    fastp \
		-w $((THREADS / 4)) \
	    -i ${WDIR}/fastq/${SAMPLE}*L00${i}_R1*.fastq.gz \
        -I ${WDIR}/fastq/${SAMPLE}*L00${i}_R2*.fastq.gz \
        --stdout \
    | bwa mem \
		-t ${THREADS} \
		-M -p -R ${RG} \
        ${REF} - \
    | samtools sort \
		-@ $((THREADS / 4)) -n \
        > bams/${SAMPLE}_L${i}.bam
	BAMSHARDS+=$(echo -n "-I ${WDIR}/bams/${SAMPLE}_L${i}.bam ")
done

echo "Finished aligning ${SAMPLE}"

gatk MarkDuplicatesSpark \
	-R ${REF} \
	${BAMSHARDS} \
	-O ${WDIR}/bams/${SAMPLE}_dd.bam \
	-M ${WDIR}/bams/metrics/${SAMPLE}_mdmetrics.txt \
	--spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished writing ${SAMPLE}"

gatk BQSRPipelineSpark \
	-R ${REF} \
	-I ${WDIR}/bams/${SAMPLE}_dd.bam \
	-O ${WDIR}/bams/${SAMPLE}.bam \
	--known-sites ${DBSNP} \
	--spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished recalibrating ${SAMPLE}"

for CHR in {1..22} X Y M
do
	GVCFSHARDS+=$(echo -n "-I ${WDIR}/gvcfs/${SAMPLE}_chr${CHR}.g.vcf.gz ")
    gatk HaplotypeCaller \
        -R ${REF} \
        -L chr${CHR} \
        -I ${WDIR}/bams/${SAMPLE}.bam \
        -O ${WDIR}/gvcfs/${SAMPLE}_chr${CHR}.g.vcf.gz \
        -dbsnp ${DBSNP} \
    	-ERC GVCF &
done

wait

gatk MergeVcfs \
    ${GVCFSHARDS} \
    -O ${WDIR}/gvcfs/${SAMPLE}.g.vcf.gz

echo "Finished calling ${SAMPLE}"

rm \
	${WDIR}/fastp* \
	${WDIR}/bams/*${SAMPLE}*_* \
	${WDIR}/gvcfs/*${SAMPLE}*chr*