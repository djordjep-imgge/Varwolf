#!/bin/bash
#
# Usage: ./varwolfa.sh <sample name>
#
# Prerequisites packages: fastp, bwa-mem, samtools, docker, GATK image, ANNOVAR, InterVar

set -uex

NAME=$1
REF=refs/hg38.fasta
SNP=db/snps/Homo_sapiens_assembly38.dbsnp138.vcf.gz
INDEL=db/indels/Homo_sapiens_assembly38.known_indels.vcf.gz
DATE=$(date +%y%m%d)
RG=$(echo "@RG\tPL:ILLUMINA\tLB:${DATE}\tSM:${NAME}\tID:1")

# Adapter and quality trimming with fastp producing processed fastq files

for i in {1..4}
do
    fastp -w 12 \
	-i fastq/${NAME}*L00${i}_R1*.fastq.gz -I fastq/${NAME}*L00${i}_R2*.fastq.gz \
	-o fastq/${NAME}_L${i}_R1.fastq.gz -O fastq/${NAME}_L${i}_R2.fastq.gz
done

for i in {1..4}
do
    bwa mem -M -t 12 -R $RG $REF \
	fastq/${NAME}_L${i}_R1.fastq.gz \
	fastq/${NAME}_L${i}_R2.fastq.gz \
	> bams/${NAME}_L${i}.sam
done

# Sorting and then compressing the SAM file using samtools producing a BAM file.

for i in {1..4}
do
    samtools sort -@ 12 -n \
	bams/${NAME}_L${i}.sam \
	> bams/${NAME}_L${i}.bam
done

# Filtering the BAMs file to only keep proper pairs using samtools.

for i in {1..4}
do
    samtools view -@ 12 -f 2 -b \
	bams/${NAME}_L${i}.bam \
	> bams/${NAME}_L${i}pp.bam
done

# Marking duplicate reads in the BAM file using GATK.
docker run -v ~/Genomics:/home -w /home gatk \
	gatk MarkDuplicatesSpark \
		-R $REF \
		-I bams/${NAME}_L1pp.bam \
		-I bams/${NAME}_L2pp.bam \
		-I bams/${NAME}_L3pp.bam \
		-I bams/${NAME}_L4pp.bam \
		-O bams/${NAME}_dd.bam \
		-M bams/${NAME}_metricsdd.txt \
		--spark-master local[12]

# Performing Base Quality Score Recalibration (BQSR) on the BAM file using GATK's BaseRecalibrator and ApplyBQSR.
docker run -v ~/Genomics:/home -w /home gatk \
	gatk BaseRecalibratorSpark \
		-R $REF \
		-I bams/${NAME}_dd.bam \
		-O bams/bqsr_reports/${NAME}_bqsr.report \
		--known-sites $SNP \
		--known-sites $INDEL \
		--spark-master local[12]

docker run -v ~/Genomics:/home -w /home gatk \
	 gatk ApplyBQSRSpark \
		-R $REF \
		-I bams/${NAME}_dd.bam \
		-O bams/${NAME}.bam \
		-bqsr bams/bqsr_reports/${NAME}_bqsr.report \
		--spark-master local[12]

# Variant calling using GATK's HaplotypeCaller to produce a VCF file.
docker run -v ~/Genomics:/home -w /home gatk \
	gatk HaplotypeCaller \
		-R $REF \
		-I bams/${NAME}.bam \
		-O vcfs/${NAME}_raw.vcf \
		-L db/manifest_hg38.bed -ip 10 \
		--dbsnp $SNP

# Hard filtering the variants using GATK's VariantFiltration.
docker run -v ~/Genomics:/home -w /home gatk \
	gatk VariantFiltration \
		-V vcfs/${NAME}_raw.vcf \
		-filter "DP < 5.0" --filter-name "DP5" \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O vcfs/${NAME}.vcf

# Left-normalizing the VCF file and converting it to .avinput for ANNOVAR.
bcftools norm -m-both vcfs/${NAME}.vcf \
	| bcftools norm -f $REF \
	> vcfs/${NAME}_norm.vcf

# VEP annotation plus InterVar ACMG	
	
docker run -v ~/Genomics:/data vep \
	vep --fork 16 --cache \
	-i vcfs/${NAME}_norm.vcf \
	-o tables/${NAME}_VEP_OMIM.vcf \
	--fasta refs/hg38.fasta \
	--dir_cache db \
	--shift_hgvs 1 \
    --pick --pick_order tsl \
    --minimal --failed 1 --vcf --no_progress --everything --xref_refseq --total_length --allele_number --check_existing --no_escape --assembly "GRCh38" --force_overwrite \
	--plugin Phenotypes,file=db/vep_data/Phenotypes.pm_homo_sapiens_109_GRCh38.gvf.gz,exclude_sources='AMDGC&Cancer_Gene_Census&ClinVar&COSMIC&dbGaP&dbVar&DDG2P&DGVa&GEFOS&GIANT&HGMD-PUBLIC&MAGIC&NHGRI-EBI_GWAS_catalog&Orphanet&Teslovich'

bcftools view -H tables/${NAME}_VEP_OMIM.vcf | tr "|" "\t" > tables/${NAME}_VEP_OMIM_tabbed.txt

bcftools view -H vcfs/${NAME}_norm.vcf | cut -f 10 | cut -d ":" -f 1,2 | tr "/" "x" | sed 's/0x1/het/' | sed 's/1x0/het/' | sed 's/1x1/hom/' > tables/${NAME}_genotype.txt

paste <(cut -f 1-8 tables/${NAME}_VEP_OMIM_tabbed.txt) tables/${NAME}_genotype.txt > tables/${NAME}_VEP_vcfinfo.txt

docker run -v ~/Genomics:/data vep \
	vep --fork 16 --cache \
	-i vcfs/${NAME}_norm.vcf \
	-o tables/${NAME}_VEP_Orphanet.vcf \
	--fasta refs/hg38.fasta \
	--dir_cache db \
	--shift_hgvs 1 \
    --pick --pick_order tsl \
    --minimal --failed 1 --vcf --no_progress --everything --xref_refseq --total_length --allele_number --check_existing --no_escape --assembly "GRCh38" --force_overwrite \
	--plugin Phenotypes,file=db/vep_data/Phenotypes.pm_homo_sapiens_109_GRCh38.gvf.gz,exclude_sources='AMDGC&Cancer_Gene_Census&ClinVar&COSMIC&dbGaP&dbVar&DDG2P&DGVa&GEFOS&GIANT&HGMD-PUBLIC&MAGIC&MIM_morbid&Teslovich'

bcftools view -H tables/${NAME}_VEP_Orphanet.vcf | tr "|" "\t" > tables/${NAME}_VEP_Orphanet_tabbed.txt

paste <(cut -f 11 tables/${NAME}_VEP_OMIM_tabbed.txt) <(cut -f 89 tables/${NAME}_VEP_OMIM_tabbed.txt | cut -d "," -f 1) <(cut -f 89 tables/${NAME}_VEP_Orphanet_tabbed.txt | cut -d "," -f 1) > tables/${NAME}_VEP_genelv.txt

cut -f 9,10,15-19,45,52,56,58,64,83 tables/${NAME}_VEP_OMIM_tabbed.txt > tables/${NAME}_VEP_variantlv.txt

docker run -v ~/Genomics:/data vep \
	vep --fork 16 --cache \
	-i vcfs/${NAME}_norm.vcf \
	-o tables/${NAME}_VEP_dbNSFP.vcf \
	--fasta refs/hg38.fasta \
    --dir_cache db \
    --shift_hgvs 1 \
    --pick --pick_order tsl \
    --minimal --failed 1 --vcf --no_progress --everything --xref_refseq --total_length --allele_number --check_existing --no_escape --assembly "GRCh38" --force_overwrite \
	--plugin dbNSFP,db/vep_data/dbNSFP4.3a_grch38.gz,pep_match=0,transcript_match=1,CADD_phred,ExAC_AF,ExAC_NFE_AF,FATHMM_pred,MetaSVM_pred,MutationTaster_pred,PROVEAN_pred,clinvar_clnsig,clinvar_hgvs,clinvar_review,clinvar_trait,clinvar_var_source,rs_dbSNP

bcftools view -H tables/${NAME}_VEP_dbNSFP.vcf | tr "|" "\t" > tables/${NAME}_VEP_dbNSFP_tabbed.txt

paste <(cut -f 47,48,89,92,93 tables/${NAME}_VEP_dbNSFP_tabbed.txt) <(cut -f 94 tables/${NAME}_VEP_dbNSFP_tabbed.txt | cut -d "&" -f 1) <(cut -f 95-100 tables/${NAME}_VEP_dbNSFP_tabbed.txt) > tables/${NAME}_VEP_pred.txt

docker run -v ~/Genomics:/data vep \
	vep --fork 16 --cache \
	-i vcfs/${NAME}_norm.vcf \
	-o tables/${NAME}_VEP_GWAS.vcf \
	--fasta refs/hg38.fasta \
	--dir_cache db \
	--shift_hgvs 1 \
    --pick --pick_order tsl \
    --minimal --failed 1 --vcf --no_progress --everything --xref_refseq --total_length --allele_number --check_existing --no_escape --assembly "GRCh38" --force_overwrite \
	--plugin Phenotypes,file=db/vep_data/Phenotypes.pm_homo_sapiens_109_GRCh38.gvf.gz,exclude_sources='AMDGC&Cancer_Gene_Census&ClinVar&COSMIC&dbGaP&dbVar&DDG2P&DGVa&GEFOS&GIANT&HGMD-PUBLIC&MAGIC&MIM_morbid&Orphanet&Teslovich'

bcftools view -H tables/${NAME}_VEP_GWAS.vcf | tr "|" "\t" > tables/${NAME}_VEP_GWAS_tabbed.txt

cut -f 89 tables/${NAME}_VEP_GWAS_tabbed.txt | cut -d "," -f 1 > tables/${NAME}_VEP_gwas.txt

cd db/InterVar

python3 ./Intervar.py -b hg38 -i ~/Genomics/vcfs/${NAME}_norm.vcf --input_type VCF -otherinfo 100 -d ../annovar/humandb/ -o ~/Genomics/tables/${NAME}_

cd ../..

cut -f 14 tables/${NAME}_.hg38_multianno.txt.intervar | tail -n +2 | sed 's/ PVS1/\tPVS1/' > tables/${NAME}_IV.txt

paste tables/${NAME}_VEP_vcfinfo.txt tables/${NAME}_VEP_genelv.txt tables/${NAME}_VEP_variantlv.txt tables/${NAME}_VEP_pred.txt tables/${NAME}_VEP_gwas.txt tables/${NAME}_IV.txt > tables/${NAME}_headerless.txt

cat db/headervcf.txt tables/${NAME}_headerless.txt > tables/${NAME}.avcf

# Remoivng temporary files

rm -rf fastq/${NAME}*fastq.gz bams/${NAME}*sam bams/${NAME}_* tables/${NAME}_* fastp.*