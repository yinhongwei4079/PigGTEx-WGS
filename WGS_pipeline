####### clean data    #################
####### trimmomatic -version 0.39 #####
trimmomatic PE -threads 4 ${i} ${id}_R1.fastq.gz ${id}_R2.fastq.gz ${id}_R1.paired.fastq.gz ${id}_R1.unpaired.fastq.gz ${id}_R2.paired.fastq.gz ${id}_R2.unpaired.fastq.gz MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20

#######  bwa align and sorted ######################
#######  bwa -version 0.7.5a-r405 ##################
#######  samtools -version 1.9 #####################
bwa mem -t 4 -M -R "@RG\tID:${id}\tSM:${id}" $ref ${id}_R1.paired.fastq.gz ${id}_R2.paired.fastq.gz | samtools view -Sb -o ${id}.unsort.bam -
samtools sort -@ 4 -O bam -o ${id}.sort.bam  ${id}.unsort.bam

#######  MarkDuplicates and generate index  ########
#######  GATK -version v4.1.4.1 ####################
gatk MarkDuplicates -I ${id}.sort.bam -M ${id}.marked_dup_metrics.txt -O ${id}.rmdup.bam
samtools index ${id}.rmdup.bam

#######  generate gvcf file for indi in chr region #######
#######  we divided the whole genome into 71 segments in 40M steps ########
gatk HaplotypeCaller --emit-ref-confidence GVCF -R /vol3/agis/likui_group/yinhongwei/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -I ${id}.rmdup.bam -L ${chr_region}.intervals -O ${id}.${chr_region}.g.vcf.gz

####### CombineGVCFs for all samples and GenotypeGVCFs in single region(more than 1000 samples) ####
gatk GenomicsDBImport -R ${ref} --sample-name-map sample_map --genomicsdb-workspace-path database_${chr_region} --intervals ${chr_region}.intervals --tmp-dir=/vol3/agis/likui_group/yinhongwei/gte_pig/gatk_tmp --reader-threads 20 
gatk GenotypeGVCFs -R ${ref} -V gendb://database_${chr_region} -O ${output}/${outnames}.${chr_region}.vcf.gz

####### combine genotype in single chr and SelectVariants  ##############
gatk SelectVariants -R $ref -V ${output}/${outnames}.${chr_region}.vcf.gz -select-type SNP -O ${output}/${outnames}.${chr_region}.raw.SNP.vcf.gz

####### mark the hard filter for snp in gatk #################################
gatk VariantFiltration -R $ref -V ${output}/${outnames}.${chr_region}.raw.SNP.vcf.gz -filter \"QD < 2.0\" --filter-name \"QD\" -filter \"MQ < 40.0\" --filter-name \"MQ\" -filter \"FS > 60.0\" --filter-name \"FS\" -filter \"SOR > 3.0\" --filter-name \"SQR\" -filter \"MQRankSum < -12.5\" --filter-name \"MQRS\" -filter \"ReadPosRankSum < -8.0\" --filter-name \"RPRS\" -O ${output}/${outnames}._${chr_region}.SNP.mark.vcf.gz

####### select snps which passed the hard filter and was the binary snp for autosome using bcftools ##############
####### bcftools -version 1.9 #####################
gatk MergeVcfs -I /vol3/agis/likui_group/yinhongwei/fst_5/chr_list/chr${j}.list -O ${output}/${outnames}.chr${j}.SNP.mark.vcf.gz
bcftools view ${output}/${outnames}.${j}.SNP.mark.vcf.gz --min-af 0.01:minor -e 'F_MISSING>0.9' -m 2 -M 2 -v snps -f PASS -O z -o ${output}/${outnames}.${j}.filter.maf0.01.missing.binary.vcf.gz

################## phasing #################################
############### beagle.18May20.d20.jar ######################
java -Xmx100g -jar /vol3/agis/likui_group/yinhongwei/software/beagle/beagle.18May20.d20.jar gt=${outnames}.${j}.filter.maf0.01.missing.binary.vcf.gz out=${outnames}.${j}.filter.maf0.01.missing.binary.phased seed=9823 nthreads=20

