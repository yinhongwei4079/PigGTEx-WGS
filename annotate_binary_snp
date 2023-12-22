############################### annotate by coovar #############################################
coovar=/vol3/agis/likui_group/yinhongwei/software/CooVar/coovar.pl
gtf=/vol3/agis/likui_group/yinhongwei/software/snpEff/data/Sscrofa11.1/genes.gtf
seq=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/03_coovar/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
for chr in {1..18}
do
   echo "cd /vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/03_coovar/" > ./codes/coovar.chr${chr}.sh
   echo "perl ${coovar} -e ${gtf} -r ${seq} -v GTE.1602.chr${chr}.indel.snp.filter.vcf.gz --feature_source -o ./res --no_contig_sum|tee ./res/coovar.chr${chr}.log" >>./codes/coovar.chr${chr}.sh
done
############################### annotate by SNPEFF ##############################################
snpeff=/vol3/agis/likui_group/yinhongwei/software/snpEff/snpEff.jar
input=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/01_binary_vcf
java -jar ${snpeff} Sscrofa11.1  ${input}/GTE.1602.chr1_18.filter.binary.vcf -interval /vol3/agis/likui_group/yinhongwei/software/snpEff/data/Sscrofa11.1/regulation.sort.bed> GTE.1602.chr1_18.filter.binary.ann.vcf
############################### annotate by VEP and UTRannotator ################################
input=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/01_binary_vcf
vep -i ${input}/GTE.1602.chr1_18.filter.binary.vcf --cache --everything --species sus_scrofa --offline --dir /vol3/agis/likui_group/yinhongwei/pig/Sscrofa11.1 \
--fasta /vol3/agis/likui_group/yinhongwei/pig/Sscrofa11.1/sus_scrofa/103_Sscrofa11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz --stats_text \
--stats_file GTE.1602.chr1_18.filter.binary.summary.html --vcf -o GTE.1602.chr1_18.filter.binary.vep.vcf --fork 20 --plugin UTRannotator
############################### annotate by pCADD ###############################################
zcat /vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/01_binary_vcf/GTE.1602.chr6.filter.binary.recode.vcf.gz| awk '{if ( $1~/^#/){print $0} else  {print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' > GTE.1602.chr6.filter.binary.recode.vcf
python2 /vol3/agis/likui_group/yinhongwei/software/pcadd/pcadd-scripts-data-master/pCADD_annotate_SNP.py -i GTE.1602.chr6.filter.binary.recode.vcf -c 1,2,4,5 -p /vol3/agis/likui_group/yinhongwei/software/pcadd/snp_data/6_pCADD-PHRED-scores.tsv.gz
