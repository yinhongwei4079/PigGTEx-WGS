ancestral=/vol3/agis/likui_group/yinhongwei/gte_reference/AIM/ancestral_2.bed.gz
hdr=/vol3/agis/likui_group/yinhongwei/gte_reference/06_prun/hdr.txt
bcftools annotate -a ${ancestral} -c CHROM,POS,REF,ALT,INFO/AA -h ${hdr} -Oz -o GTE.1602.chr1_18.filter.binary.vcf.AA.gz GTE.1602.chr1_18.filter.binary.vcf.gz --threads 40
