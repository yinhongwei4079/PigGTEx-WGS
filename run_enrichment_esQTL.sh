cd /vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/24_basenji/qtl_snp/eqtl

for i in Adipose Artery Blastocyst Blastomere Blood Brain Cartilage Colon Duodenum Embryo Fetal_thymus Frontal_cortex Heart Hypothalamus Ileum Jejunum Kidney Large_intestine Liver Lung Lymph_node Macrophage Milk Morula Muscle Oocyte Ovary Pituitary Placenta Small_intestine Spleen Synovial_membrane Testis Uterus
do

true_qtl=/vol3/agis/likui_group/yinhongwei/pig/eQTL/eQTL_T/${i}.cis_qtl_true.txt
for chr in {1..18}
do
        qtl_file=/vol3/agis/likui_group/yinhongwei/pig/eQTL/${i}.cis_qtl_pairs.${chr}.txt.gz
        zcat ${qtl_file}|awk '{if($1 !="phenotype_id"){print $2}}'>>${i}.tmp_id
done

cat ${i}.tmp_id|sort|uniq|awk -F '_' '{print "chr"$1"\t"$2"\t"$2}' > ${i}.tmp_id.bed
awk '{print $2}' ${true_qtl}|sort|uniq|awk -F '_' '{print "chr"$1"\t"$2"\t"$2}'> ${i}.true_qtl.bed
all_qtl_snp=`cat ${i}.tmp_id.bed|wc -l`
sig_qtl_snp=`cat ${i}.true_qtl.bed|wc -l`


for l in Adipose Cecum Cerebellum Colon Cortex Duodenum Hypothalamus Ileum Jejunum Liver Lung Muscle Spleen Stomach
do
for j in `seq 0.10 0.10 1.00`
do
sad_region_bed=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/24_basenji/sad_cat/${l}.aver.sad.${j}.bed
sad_region_bed_all_nums=`bedtools intersect -a ${sad_region_bed} -b ${i}.tmp_id.bed -wa -f 1.0|sort|uniq|wc -l`
sad_region_bed_sig_nums=`bedtools intersect -a ${sad_region_bed} -b ${i}.true_qtl.bed -wa -f 1.0|sort|uniq|wc -l`

echo ${i} ${l} ${j}  ${all_qtl_snp} ${sad_region_bed_all_nums} ${sig_qtl_snp} ${sad_region_bed_sig_nums} >>eqtl.snp.enrich.txt
done
done
