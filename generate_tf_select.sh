for i in Adipose Cecum Colon Cortex Duodenum Hypothalamus Ileum Jejunum Liver Lung Muscle Spleen Stomach
do
bed_file=/home/yinhongwei/pig_yin/merge_peak/ATAC_${i}_Peaks.merge100.bed
findMotifsGenome.pl ${bed_file} /home/yinhongwei/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa ${i}_motif -len 8,10,12 -p 30
mkdir -p ${i}_tf
for j in ${i}_motif/homerResults/motif[0-9].motif ${i}_motif/homerResults/motif[0-9][0-9].motif
do
knowfie=`basename $j|cut -d . -f 1`    
scanMotifGenomeWide.pl $j /home/yinhongwei/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -p 30 > ${i}_tf/${i}.${knowfie}.tf_sites
done
done

for i in Adipose Cecum Colon Cortex Duodenum Hypothalamus Ileum Jejunum Liver Lung Muscle Spleen Stomach
do
for j in ${i}_tf/${i}.motif*.tf_sites
do
awk -F '\t' '{if($6>=5.0){print $3"\t"$4"\t"$6"\t"$2}}' ${j}|awk '{print $4"\t"$1"\t"$2"\t"$3}' >>${i}.motif.bed
done
bed_file=/home/yinhongwei/pig_yin/merge_peak/ATAC_${i}_Peaks.merge100.bed
bedtools intersect -a ${i}.motif.bed -b ${bed_file} -f 1.0 -wa > ${i}.OR5.motif.in.OCR
rm ${i}.motif.bed
done
