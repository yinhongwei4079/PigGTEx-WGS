for i in Adipose Cecum Cerebellum Colon Cortex Duodenum Hypothalamus Ileum Jejunum Liver Lung Muscle Spleen Stomach
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
