cd /vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/24_basenji/gwas
sad_bed=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/24_basenji/sad_cat/Duodenum.aver.sad.0.30.bed
gwas_bed=/vol3/agis/likui_group/yinhongwei/gte_reference/losfunction/22_gwas/MainClass/M_GD.fgwas.bed
bedtools intersect -a ${gwas_bed} -b ${sad_bed} -f 1.0 -loj |awk '{if($9 == -1){print $4"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7"\t"0} else {print $4"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7"\t"1}}' >fgwas/M_GD.0.30.Duodenum.txt
sed -i "1i SNPID\tCHR\tPOS\tF\tZ\tN\tcoding_exon" fgwas/M_GD.0.30.Duodenum.txt
gzip -f fgwas/M_GD.0.30.Duodenum.txt
fgwas -i  fgwas/M_GD.0.30.Duodenum.txt.gz  -w coding_exon -o fgwas_res/M_GD.0.30.Duodenum
