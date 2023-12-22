# lightning
###generate two file for mapping and coverage
DIR_IN=$1
DIR_OUT=$2

samtools_path=/vol3/agis/likui_group/yinhongwei/software/miniconda/bin


OUT=$3
filename=${DIR_IN}/${OUT}/${OUT}.rmdup.bam

if [ ! -s ${DIR_OUT}/$OUT.coverage.more ]; then
	echo -en "BAM_file: ${OUT}\n" > ${DIR_OUT}/$OUT.coverage.more
	genome_size=$($samtools_path/samtools view -H $filename | grep -P '^@SQ' | cut -f 3 -d ":" | awk '{sum += $1}END{print sum}')
	echo -en "genome_size: ${genome_size}\n" >> ${DIR_OUT}/$OUT.coverage.more
	eval $($samtools_path/samtools depth -Q 1 $filename | awk '{
		sum += $3; 
		if($3 > 0){
			base += 1
		}; 

		if($3 >= 3){
			base_3+=1
		};

		if($3 >= 5){
			base_5 += 1
		};

		if($3 >= 7){
			base_7 += 1
		};

		if($3 >= 9){
			base_9 += 1
		}; 

		if($3 >= 11){
		base_11+=1
		};
	}

	END{
		printf("total_size=%f; total_base=%f; total_base_3=%f; total_base_5=%f; total_base_7=%f; total_base_9=%f; total_base_11=%f;",sum, base, base_3, base_5, base_7, base_9, base_11)
	}'
	)
	echo -en "total_size: ${total_size}\n">> ${DIR_OUT}/$OUT.coverage.more
	mean_depth=$(echo "$total_size/$genome_size" | bc -l)
	coverage=$(echo "$total_base/$genome_size" | bc -l)
	echo -en "mean_depth: ${mean_depth}\n">> ${DIR_OUT}/$OUT.coverage.more
	echo -en "coverage: ${coverage}\n">> ${DIR_OUT}/$OUT.coverage.more

	coverage_3=$(echo "$total_base_3/$genome_size" | bc -l)
	coverage_5=$(echo "$total_base_5/$genome_size" | bc -l)
	coverage_7=$(echo "$total_base_7/$genome_size" | bc -l)
	coverage_9=$(echo "$total_base_9/$genome_size" | bc -l)
	coverage_11=$(echo "$total_base_11/$genome_size" | bc -l)

	echo -en "coverage_3: ${coverage_3}\n">> ${DIR_OUT}/$OUT.coverage.more
	echo -en "coverage_5: ${coverage_5}\n">> ${DIR_OUT}/$OUT.coverage.more
	echo -en "coverage_7: ${coverage_7}\n">> ${DIR_OUT}/$OUT.coverage.more
	echo -en "coverage_9: ${coverage_9}\n">> ${DIR_OUT}/$OUT.coverage.more
	echo -en "coverage_11: ${coverage_11}\n">> ${DIR_OUT}/$OUT.coverage.more
fi &

if [ ! -s ${DIR_OUT}/$OUT.mappingsum  ]; then
	$samtools_path/samtools flagstat $filename > ${DIR_OUT}/$OUT.mappingsum
fi &
wait

