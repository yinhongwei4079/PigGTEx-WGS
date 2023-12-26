basenji_data.py -d 0.1 --local -o bji_Adipose_13kb_${i} -p 12 -r 64 -w 128 --crop 16384 -l 131072 -t 0.10 -v 0.10  susScr11.fa Adipose.pig.txt
basenji_train.py -o models/Adipose_13kb_${i} models/params_13kb_0.1.json bji_Adipose_13kb_${i}
basenji_test.py --ai 0,1,2 -o output/bji_Adipose_13kb_${i} --rc --shifts "1,0,-1" models/params_13kb_0.1.json models/Adipose_13kb_${i}/model_best.h5 ./bji_Adipose_13kb_${i}
