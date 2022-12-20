import sys
chr = sys.argv[1]
def get_dict(filename):
	f = open(filename,"r")
	f_dict = {}
	for i in f:
		if i.startswith("CHROM"):
			continue
		else:
			i = i.strip().split()
			f_dict[i[0] + "_" + i[1]] = i[4] + " " + i[5]
	return f_dict

Babyrousa = get_dict("Babyrousa.group.chr" + chr + ".frq")
Phacochoerus = get_dict("Phacochoerus.group.chr" + chr + ".frq")
Porcula = get_dict("Porcula.group.chr" + chr + ".frq")
Potamochoerus = get_dict("Potamochoerus.group.chr" + chr + ".frq")
Sus = get_dict("Sus.group.chr" + chr + ".frq")
ancestel = open(sys.argv[2] + chr + ".bed" ,"w+")

for chrpos in Sus.keys():
	if Sus[chrpos].split(" ")[0].split(":")[1] == "0" or "1":
		if Sus[chrpos].split(" ")[0].split(":")[1] == "1":
			ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Sus[chrpos].split(" ")[0].split(":")[0] + "\n")
		else:
			ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Sus[chrpos].split(" ")[1].split(":")[0] + "\n")
	else:
		if Porcula[chrpos].split(" ")[0].split(":")[1] == "0" or "1":
			if Porcula[chrpos].split(" ")[0].split(":")[1] == "1":
				ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Porcula[chrpos].split(" ")[0].split(":")[0] + "\n")
			else:
				ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Porcula[chrpos].split(" ")[1].split(":")[0] + "\n")
		else:
			if Phacochoerus[chrpos].split(" ")[0].split(":")[1] == "0" or "1":
				if Phacochoerus[chrpos].split(" ")[0].split(":")[1] == "1":
					ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Phacochoerus[chrpos].split(" ")[0].split(":")[0] + "\n")
				else:
					ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Phacochoerus[chrpos].split(" ")[1].split(":")[0] + "\n")
			else:
				if Potamochoerus[chrpos].split(" ")[0].split(":")[1] == "0" or "1":
					if Potamochoerus[chrpos].split(" ")[0].split(":")[1] == "1":
						ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Potamochoerus[chrpos].split(" ")[0].split(":")[0] + "\n")
					else:
						ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Potamochoerus[chrpos].split(" ")[1].split(":")[0] + "\n")
				else:
					if Babyrousa[chrpos].split(" ")[0].split(":")[1] == "0" or "1":
						if  Babyrousa[chrpos].split(" ")[0].split(":")[1] == "1":
							ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Babyrousa[chrpos].split(" ")[0].split(":")[0] + "\n")
						else:
							ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + Babyrousa[chrpos].split(" ")[1].split(":")[0] + "\n")
					else:
						ancestel.write(chrpos.split("_")[0] + "\t" + chrpos.split("_")[1] + "\t" + "unknow" + "\n")
