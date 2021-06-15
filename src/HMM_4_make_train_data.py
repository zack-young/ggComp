# Python 2
import os
import sys
if len(sys.argv) == 3:
	opath = sys.argv[1]
	timetag = sys.argv[2]
else:
	sys.stderr.write('Parameter error\n')
	exit(1)

path = opath + "/" + timetag
file_backend = "homo_undefined_snp_level.sorted"

folders= os.listdir(path)

for folder in folders:
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c= os.listdir(path_c)
		for folder_c in folders_c:
			if folder_c[4:] == "A."+file_backend or folder_c[4:] == "B."+file_backend or folder_c[4:] == "D."+file_backend:
				f = open(path_c + "/" + folder_c)
				lines = f.readlines()
				f.close()
				string_lv = ""
				for line in lines:
					line_ele = line.strip().split("\t")
					if  line_ele[3] == "PHR":
						string_lv += 'B'
					elif  line_ele[3] == "SGR":
						string_lv += 'A'
					else:
						string_lv += 'D'
				f = open(path + "/strings-raw.txt", 'a')
				f.write(string_lv + "\n")
				f.close()