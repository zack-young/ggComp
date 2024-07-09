# Python 2
import os
import sys
from multiprocessing import Pool, Process, Manager

if len(sys.argv) == 4:
	opath = sys.argv[1]
	timetag = sys.argv[2]
	maxprocess = sys.argv[3]
else:
	sys.stderr.write('Parameter error\n')
	exit(1)

path = opath + "/" + timetag

#file1_backend = "homo_undefined_snp_level"

folders= os.listdir(path)

def process_data(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c = os.listdir(path_c)
		for files_c in folders_c:
			print(files_c)
			shell = "sort -n -k 2 " + path_c + "/" + files_c + " > " + path_c + "/" + files_c + ".sorted"
			os.system(shell)
			f = open(path_c + "/" + files_c + ".sorted")
			lines1 = f.readlines()
			f.close()
		
#		for i in range(0,7):
#			for j in ["A", "B", "D"]:
#				shell = "sort -n -k 2 " + path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + " > " + path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + ".sorted"
#				os.system(shell)
#				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + ".sorted")
#				lines1 = f.readlines()
#				f.close()
	return 0


proce_pool = Pool(processes = int(maxprocess))
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(process_data, args=(folder,)))

proce_pool.close()
proce_pool.join()
