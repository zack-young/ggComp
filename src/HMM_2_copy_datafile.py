# Python 2
import os
import sys
from multiprocessing import Pool, Process, Manager

if len(sys.argv) == 6:
	path = sys.argv[1]
	opath = sys.argv[2]
	fpath = sys.argv[3]
	timetag = sys.argv[4]
	maxprocess = sys.argv[5]
else:
	sys.stderr.write('Parameter error\n')
	exit(1)

with open(fpath) as f:
	folders = f.readlines()

os.system("mkdir -p " + opath + "/" + timetag)


def copy_file(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c = os.listdir(path_c)
		os.makedirs(opath + "/" + timetag + "/" + folder)
		for folder_c in folders_c:
			#if folder_c[5:] == ".homo_undefined_snp_level":
			shell = "ln -s " + path_c + "/" + folder_c + " " + opath + "/" + timetag + "/" + folder + "/"
			os.system(shell)
	return 0


proce_pool = Pool(processes = int(maxprocess))
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(copy_file, args=(folder.strip(),)))

proce_pool.close()
proce_pool.join()
				
