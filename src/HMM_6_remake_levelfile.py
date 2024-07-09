#python3
import os
import sys
import numpy
from hmmlearn import hmm
from sklearn.preprocessing import normalize
from multiprocessing import Pool, Process, Manager

if len(sys.argv) == 4:
	opath = sys.argv[1]
	timetag = sys.argv[2]
	maxprocess = sys.argv[3]
else:
	sys.stderr.write('Parameter error\n')
	exit(1)


path = opath + "/" + timetag

file_backend = "sorted"

new_file_backend = "smoothed"

states = ["SGR", "PHR", "CNV"]
n_states = 3
n_obs = 3

model = hmm.MultinomialHMM(n_components=n_states, n_iter=20, tol=0.001)
model.startprob_ = numpy.array([0.333,0.333,0.334])
model.transmat_ = numpy.array([[0.96, 0.03,0.006],[0.0009, 0.99, 0.01],[0.003, 0.16, 0.84]])
model.emissionprob_ = numpy.array([[0.96, 0.03, 0], [0.0009, 0.99, 0], [0, 0, 1]])

model.transmat_ = model.transmat_ / model.transmat_.sum(axis=1)[:, numpy.newaxis]
model.emissionprob_ = model.emissionprob_ / model.emissionprob_.sum(axis=1)[:, numpy.newaxis]

folders = os.listdir(path)

def mklevelfile(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c = os.listdir(path_c)
		sorted_files = [file for file in folders_c if file.endswith('.sorted')]
		for folder_c in sorted_files:
			#if folder_c[6:] == file_backend:
			f = open(path_c + "/" + folder_c)
			lines = f.readlines()
			f.close()

			raw_seq = []
			lines_id = []
			for i in range(len(lines)):
				line = lines[i].strip().split("\t")
				if line[3] == "SGR":
					lines_id.append(i)
					raw_seq.append(0)
				elif line[3] == "PHR":
					lines_id.append(i)
					raw_seq.append(1)
				else:
					lines_id.append(i)
					raw_seq.append(2)
			
			raw_seq = numpy.array([raw_seq]).T
			logprob, dec = model.decode(raw_seq, algorithm="viterbi")
			#print(path_c + "/" + folder_c)
			#print(logprob)
			count_lev_to_CNV = 0
			for i in range(len(lines_id)):
				state = states[dec[i]]
				line = lines[lines_id[i]].strip().split("\t")
				if line[3] == "SGR":
					if state == "CNV":
						count_lev_to_CNV += 1
					else:
						line[3] = state
				elif line[3] == "PHR":
					if state == "CNV":
						count_lev_to_CNV += 1
					else:
						line[3] = state
				line = "\t".join(line)
				line += "\n"
				lines[lines_id[i]] = line
			#print(count_lev_to_CNV)
			f = open(path_c + "/" + folder_c + "." + new_file_backend, "w")
			f.writelines(lines)
			f.close()

proce_pool = Pool(processes = int(maxprocess))
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(mklevelfile, args=(folder,)))
	#mklevelfile(folder)

proce_pool.close()
proce_pool.join()
