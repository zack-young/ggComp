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

file_backend = "homo_undefined_snp_level.sorted"

new_file_backend = "homo_undefined_snp_level.sorted.smoothed"

states = ["SGR", "PHR", "CNV"]
n_states = 3
n_obs = 3

model = hmm.MultinomialHMM(n_components=n_states, n_iter=20, tol=0.001)
model.startprob_ = numpy.array([0.333,0.333,0.334])
