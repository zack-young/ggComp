import sys
import numpy
from hmmlearn import hmm
from multiprocessing import Pool, Process, Manager

if len(sys.argv) == 4:
	opath = sys.argv[1]
	timetag = sys.argv[2]
	n_itera = sys.argv[3]
else:
	sys.stderr.write('Parameter error\n')
	exit(1)
#path = sys.argv[1]

all_lines = []

f = open(opath + "/" + timetag +  "/strings-raw.txt")
all_lines.extend(f.readlines())
f.close()

count = 0.0
count_AA = 0.0
count_AB = 0.0
count_AD = 0.0
count_BA = 0.0
count_BB = 0.0
count_BD = 0.0
count_DA = 0.0
count_DB = 0.0
count_DD = 0.0
count_strange = 0.0

for line in all_lines:
	line = line.strip()
	for i in range(0, len(line)-1):
		if line[i:i+2] == "AA":
			count += 1
			count_AA += 1
		elif line[i:i+2] == "AB":
			count += 1
			count_AB += 1
		elif line[i:i+2] == "AD":
			count += 1
			count_AD += 1
		elif line[i:i+2] == "BA":
			count += 1
			count_BA += 1
		elif line[i:i+2] == "BB":
			count += 1
			count_BB += 1
		elif line[i:i+2] == "BD":
			count += 1
			count_BD += 1
		elif line[i:i+2] == "DA":
			count += 1
			count_DA += 1
		elif line[i:i+2] == "DB":
			count += 1
			count_DB += 1
		elif line[i:i+2] == "DD":
			count += 1
			count_DD += 1
		else:
			count_strange += 1

Prob_AA = count_AA / (count_AA+count_AB+count_AD)
Prob_AB = count_AB / (count_AA+count_AB+count_AD)
Prob_AD = count_AD / (count_AA+count_AB+count_AD)

Prob_BA = count_BA / (count_BA+count_BB+count_BD)
Prob_BB = count_BB / (count_BA+count_BB+count_BD)
Prob_BD = count_BD / (count_BA+count_BB+count_BD)

Prob_DA = count_DA / (count_DA+count_DB+count_DD)
Prob_DB = count_DB / (count_DA+count_DB+count_DD)
Prob_DD = count_DD / (count_DA+count_DB+count_DD)

def render_output_prob(curr_p, all_p):
	if curr_p == min(all_p):
		curr_p = 1.0 - (float(str(all_p[0])[:6]) + float(str(all_p[1])[:6]) + float(str(all_p[2])[:6]) - float(str(min(all_p))[:6]))
		return float(str(curr_p))
	else:
		return float(str(curr_p)[:6])

AA = render_output_prob(Prob_AA, [Prob_AA, Prob_AB,Prob_AD])
AB = render_output_prob(Prob_AB, [Prob_AA, Prob_AB, Prob_AD])
AD = render_output_prob(Prob_AD, [Prob_AA, Prob_AB, Prob_AD])
BA = render_output_prob(Prob_BA, [Prob_BA, Prob_BB, Prob_BD])
BB = render_output_prob(Prob_BB, [Prob_BA, Prob_BB, Prob_BD])
BD = render_output_prob(Prob_BD, [Prob_BA, Prob_BB, Prob_BD])
DA = render_output_prob(Prob_DA, [Prob_DA, Prob_DB, Prob_DD])
DB = render_output_prob(Prob_DB, [Prob_DA, Prob_DB, Prob_DD])
DD = render_output_prob(Prob_DD, [Prob_DA, Prob_DB, Prob_DD])

def trans_data(seq):
	seq = seq.strip()
	seq_data = []
	for i in range(len(seq)):
		if seq[i] == "A":
			seq_data.append([0])
		if seq[i] == "B":
			seq_data.append([1])
		if seq[i] == "D":
			seq_data.append([2])
	return seq_data

states = ["a", "b", "d"]
n_states = 3
obs = ["A", "B", "D"]
n_obs = 3

start_arr = numpy.array([0.3333,0.3333,0.3334])
trans_mat = numpy.array([[AA,AB,AD],[BA,BB,BD],[DA,DB,DD]])
emiss_mat = numpy.array([[AA,AB,AD],[BA,BB,BD],[DA,DB,DD]])

#print(trans_mat)

model = hmm.MultinomialHMM(n_components = n_states, verbose = True, n_iter = int(n_itera), tol = 0.001, init_params = "")
model.startprob_ = start_arr
model.transmat_ = trans_mat
model.emissionprob_ = emiss_mat

model.transmat_ = model.transmat_ / model.transmat_.sum(axis=1)[:, numpy.newaxis]
model.emissionprob_ = model.emissionprob_ / model.emissionprob_.sum(axis=1)[:, numpy.newaxis]

#path = opath + "/" + timetag +  "/strings-raw.txt"

#print(path)

#print(type(emiss_mat))

#f = open(path)
lines = all_lines
#f.close()

train_data = []

proce_pool = Pool()
mulres = []
for line in lines:
	mulres.append(proce_pool.apply_async(trans_data, args=(line,)))

proce_pool.close()
proce_pool.join()

for res in mulres:
	train_data.append(res.get())
	#print res.get()[0]


del proce_pool
del mulres

#print(type(train_data))
#print(type(train_data[0]))
train_data = numpy.concatenate(train_data)
#print(type(train_data))
#print(train_data)

model.fit(train_data)

T = model.transmat_.tolist()
E = model.emissionprob_.tolist()
E[0][2] = 0
E[1][2] = 0
E[2][0] = 0
E[2][1] = 0
print('''model.transmat_ = numpy.array(''' + str(T) + ''')''')
print('''model.emissionprob_ = numpy.array(''' + str(E) + ''')''')



#print("startprob")
#print(model.startprob_)
#print("transmat_")
#print(model.transmat_)
#print("emissionprob_")
#print(model.emissionprob_)
#print("train_data_score")
#print(model.score(train_data))