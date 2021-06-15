
model.transmat_ = model.transmat_ / model.transmat_.sum(axis=1)[:, numpy.newaxis]
model.emissionprob_ = model.emissionprob_ / model.emissionprob_.sum(axis=1)[:, numpy.newaxis]

folders= os.listdir(path)

def mklevelfile(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c= os.listdir(path_c)
		for folder_c in folders_c:
			if folder_c[6:] == file_backend:
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
						raw_seq.append(3)

				raw_seq = numpy.array([raw_seq]).T
				logprob, dec = model.decode(raw_seq, algorithm="viterbi")
				print(path_c + "/" + folder_c)
				print(logprob)
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
				print(count_lev_to_CNV)
				f = open(path_c + "/" + folder_c[:5] + "." + new_file_backend, "w")
				f.writelines(lines)
				f.close()

proce_pool = Pool(processes = int(maxprocess))
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(mklevelfile, args=(folder,)))

proce_pool.close()
proce_pool.join()
