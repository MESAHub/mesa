import numpy as np

fi = open('ion.data','r')


data = np.zeros(28) # H-Ni
for i,line in enumerate(fi):
	if i <= 2:
		continue
	line = line.strip().split(' ')
	line = list(float(l) for l in line if len(l) > 0)
	if line[1] > data[int(line[0])-1]:
		data[int(line[0])-1] = line[1] # Pull out max binding energy

for d in data:
	print(d)