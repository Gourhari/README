import numpy as np
import random
import math
import re
import sys
import os
import subprocess
import time

atomCount = 5

def runGaussian(filename):
	cmd = "g09 " + filename + ".gjf"
	os.system(cmd)
	#returned_value = subprocess.call(cmd, shell=True)

def loadArray(filename):
	readFile = open(filename + ".txt", "r")
	rline = readFile.readlines()
	print("Range: " + str(len(rline)))
	struct = []
	for i in range(atomCount):
		struct.append([])
		words = rline[i].split()
		struct[i].append(float(words[0]))
		struct[i].append(float(words[1]))
		struct[i].append(float(words[2]))
	readFile.close()
	return struct
def arrayToTxt(struct, filename):
	writeFile = open(filename + ".txt", "w")
	for i in range(atomCount):
		for j in range(3):
			writeFile.write(str(struct[i][j]))
			writeFile.write("\t")
		writeFile.write("\n")
	writeFile.close()

def genGJF(filename):
	readFile = open(filename + ".txt", "r")
	writeFile = open(filename + ".gjf", "w")
	lines_of_text = ["%nprocshared=2\n", "%mem=2GB\n", "%nosave\n", "# b3lyp/6-311+g(d,p)\n", "\n", "B5 struct\n", "\n", "0 2\n"]
	writeFile.writelines(lines_of_text)
	rline = readFile.readlines()
	for line in rline:
		writeFile.write("B\t")
		writeFile.write(line)
	writeFile.write("\n\n")
	readFile.close()
	writeFile.close()
def readEnergy(filename):
	readFile = open(filename + ".log", "r")
	rline = readFile.readlines()
	for i in range (len(rline)):
		if "E(UB3LYP)" in rline[i]:
			words = rline[i].split()
			readFile.close()
			return float(words[4])
	readFile.close()
	return 9999.00
def readTimeTaken(filename):
	readFile = open(filename + ".log", "r")
	rline = readFile.readlines()
	for i in range (len(rline)):
		if "Job cpu time" in rline[i]:
			words = rline[i].split()
			readFile.close()
			return (float(words[7]) * 60 + float(words[9]))
	readFile.close()

filename = sys.argv[1]
genGJF(filename)
runGaussian(filename)
#print(readEnergy(filename))
struct = loadArray(filename)

iter = 30
trailsPerIter = 20


na = 0.0
p1 = 0.7
p50 = 0.001

# Initial temperature
t1 = -1.0/math.log(p1)
#print("t1:", t1)
# Final temperature
t50 = -1.0/math.log(p50)
#print("t50:", t50)


# Fractional reduction every cycle
frac = (t50/t1)**(1.0/(iter-1.0))

na = na + 1.0
structNew = []
for i in range(atomCount):
	structNew.append([])
	structNew[i].append(struct[i][0])
	structNew[i].append(struct[i][1])
	structNew[i].append(struct[i][2])

structBest = []
for i in range(atomCount):
	structBest.append([])
	structBest[i].append(struct[i][0])
	structBest[i].append(struct[i][1])
	structBest[i].append(struct[i][2])

initialCost = readEnergy(filename)
costBest = initialCost
newCost = initialCost

t = t1
timeTotal = 0.0
DeltaE_avg = 0.0
for i in range(iter):
	print('Cycle: ' + str(i) + ' with Temperature: ' + str(t))
	for j in range(trailsPerIter):
		print('\tTrial number: ' + str(j))
		# Generate new trial points
		for k in range(atomCount):
			structNew[k][0] = structNew[k][0] + random.random() - 0.5
			structNew[k][1] = structNew[k][1] + random.random() - 0.5
			structNew[k][2] = structNew[k][2] + random.random() - 0.5
		arrayToTxt(structNew, filename + "new")
		genGJF(filename + "new")
		runGaussian(filename + "new")
		print("\tNew Energy:\t" + str(readEnergy(filename + "new")))
		timeTotal = timeTotal + readTimeTaken(filename + "new")
		print("\tTime Taken:\t" + str(readTimeTaken(filename + "new")))
		DeltaE = abs(readEnergy(filename + "new") - costBest)
		if(readEnergy(filename + "new") > costBest):
			if (i==0 and j==0): 
				DeltaE_avg = DeltaE
			p = math.exp(-DeltaE/(DeltaE_avg * t))
			if (random.random() < p):
				accept = True
			else:
				accept = False
		else:
			accept = True
		if(accept == True):
			for k in range(atomCount):
				structBest[k][0] = structNew[k][0]
				structBest[k][1] = structNew[k][1]
				structBest[k][2] = structNew[k][2]
				costBest = readEnergy(filename + "new")
				na = na + 1.0
				DeltaE_avg = (DeltaE_avg * (na - 1.0) + DeltaE) / na
	t = frac * t
print("Time Taken: " + str(timeTotal) + "\n")
print("Best Structure: \n")
print(structBest)
print("Best Cost: \n")
print(costBest)

#
