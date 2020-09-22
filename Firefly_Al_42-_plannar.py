import threading
import numpy as np
import time
import os
import math
import os.path
charge=-2
spin=1
molname="Al"
numtimes=100
count=0
norm_energies=0
best_energy=0
num_calls=0
time_in_secs=time.time()

#Checks if the termination is normal(returns 1) or error(returns 2) or incomplete(returns 0)
def normal_term(file_to_open):
	file=open(file_to_open,'r+')
	X=file.readlines()
	for i in range(len(X)):
		Y=[]
		Y=X[i].split(" ")
		if len(Y)>=3 and Y[2]=="termination":
			if Y[1]=='Normal':
				file.close()
				return 1
			elif Y[1]=="Error":
				file.close()
				return 2
	
	file.close()
	return 0
#Returns the energy. Also appends the different values.
def energy(name,X,partnum):
		file_to_open=name[:-3]+"log"
		if(normal_term(file_to_open)==2):
			return X
		file=open(file_to_open,'r+')
		Y=file.readlines()
		for i in range(len(Y)):
			Z=[]
			Z=Y[i].split(" ")
			for j in range(len(Z)):
				if "E(" in Z[j]:
					X=Z[j+3]
					file.close()
					lock.acquire()
					valuesnow[partnum]=float(X)
					global norm_energies
					norm_energies +=1
					lock.release()
					return X
		file.close()


#Creates the required gjf format with dynamic threads
def gjfconvert(X,iternum,partnum,threadnums,whereto):
	molname="Al"
	name="str"+str(partnum)+"__"+str(whereto)+"_"+str(iternum)+".gjf"
	file=open(name,"w+")    
	file.write("%nprocshared="+str(int(30/threadnums))+"\n")
	file.write("%mem=2GB\n")
	file.write("%nosave\n")
	file.write("# b3lyp/6-311+g*\n\n")#CCSD(T)/6-311g*
	file.write("Al struct\n\n")
	file.write(str(charge)+" "+str(spin)+"\n")
	for i in range(int(len(X)/2)):
		string=" "+molname+"\t\t\t\t\t"+str(round(X[i*2],8))+"\t"+str(round(X[i*2+1],8))+"\t0.00000000\n"
		file.write(string)
		
	file.write("\n\n\n\n")
	file.close()
	return name

#Checks if the files are created or not and removes the gjf files
def function(position,iternum,partnum,threadnums,whereto):
		name=gjfconvert(position,iternum,partnum,threadnums,whereto)
		os.system("g09 "+name+" &")
		#time.sleep(3)
		while (os.path.isfile(name[:-3]+"log")) is False:
		   pass
		while normal_term(name[:-3]+"log")==0: 
			pass
		energy(name,valuesnow[partnum],partnum)
		os.remove(name)
		os.remove(name[:-3]+"log")

#hyperparameters
'''
num_feat=int(input("Enter the number of particles each cluster has:"))
num_units=int(input("Enter the number of particle clusters:"))
c1=float(input("Enter positional dependence:"))
c2=float(input("Enter global dependence:"))
n=int(input("Enter the number of times it will run:"))
'''
num_feat=4
num_units=4
c1=0.2
c2=0.9/8
a=0.1
p=0.97
n=1000
num_feat=2*num_feat
#initialisation
position=(np.random.random((num_units,num_feat))-0.5)*6
#Checks if any Present values are there. If present,it uses them
if(os.path.isfile("PresentValues.txt")):
	position=position.reshape((-1,))
	file=open("PresentValues.txt",'r+')
	count=0
	X=file.readlines()
	for i in range(len(X)):
		X[i]=X[i][:-1]
		Y=X[i].split("\t")
		for j in range(len(Y)):
			try:
				position[count]=float(Y[j])
				count +=1
			except:
				pass
	position=position.reshape((num_units,num_feat))
	file.close()
best_individual=np.copy(position)
valuesnow=np.zeros((num_units,1),dtype=float)
#Creates the Threads

print("Iteration starts at "+time.ctime(time.time()))
t=[]
lock=threading.Lock()
for i in range(num_units):
	t.append(threading.Thread(target=function,args=(position[i],0,i,num_units,i)))
	t[i].start()
for i in range(num_units):
	t[i].join()
best_pos=np.copy(position[np.argmin(valuesnow),:].reshape(1,num_feat))
temp=np.zeros((num_units,1),dtype=int)
indibestnow=np.copy(valuesnow)
best_energy=min(valuesnow)
#Creates a Folder to save the files
foldertosave="valueAl_"
u=0
while(os.path.exists(foldertosave+str(u))):
	u =u+1
os.mkdir(foldertosave+str(u))
foldertosave="valueAl_"+str(u)
#Saves the values
k=0
init=open(foldertosave+"/initialvalues.txt","w")
for i in range(num_units):
	init.write("\n\n"+str(i)+"\n")
	for j in range(int(num_feat/2)):
		init.write(" "+molname+"\t\t\t\t\t"+str(position[i][j*2])+"\t"+str(position[i][j*2+1])+"\t0.00000000\t\n")
init.close()
count=0
#Creates the main part of Firefly
num_calls +=num_units
k=0
print("\n............\n.............\nThe best pos of iteration "+str(k) +"is:\n")
for i in range(int(num_feat/2)):
	string=str(round(best_pos[0][i*2],8))+"\t"+str(round(best_pos[0][i*2+1],8))+"\t0.00000000\n"
	print(string)
print("The best energy is"+str(min(valuesnow)))
print("\nThe time from start is: ")
print(time.time()-time_in_secs)

for k in range(1,n):
	p=p*0.999
	for i in range(num_units):
		toedit=[]
		for j in range(num_units):
			if(valuesnow[i]<valuesnow[j]):
				dist=np.dot((position[j]-position[i]).T,(position[j]-position[i]))
				position[j]= position[j]+c1*(position[i]-position[j])*((math.e)**(-c2*dist))+0.1*(np.random.random((1,num_feat))-0.5)*p

				toedit.append(j)
		t=[]
		num_calls +=len(toedit)
		for j in range(len(toedit)):
			t.append(threading.Thread(target=function,args=(position[toedit[j]],k,toedit[j],len(toedit),i)))
			t[j].start()
		for j in range(len(t)):
			t[j].join()

	best_pos=np.copy(best_individual[np.argmin(indibestnow),:].reshape(1,num_feat))

	print("\n............\n.............\nThe best pos of iteration "+str(k) +" is:\n")
	for i in range(int(num_feat/2)):
		string=" Al\t\t\t\t"+str(round(best_pos[0][i*2],8))+"\t"+str(round(best_pos[0][i*2+1],8))+"\t0.00000000\n"
		print(string)
	print("The best energy is"+str(min(valuesnow)))
	print("\nThe time from start is :")

	print(time.time()-time_in_secs)

	if(float(min(valuesnow))<(float(best_energy)-0.00001)):
		count =0
		best_energy=min(valuesnow)

	else:
		count +=1
		
	if(count>=numtimes):
		print("No improvement. Program closed.")
		break


print("--------------\n------------------\nThe number of calls made is "+str(num_calls)+"\nThe number of normal termination is "+str(norm_energies))
os.system("mv output.txt "+foldertosave+"/output.txt")


time.sleep(3)
#uncomment to run again
os.system("nohup python -u FireflyAl.py &> output.txt &")
