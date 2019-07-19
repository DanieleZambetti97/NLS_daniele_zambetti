import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

os.system("rm risultati/h_0.0/*")
os.system("rm risultati/*")

ene = open("risultati/h_0.0/temp_ene.dat", "a")
heat = open("risultati/h_0.0/temp_heat.dat", "a")
susce = open("risultati/h_0.0/temp_susce.dat", "a")
magn = open("risultati/h_0.0/temp_magn.dat", "a")

metro= [0,1]
temp = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]
j=1.0
h=0.0
nblk=20
nstep=1000
i=0;

for m in metro:
	for t in temp:
		os.system("./Monte_Carlo_ISING_1D.exe 0 "+str(m)+" "+str(t)+" "+str(j)+" "+str(h)+" "+str(nblk)+" "+str(nstep)+" > eliminami.dat")
		print("Restart=0		Metro="+str(m)+"		temp="+str(t)+"		j="+str(j)+"		h="+str(h)+"	"+str(nblk)+" "+str(nstep))
		inutile, y, yerror = np.loadtxt("./risultati/DataBlocking_Ene.dat", usecols=(0,1,2), delimiter='	', unpack='true', skiprows=int( nblk-1 + nblk*i ), max_rows=1 )
		ene.write(str(t)+"	"+str(y)+"	"+str(yerror)+"\n")

		inutile, y, yerror = np.loadtxt("./risultati/DataBlocking_Heat.dat", usecols=(0,1,2), delimiter='	', unpack='true', skiprows=int( nblk-1 + nblk*i ), max_rows=1 )
		heat.write(str(t)+"	"+str(y)+"	"+str(yerror)+"\n")

		inutile, y, yerror = np.loadtxt("./risultati/DataBlocking_Susce.dat", usecols=(0,1,2), delimiter='	', unpack='true', skiprows=int( nblk-1 + nblk*i ), max_rows=1 )
		susce.write(str(t)+"	"+str(y)+"	"+str(yerror)+"\n")

		inutile, y, yerror = np.loadtxt("./risultati/DataBlocking_Magn.dat", usecols=(0,1,2), delimiter='	', unpack='true', skiprows=int( nblk-1 + nblk*i ), max_rows=1 )
		magn.write(str(t)+"	"+str(y)+"	"+str(yerror)+"\n")
		i=i+1

ene.close()
heat.close()
susce.close()
magn.close()
