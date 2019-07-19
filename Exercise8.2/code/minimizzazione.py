import os
import numpy as np
import string

os.system('make clean')
os.system('make')

mu_min=0.7
mu_max=0.9
sigma_min=0.5
sigma_max=0.7

for mu in np.arange(mu_min, mu_max, 0.05):
	for sigma in np.arange(sigma_min, sigma_max, 0.05):
		print ("{:1.4}".format(mu),"	{:1.4}".format(sigma))
		os.system('./main.exe'+' '+format(mu)+' '+format(sigma))


mu,sigma,ene= np.loadtxt("Energie.dat", usecols=(0,1,2), delimiter='	', unpack='true')
ind = np.argmin(ene)
print("\n\nBest set of parameters:\nmu=",mu[ind]," sigma=",sigma[ind],"\npotential energy=",ene[ind],"\n\n")


os.system('make clean')
os.system('make isto')
os.system('./main_isto.exe'+' '+format(mu[ind])+' '+format(sigma[ind]))
