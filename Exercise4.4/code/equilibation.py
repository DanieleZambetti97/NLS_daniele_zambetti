import os

os.system("make clean")
os.system("make 2>eliminami.dat")

desire_temp=float(open('input.dat', 'r').readline())
print("\n \nDisidered temperature = ",desire_temp," +- 0.01")
j=0
final_temp=0

os.system("./MolDyn_NVE.x 0 >eliminami.dat")

with open("./risultati/ave_temp.dat", "r") as file:
	lines = file.read().splitlines()
i,final_temp,error=map(float, lines[-1].split(' '))
print("Final temperature of the starting try = ",final_temp)


while(final_temp<desire_temp-0.01 or final_temp>desire_temp+0.01):
    j += 1
    os.system("./MolDyn_NVE.x 1 >eliminami.dat")
    with open("./risultati/ave_temp.dat", "r") as file:
        lines = file.read().splitlines()
    i,final_temp,error=map(float, lines[-1].split(' '))
    print("Final temperature of the ",j,"-th try = ",final_temp)
    if(j==50):
       break
