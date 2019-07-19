import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import string
fig=plt.figure(figsize=(10, 5))

i, x, y = np.loadtxt("BestTrip.dat", usecols=(0,1,2), delimiter='	', unpack='true')

title = str(np.loadtxt("BestLength.dat"))
plt.plot(x,y,'.-b')
plt.title('Length: '+title)
plt.grid()

plt.tight_layout()
#plt.savefig("./pictures/BestTrips.png")

plt.show()
