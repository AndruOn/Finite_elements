import numpy as np 
import matplotlib.pyplot as plt


dataA = np.genfromtxt("../data/mesh4424A.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
dataB = np.genfromtxt("../data/mesh4424B.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
dataC = np.genfromtxt("../data/mesh4424C.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
dataALL = np.genfromtxt("../data/mesh4424_0,23.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])



time = dataA["time"]
omegaA = dataA["omega"]
torqueA = dataA["couple"]
theta = dataA["theta"]

torqueB = dataB["couple"]
torqueC = dataC["couple"]
torqueALL = dataALL["couple"]

plt.plot(theta,torqueA,label='Paire de bobines A')
plt.ylabel('Couple (Nm)')
plt.xlabel('Theta (rad)')

plt.plot(theta,torqueB,label='Paire de bobines B')
plt.ylabel('Couple (Nm)')
plt.xlabel('Theta (rad)')

plt.plot(theta,torqueC,label='Paire de bobines C ')
plt.ylabel('Couple (Nm)')
plt.xlabel('Theta (rad)')

plt.plot(theta,torqueALL,label='Toutes les bobines (moyenne du couple='+ str(round(sum(torqueALL/300),5)) +')',)
plt.ylabel('Couple (Nm)')
plt.xlabel('Theta (rad)')


plt.title('Couple généré par chaque bobine en fonction de theta')
plt.legend()

plt.savefig("../images/coupleCoils.png")
plt.show()


#plt.loglog(h, error)