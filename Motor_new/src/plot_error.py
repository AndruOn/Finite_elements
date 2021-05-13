import numpy as np 
import matplotlib.pyplot as plt
import math


data_14000 = np.genfromtxt("../data/mesh14608.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])

data_400 = np.genfromtxt("../data/mesh400.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_800 = np.genfromtxt("../data/mesh838.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_1600 = np.genfromtxt("../data/mesh1667.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_4000 = np.genfromtxt("../data/mesh4424.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])



time = data_14000["time"]


list_data = [data_400,data_800,data_1600,data_4000,data_14000]

noeuds = [400,838,1667,4424,14608]
noeuds_log = [math.log(400, math.e),math.log(838, math.e),math.log(1667, math.e),math.log(4424, math.e),math.log(14608, math.e)]


cputime = []
error = []
i = 0
for data in list_data:
    cputime.append( math.log(sum(data["cputime"])/301,math.e) )
    delta = abs(data["couple"][0] - data_14000["couple"][0])
    print("delta",delta)
    if delta != 0.0 :
        error.append( math.log(delta, math.e) )
    else:
        error.append( 0 )


plt.scatter(noeuds,cputime,label="log(moyenne de temps cpu) (s)")
plt.scatter(noeuds[:-1],error[:-1],label="log(erreur) (-)")

plt.xlabel('nombres de noeuds (-)')
plt.xscale('log') 





plt.title('Erreur et temps cpu pour chacun des mesh')
plt.legend()

plt.savefig("../images/error_cputime.png")
plt.show()


#plt.loglog(h, error)