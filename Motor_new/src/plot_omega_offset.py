import numpy as np 
import matplotlib.pyplot as plt


data_17 = np.genfromtxt("../data/mesh4424_0,17.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_2 = np.genfromtxt("../data/mesh4424_0,2.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_23 = np.genfromtxt("../data/mesh4424_0,23.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])



time = data_17["time"]

omega_17 = data_17["omega"]
omega_2 = data_2["omega"]
omega_23 = data_23["omega"]


plt.plot(time,omega_17,label='offset = 0.17')
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_2,label='offset = 0.2')
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_23,label='offset = 0.23')
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')


plt.title('Evolution de omega en fonction du offset')
plt.legend()

plt.savefig("../images/omega_offset.png")
plt.show()


#plt.loglog(h, error)