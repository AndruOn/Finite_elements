import numpy as np 
import matplotlib.pyplot as plt


data_14000 = np.genfromtxt("../data/mesh14608.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])

data_400 = np.genfromtxt("../data/mesh400.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_800 = np.genfromtxt("../data/mesh838.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_1600 = np.genfromtxt("../data/mesh1667.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])
data_4000 = np.genfromtxt("../data/mesh4424.csv", delimiter=",", names=["cputime","time","couple","theta","omega","coil","offset"])



time = data_14000["time"]

omega_14000 = data_14000["omega"]
omega_400 = data_400["omega"]
omega_800 = data_800["omega"]
omega_1600 = data_1600["omega"]
omega_4000 = data_4000["omega"]


plt.plot(time,omega_14000,label="mesh14608")
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_400,label="mesh400")
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_800,label="mesh838")
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_1600,label="mesh1667")
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

plt.plot(time,omega_4000,label="mesh4424")
plt.ylabel('omega (rad/s)')
plt.xlabel('time (s)')

print("cpu_time" + "mesh400" + "cpu_time:" + str(round(sum(data_400["cputime"])/301,6)))
print("cpu_time" + "mesh838" + "cpu_time:" + str(round(sum(data_800["cputime"])/301,6)))
print("cpu_time" + "mesh1667" + "cpu_time:" + str(round(sum(data_1600["cputime"])/301,6)))
print("cpu_time" + "mesh4424" + "cpu_time:" + str(round(sum(data_4000["cputime"])/301,6)))
print("cpu_time" + "mesh14608" + "cpu_time:" + str(round(sum(data_14000["cputime"])/301,6)))



plt.title('Evolution de omega pour chacun des mesh')
plt.legend()

plt.savefig("../images/omega_mesh.png")
plt.show()


#plt.loglog(h, error)