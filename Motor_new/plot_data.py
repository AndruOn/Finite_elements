import numpy as np 
import matplotlib.pyplot as plt

filename = "build/mmesh400.txt"
X = np.loadtxt(filename)

t = X[:,1]
omega = X[:,4]
torque = X[:,2]


# t = np.arange(len(X))

fig, ax = plt.subplots(1,2)
ax[0].plot(t,omega)
ax[0].set_title("Omega")
ax[1].plot(t,torque)
ax[1].set_title("Couple")


plt.figure();
plt.plot(omega*180/np.pi,torque);



plt.show()

#plt.loglog(h, error)