import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = np.loadtxt('points', skiprows=1)
theta = data[:, 0]
phi = data[:, 1]
r = theta*0+1

x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c='r', marker='o')
plt.show()
