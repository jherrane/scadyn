import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = np.loadtxt('points', skiprows=1)
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c='r', marker='o')
plt.show()
