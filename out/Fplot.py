import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = np.loadtxt('F', skiprows=1)
xi = data[:, 0]
psi = data[:, 1]
F = data[:, 2]

psi_unique = np.unique(psi)

fig = plt.figure()

for p in psi_unique: 
 angle = str(round(np.arccos(p) * 180 / np.pi)).split('.')[0]
 label = r'$\psi = ' + angle + '^{\circ}$'
 
 plt.plot(xi[psi == p], F[psi == p], '-', label = label)
 plt.xlabel(r'\cos \xi')
 plt.ylabel(r'F', fontsize = 16)

ax = fig.add_subplot(111)
ax.set_xlim((-1, 1))
plt.legend()
plt.show()
