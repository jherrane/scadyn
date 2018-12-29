import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = np.loadtxt('q', skiprows=1)
q = data[:, 0]
t = data[:, 2]/3600/24/365

fig = plt.figure()

plt.plot(t, q, 'k-', label=r'$q(t)$')
plt.xlabel(r'$t (years)$', fontsize = 16)
plt.ylabel(r'$q$', fontsize = 16)

ax = fig.add_subplot(111)
ax.set_ylim(bottom=1)
ax.set_xlim(left=0)
plt.legend()
plt.show()
