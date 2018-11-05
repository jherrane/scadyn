import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fig = plt.figure(figsize=(12,8))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'weight' : 'bold', 'size' : 22}
mpl.rc('font',**font)

data = np.loadtxt('Q', skiprows=1)
cos = data[:, 0]
Q1 = data[:, 1]
Q2 = data[:, 2]
Q3 = data[:, 3]

Q1mx = max(np.absolute(Q1))
Q2mx = max(np.absolute(Q2))
print Q1mx/Q2mx

plt.plot(cos, Q1, 'k-', label=r'$Q_{e_1}$')
plt.plot(cos, Q2, 'g:', label=r'$Q_{e_2}$')
#plt.plot(cos, Q3, 'y--', label='$Q_{e_3}$')
plt.xlabel(r'$\cos \theta$')
plt.ylabel(r'$Q_{e_i}$')

ax = fig.add_subplot(111)
ax.set_xlim((-1,1))
plt.legend()
plt.show()
