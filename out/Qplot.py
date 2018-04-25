import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = np.loadtxt('Q.out', skiprows=1)
cos = data[:, 0]
Q1 = data[:, 1]
Q2 = data[:, 2]
Q3 = data[:, 3]

Q1mx = max(np.absolute(Q1))
Q2mx = max(np.absolute(Q2))
print Q1mx/Q2mx

fig = plt.figure()

plt.plot(cos, Q1, 'k-', label=r'$Q_{e_1}$')
plt.plot(cos, Q2, 'g:', label=r'$Q_{e_2}$')
#plt.plot(cos, Q3, 'y--', label='$Q_{e_3}$')
plt.xlabel(r'$\cos \theta$', fontsize = 16)
plt.ylabel(r'$Q_{e_i}$', fontsize = 16)

ax = fig.add_subplot(111)
ax.set_xlim((-1,1))
plt.legend()
plt.show()
