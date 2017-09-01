import matplotlib.pyplot as plt
import numpy as np
import math

log = open('Q.out','r')

cos = []
Q = []

lines = log.readlines()[1:]

for line in lines:
 data = map(float,line.split())
 cos.append(data[0])
 Q.append(data[1:4])
log.close()

Q1, Q2, Q3 = zip(*Q)

Q1mx = max(np.absolute(Q1))
Q2mx = max(np.absolute(Q2))
print Q1mx/Q2mx
Q3mx = max(np.absolute(Q3))
#Q1 = tuple([x/Q1mx for x in Q1])
#Q2 = tuple([x/Q2mx for x in Q2])
#Q3 = tuple([x/Q3mx for x in Q3])

fig = plt.figure()

Qone, = plt.plot(cos, Q1, 'k-', label='Q_e1')
Qtwo, = plt.plot(cos, Q2, 'g:', label='Q_e2')
#Qtri, = plt.plot(cos, Q3, 'y--', label='Q_e3')

ax = fig.add_subplot(111)
ax.set_xlim((-1,1))
plt.legend(handles=[Qone,Qtwo])
#plt.plot(cos, Q1, 'k-', cos, Q2, 'g--')
#plt.plot(cos, Q3, 'b--')
plt.show()
