import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, getopt

Qfile = 'Q'
inputfile = ''
pth = ''

argv = sys.argv[1:]
try:
   opts, args = getopt.getopt(argv,"hi:p:")
except getopt.GetoptError:
   print 'python Fplot.py -h'
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h' or opt == '-help':
      print 'Fplot.py -i <inputfile>'
      print 'Fplot.py -i <inputfile>'
      sys.exit()
   elif opt in ("-i"):
      inputfile = arg
   elif opt in ("-p"):
      pth = arg+'/'

fig = plt.figure(figsize=(12,8))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'weight' : 'bold', 'size' : 22}
mpl.rc('font',**font)

Qfile = pth+Qfile+inputfile
print Qfile
data = np.loadtxt(Qfile, skiprows=1)

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

ax = plt.gca()
ax.set_xlim((-1,1))
plt.legend()
plt.savefig(Qfile+'.png')
