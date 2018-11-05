import matplotlib.pyplot as plt
import numpy as np
import sys, getopt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Ffile = 'F'
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

Ffile = pth+Ffile+inputfile

data = np.loadtxt(Ffile, skiprows=1)
xi = data[:, 0]
psi = data[:, 1]
F = data[:, 2]
H = data[:,3]

psi_unique = np.unique(psi)

fig = plt.figure()

for p in psi_unique: 
 angle = str(round(np.arccos(p) * 180 / np.pi)).split('.')[0]
 label = r'$\psi = ' + angle + '^{\circ}$'
 
 plt.plot(xi[psi == p], F[psi == p], '-', label = 'F('+label+')')
 plt.plot(xi[psi == p], H[psi == p], label = 'H('+label+')')
 plt.xlabel(r'$\cos \xi$')
 plt.ylabel(r'$F, H$', fontsize = 16)

ax = plt.gca()
ax.set_xlim((-1, 1))
plt.legend()
plt.savefig(Ffile+'.png')
