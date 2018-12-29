# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from cd import *
from arrow import *
import sys, getopt
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math
from itertools import islice
from io import StringIO
import codecs
from matplotlib import rc

fileQ = "Qt"

def plot_orbit(orbit1, orbit2, orbit3,w, k,title):
   MAP = 'winter'
   COLOR1 = 'cyan'
   COLOR2 = 'magenta'
   COLOR3 = 'blue'
   COLOR4 = 'black'
   x1, x2, x3 = zip(*orbit1)
   y1, y2, y3 = zip(*orbit2)
   z1, z2, z3 = zip(*orbit3)
   w1, w2, w3 = zip(*w)
   N = len(x1)
   font = {'weight' : 'bold', 'size' : 40}
   fig = plt.figure(figsize=(12,12))
   ax = fig.add_subplot(111, projection='3d')
   
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   ax.set_aspect('equal')
   step = N/50 + 1
   ax.plot(x1, x2, x3, color=COLOR1)
   ax.plot(y1, y2, y3, color=COLOR2)
   ax.plot(z1, z2, z3, color=COLOR3)
   ax.plot(w1, w2, w3, color=COLOR4)
   a = Arrow3D([0,k[0]],[0,k[1]],[0,k[2]], mutation_scale=20, lw=3, arrowstyle="-|>", color="k")
   ax.add_artist(a)

   ax.set_xlim(-1,1)
   ax.set_ylim(-1,1)
   ax.set_zlim(-1,1)
   ax.set_title(title, fontweight='bold', fontsize=50)
   ax.xaxis.set_ticklabels([])
   ax.yaxis.set_ticklabels([])
   ax.zaxis.set_ticklabels([])

   a1_patch = mpatches.Patch(color=COLOR1,label=r'$\hat{\mathbf{a}}_{1}$')
   a2_patch = mpatches.Patch(color=COLOR2,label=r'$\hat{\mathbf{a}}_{2}$')
   a3_patch = mpatches.Patch(color=COLOR3,label=r'$\hat{\mathbf{a}}_{3}$')
   w_patch = mpatches.Patch(color=COLOR4,label=r'$\hat{\mathbf{J}}$')
   plt.legend(handles = [a1_patch, a2_patch, a3_patch, w_patch],prop={'size':44},loc='lower right')
   plt.savefig(fileQ + '.png',bbox_inches='tight',pad_inches=0.05)

def plot_R(R_list,I,w,Q,k,title):
   I3 = I[2]
   I2 = I[1]
   I1 = I[0]
   v1 = []
   v2 = []
   v3 = []
   for R1 in R_list: 
      R = np.asarray(R1).reshape((3,3))
      RQ = np.dot(R,np.transpose(Q))
      v1.append(RQ[0,:]*0.33)
      v2.append(RQ[1,:]*0.66)
      v3.append(RQ[2,:])
   plot_orbit(v1,v2,v3,w,k,title)

if __name__ == "__main__":
   inputfile = 'log'
   pth = 'figs'
   skip = 22
   last = 0 # Draw only last n lines, 0 if all
   first = 0 # Draw only first n lines, 0 if all
   #plt.xkcd()
   #plt.rc('font',family='Times New Roman')
   
   argv = sys.argv[1:]
   try:
      opts, args = getopt.getopt(argv,"hi:p:s:e:f:")
   except getopt.GetoptError:
      print 'python evo.py -h'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'evo.py -i <inputfile>'
         print '       -p <path>'
         print '       -s n'
         print '       -e n'
         sys.exit()
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-p"):
         pth = arg
      elif opt in ("-s"):
         skip = int(arg)
      elif opt in ("-e"):
         last = int(arg)
      elif opt in ("-f"):
         first = int(arg)
   
   splt = inputfile.split("log")
   fileQ = fileQ + splt[-1]

   log = open(inputfile,'r')
   num_lines = sum(1 for line in log)
   log.seek(0)
   if last != 0:
      skip = num_lines-last
   inputf = codecs.open(inputfile, encoding='utf-8').read()
   inputf = inputf.replace('|','')
   
   lines = np.loadtxt(StringIO(inputf), skiprows=skip)
   if first != 0:
      lines = lines[0:first]
      num_lines = 22+first
   
   string = log.readlines()[1]
   k = [float(s) for s in string.split()[2:5]]
   log.seek(0)
   Q = np.genfromtxt(islice(log,17,20))
   log.seek(0)
   I = np.genfromtxt(islice(log,15,16))
   I = np.diag(I)
   log.seek(0)
   a = np.genfromtxt(islice(log,10,11))
   a = a[3]*1e-9
   log.close()
   title = 'Evolution of rotation between\nsteps ' + str(skip-22+1) + '-' +str(num_lines-22) 

   t = lines[:,1]
   x = lines[:,2:5]/a
   w = lines[:,5:8]
   v = lines[:,8:11]
   J = 0*w
   N = lines[:,11:14]
   F = lines[:,14:17]
   R = lines[:,17:27]
   q = 0*t
      
   for i in range(0,w.shape[0]):
      J[i,:] = np.matmul(I,w[i,:])
      w[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),w[i,:]  ))
      J[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),J[i,:]  ))
      q[i] = I[2,2]*np.matmul(J[i,:],w[i,:])/np.matmul(J[i,:],J[i,:])
      w[i,:] = 1.1*w[i,:]/np.sqrt(np.sum(np.power(w[i,:],2)))
      J[i,:] = 1.1*J[i,:]/np.sqrt(np.sum(np.power(J[i,:],2)))
   
   I = [I[0,0], I[1,1], I[2,2]]
   I = I/np.linalg.norm(I)
   with cd(pth):
		plot_R(R,I,J,Q,k,title)
