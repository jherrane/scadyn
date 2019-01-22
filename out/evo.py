import math, codecs, os, errno, sys, getopt, numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, patches as mpatches, rc
from itertools import islice
from io import StringIO
from mpl_toolkits.mplot3d import proj3d

class cd:
   def __init__(self, newPath):
      self.make_sure_path_exists(newPath)
      self.newPath = os.path.expanduser(newPath)
   def __enter__(self):
      self.savedPath = os.getcwd()
      os.chdir(self.newPath)
   def __exit__(self, etype, value, traceback):
      os.chdir(self.savedPath)
   def make_sure_path_exists(self,path):
      try:
         os.makedirs(path)
      except OSError as exception:
         if exception.errno != errno.EEXIST:
            raise

file = "evo"

def plot_orbit(orbit1, orbit2, orbit3, J, title):
   MAP = 'winter'
   COLOR1 = 'cyan'
   COLOR2 = 'magenta'
   COLOR3 = 'blue'
   COLOR4 = 'black'
   x1, x2, x3 = zip(*orbit1)
   y1, y2, y3 = zip(*orbit2)
   z1, z2, z3 = zip(*orbit3)
   J1, J2, J3 = zip(*J)
   N = len(x1)
   fig = plt.figure(figsize=(16,14))
   ax = fig.add_subplot(111, projection='3d')
  
   ax.set_aspect('equal')
   step = N/50 + 1
   ax.plot(x1, x2, x3, color=COLOR1)
   ax.plot(y1, y2, y3, color=COLOR2)
   ax.plot(z1, z2, z3, color=COLOR3)
   ax.plot(J1, J2, J3, color=COLOR4)

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
   J_patch  = mpatches.Patch(color=COLOR4,label=r'$\hat{\mathbf{J}}$')
   plt.legend(handles = [a1_patch, a2_patch, a3_patch, J_patch],prop={'size':44},loc='lower right')
   plt.savefig(file + '.png',bbox_inches='tight',pad_inches=0.05)

def plot_R(R_list,J,Q,title):
   v1 = []
   v2 = []
   v3 = []
   for R1 in R_list: 
      R = np.asarray(R1).reshape((3,3))
      RQ = np.dot(R,np.transpose(Q))
      v1.append(RQ[0,:]*0.33)
      v2.append(RQ[1,:]*0.66)
      v3.append(RQ[2,:])
   plot_orbit(v1,v2,v3,J,title)

if __name__ == "__main__":
   inputfile = 'log'
   skip = 22
   last = 0 # Draw only last n lines, 0 if all
   first = 0 # Draw only first n lines, 0 if all

   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   font = {'weight' : 'bold', 'size' : 48}
   mpl.rc('font',**font)
   
   argv = sys.argv[1:]
   try:
      opts, args = getopt.getopt(argv,"hi:e:f:")
   except getopt.GetoptError:
      print 'python evo.py -h'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'evo.py -i <inputfile>'
         print '       -f n'
         print '       -e n'
         sys.exit()
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-e"):
         last = int(arg)
      elif opt in ("-f"):
         first = int(arg)
   
   splt = inputfile.split("log")
   file = file + splt[-1]

   log = open(inputfile,'r')
   num_lines = sum(1 for line in log)
   log.seek(0)
   
   inputf = codecs.open(inputfile, encoding='utf-8').read()
   inputf = inputf.replace('|','')
   lines = np.loadtxt(StringIO(inputf), skiprows=skip)
   
   if first != 0:
      lines = lines[0:first]
      num_lines = 22+first
   if last != 0:
      skip = num_lines-last
      
   Q = np.genfromtxt(islice(log,17,20))
   log.seek(0)
   I = np.genfromtxt(islice(log,15,16))
   I = np.diag(I)
   log.close()
   title = 'Rotations during steps ' + str(skip-22+1) + '-' +str(num_lines-22) 

   w = lines[:,5:8]
   J = 0*w
   R = lines[:,17:27]
      
   for i in range(0,w.shape[0]):
      J[i,:] = np.matmul(I,w[i,:])
      J[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),J[i,:]  ))
      J[i,:] = 1.1*J[i,:]/np.sqrt(np.sum(np.power(J[i,:],2)))

   with cd('figs'):
		plot_R(R,J,Q,title)
