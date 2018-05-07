# Module to plot the results of T-VIEDYN integration (WIP)
#from cd import cd, make_sure_path_exists
import StringIO
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import errno
from itertools import islice
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
   def __init__(self, xs, ys, zs, *args, **kwargs):
      FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
      self._verts3d = xs, ys, zs

   def draw(self, renderer):
      xs3d, ys3d, zs3d = self._verts3d
      xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
      self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
      FancyArrowPatch.draw(self, renderer)

class cd:
   """Context manager for changing the current working directory"""
   def __init__(self, newPath):
      self.newPath = os.path.expanduser(newPath)

   def __enter__(self):
      self.savedPath = os.getcwd()
      os.chdir(self.newPath)

   def __exit__(self, etype, value, traceback):
      os.chdir(self.savedPath)

def make_sure_path_exists(path):
   try:
      os.makedirs(path)
   except OSError as exception:
      if exception.errno != errno.EEXIST:
         raise

def Round_To_n(x, n):
   return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def plot_orbit(orbit, title, figname):
   MAP = 'winter'
   COLOR = 'blue'
   x1, x2, x3 = zip(*orbit)
   max_range = np.array([max(x1)-min(x1), max(x2)-min(x2), max(x3)-min(x3)]).max() / 2.0

   mid_x = (max(x1)+min(x1)) * 0.5
   mid_y = (max(x2)+min(x2)) * 0.5
   mid_z = (max(x3)+min(x3)) * 0.5
   NPOINTS = len(x1)
   fig = plt.figure(figsize=(12,12))
   ax = fig.add_subplot(111, projection='3d')
   ax.set_aspect('equal')
   ax.plot(x1,x2,x3)
   ax.plot(x1,x2,(mid_z - max_range)*np.ones(len(x3)),alpha=0.1,color='k')
   ax.plot(x1,np.ones(len(x2))*(mid_y + max_range),x3,alpha=0.1,color='k')
   ax.plot((mid_x - max_range)*np.ones(len(x1)),x2,x3,alpha=0.1,color='k')

   ax.set_xlim(mid_x - max_range, mid_x + max_range)
   ax.set_ylim(mid_y - max_range, mid_y + max_range)
   ax.set_zlim(mid_z - max_range, mid_z + max_range)

   ax.set_title(title,fontweight='bold')
   plt.locator_params(nbins=4)
   plt.savefig(figname)

def plot_fig(x, y, markevery, title, xlabel, ylabel, figname):
   font = {'weight' : 'bold',
			   'size'   : 18}
   fig = plt.figure(figsize=(12,12))
   plt.autoscale(enable=True, axis='x', tight=True)

   mpl.rc('font',**font)
   y1,y2,y3 = zip(*y)
   plt.plot(x,y1,label=r'$'+ylabel+'_{1}$',lw=3.0,markevery=markevery)
   plt.plot(x,y2,label=r'$'+ylabel+'_{2}$',lw=3.0,markevery=markevery)
   plt.plot(x,y3,label=r'$'+ylabel+'_{3}$',lw=3.0,markevery=markevery)
   plt.legend()
   plt.xlabel(r'$'+xlabel+'$')
   plt.ylabel(r'$'+ylabel+'$')
   plt.title(title,fontweight='bold')
   plt.savefig(figname)
   plt.close(fig)

def runningMean(x, N):
   cumsum = np.cumsum(np.insert(x, 0, 0)) 
   return (cumsum[N:] - cumsum[:-N]) / N

def plot_mav(x,y,n,title,xlabel,ylabel,figname):
   font = {'weight' : 'bold',
			   'size'   : 18}
   fig = plt.figure(figsize=(12,12))
   plt.autoscale(enable=True, axis='x', tight=True)
   miny = min(np.ndarray.flatten(np.asarray(y)))
   maxy = max(np.ndarray.flatten(np.asarray(y)))
   plt.ylim(1.3*miny,1.3*maxy)
   mpl.rc('font',**font)
   ym1,ym2,ym3 = zip(*y)
   y1 = runningMean(ym1,n)
   y2 = runningMean(ym2,n)
   y3 = runningMean(ym3,n)
   plt.plot(x[0:-n+1],y1,label=r'$'+ylabel+'_{1}$',lw=3.0)
   plt.plot(x[0:-n+1],y2,label=r'$'+ylabel+'_{2}$',lw=3.0)
   plt.plot(x[0:-n+1],y3,label=r'$'+ylabel+'_{3}$',lw=3.0)
   plt.legend()
   plt.xlabel(r'$'+xlabel+'$')
   plt.ylabel(r'$'+ylabel+'$')
   plt.title(title,fontweight='bold')
   plt.savefig(figname)
   plt.close(fig)

def plot_R(R_list,Q,k,title,figname):
   font = {'weight' : 'bold', 'size' : 18}
   fig = plt.figure(figsize=(12,12))
   # Set up test vector
   #v_test = (Q[:,0]+Q[:,1]+Q[:,2])/np.sqrt(3)
   v1 = []
   v2 = []
   v3 = []
   for R1 in R_list: 
    R = np.asarray(R1).reshape((3,3))
    RQ = np.dot(R,np.transpose(Q))
    v1.append(RQ[0,:])
    v2.append(RQ[1,:])
    v3.append(RQ[2,:])
   plot_orbit2(v1,k,title+' a_1' ,'a_1'+figname)
   plot_orbit2(v2,k,title+' a_2' ,'a_2'+figname)
   plot_orbit2(v3,k,title+' a_3' ,'a_3'+figname)

if __name__ == "__main__":
   #This block will be evaluated if this script is called as the main routine
   #and will be ignored if this file is imported from another script.

   pth = 'figs'
   make_sure_path_exists(pth)
   
   #plt.xkcd()
   
   log = open('log','r')
   lines = np.loadtxt('log', skiprows=22)
      
   magtol = 1e-3

   string = log.readlines()[3]
   Nmax = [int(s) for s in string.split() if s.isdigit()]
   Nmax = Nmax[0]
   log.seek(0)
   string = log.readlines()[1]
   k = [float(s) for s in string.split()[2:5]]
   log.seek(0)
   Q = np.genfromtxt(islice(log,17,20))
   log.close()
   
   markevery = round(Nmax/1000)
   mav = int(round(Nmax/100))
   if mav < 2: mav = 2
   if markevery < 1: markevery = 1

   x = lines[:,1:4]
   v = lines[:,4:7]
   w = lines[:,7:10]
   J = lines[:,10:13]
   N = lines[:,13:16]
   F = lines[:,16:19]
   t = lines[:,19]
   R = lines[:,20:30]
   Jb = np.zeros(J.shape)
   wb = np.zeros(w.shape)
   for i in np.arange(R.shape[0]):
      wlen = np.linalg.norm(w[i,:])
      RR = R[i,:].reshape(3,3)
      if wlen != 0:
         wb[i,:] = (RR.dot(Q)).T.dot(w[i,:])/wlen
         Jb[i,:] = (RR.dot(Q)).T.dot(J[i,:])/np.linalg.norm(J[i,:])        
      else:
         wb[i,:] = [0, 0, 0]
         Jb[i,:] = [0, 0, 0]
         
   with cd(pth):
      plot_fig(t,w,markevery,'Angular velocity vs. Time','t (s)','\omega','w.png')
      plot_fig(t,J,markevery,'Angular momentum vs. Time','t (s)','J (Nms)','J.png')
      plot_fig(t,F,markevery,'Force vs. Time','t (s)','F (N)','F.png')
      plot_fig(t,N,markevery,'Torque vs. Time','t (s)','N (Nm)','N.png')
      plot_orbit(w,'Spin of angular velocity','worbit.png')
      plot_orbit(wb,'Spin of angular velocity in body frame','wborbit.png')
      plot_orbit(J,'Time evolution of angular momentum','Jorbit.png')
      plot_orbit(Jb,'Angular momentum direction in body frame','Jborbit.png')
      plot_R(R,Q,k,'Time evolution of','.png')
      plot_mav(t,F,mav,'Running average of Force vs. Time','t (s)','F (Nm)','Fav.png')
      plot_mav(t,w,mav,'Running average of angular velocity vs. Time','t (s)','\omega (rad/s)','wav.png')
      plot_mav(t,N,mav,'Running average of torque vs. Time','t (s)','N (Nm)','Nav.png')

