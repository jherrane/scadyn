# Module to plot the results of T-VIEDYN integration (WIP)
#from cd import cd, make_sure_path_exists
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
# cm = plt.get_cmap(MAP)
# ax.set_color_cycle([cm(1.*i/(NPOINTS-1)) for i in range(NPOINTS-1)])
# for i in range(NPOINTS):
#  ax.plot(x1[i:i+2],x2[i:i+2],x3[i:i+2])
#  ax.plot(x1[i:i+2],x2[i:i+2],x3[i:i+2],alpha=float(i)/(NPOINTS-1),color=COLOR)

 ax.set_xlim(mid_x - max_range, mid_x + max_range)
 ax.set_ylim(mid_y - max_range, mid_y + max_range)
 ax.set_zlim(mid_z - max_range, mid_z + max_range)
# ax.set_xlim(-1,1)
# ax.set_ylim(-1,1)
# ax.set_zlim(-1,1)
 ax.set_title(title,fontweight='bold')
 plt.locator_params(nbins=4)
 plt.savefig(figname)
# plt.show()
# plt.close(fig)

def plot_orbit2(orbit, k, title, figname):
 MAP = 'winter'
 COLOR = 'blue'
 x1, x2, x3 = zip(*orbit)
 NPOINTS = len(x1)
 fig = plt.figure(figsize=(12,12))
 ax = fig.add_subplot(111, projection='3d')
 ax.set_aspect('equal')
 ax.plot(x1,x2,x3)
 a = Arrow3D([0,k[0]],[0,k[1]],[0,k[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
 ax.add_artist(a)
# cm = plt.get_cmap(MAP)
# ax.set_color_cycle([cm(1.*i/(NPOINTS-1)) for i in range(NPOINTS-1)])
# for i in range(NPOINTS):
#  ax.plot(x1[i:i+2],x2[i:i+2],x3[i:i+2])
#  ax.plot(x1[i:i+2],x2[i:i+2],x3[i:i+2],alpha=float(i)/(NPOINTS-1),color=COLOR)
 max_range = np.array([max(x1)-min(x1), max(x2)-min(x2), max(x3)-min(x3)]).max() / 2.0

 mid_x = (max(x1)+min(x1)) * 0.5
 mid_y = (max(x2)+min(x2)) * 0.5
 mid_z = (max(x3)+min(x3)) * 0.5
 ax.set_xlim(-1,1)
 ax.set_ylim(-1,1)
 ax.set_zlim(-1,1)
 ax.set_title(title,fontweight='bold')
 
 plt.savefig(figname)
# plt.show()
# plt.close(fig)

def plot_fig(x, y, markevery, title, xlabel, ylabel, figname):
	font = {'weight' : 'bold',
				'size'   : 18}
	fig = plt.figure(figsize=(12,12))
	plt.autoscale(enable=True, axis='x', tight=True)
#	miny = min(np.ndarray.flatten(np.asarray(y)))
#	maxy = max(np.ndarray.flatten(np.asarray(y)))
#	plt.ylim(1.3*miny,1.3*maxy)
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

def plot_norm(x,y,markevery,title,xlabel,ylabel,figname):
	font = {'weight' : 'bold',
				'size'   : 18}
	fig = plt.figure(figsize=(12,12))
	plt.autoscale(enable=True, axis='x', tight=True)
	miny = min(np.ndarray.flatten(np.asarray(y)))
	maxy = max(np.ndarray.flatten(np.asarray(y)))
	plt.ylim(miny,maxy)
	mpl.rc('font',**font)
	plt.plot(x,y,label=r'$|'+ylabel+'|$',lw=3.0,markevery=markevery)
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

def plot_N(x,y,mag,title,xlabel,ylabel,figname):
	font = {'weight' : 'bold',
				'size'   : 18}
	fig = plt.figure(figsize=(12,12))
	plt.autoscale(enable=True, axis='x', tight=True)
	plt.ylim(-1.3*mag,1.3*mag)
	mpl.rc('font',**font)
	y1,y2,y3 = zip(*y)
	plt.axhline(mag,lw=1.0,ls='--')
	plt.axhline(-mag,lw=1.0,ls='--')
	plt.plot(x,y1,label=r'$'+ylabel+'_{1}$',lw=3.0)
	plt.plot(x,y2,label=r'$'+ylabel+'_{2}$',lw=3.0)
	plt.plot(x,y3,label=r'$'+ylabel+'_{3}$',lw=3.0)
	plt.legend()
	plt.xlabel(r'$'+xlabel+'$')
	plt.ylabel(r'$'+ylabel+'$')
	plt.title(title,fontweight='bold')
	plt.savefig(figname)
	plt.close(fig)

if __name__ == "__main__":
 #This block will be evaluated if this script is called as the main routine
 #and will be ignored if this file is imported from another script.

 pth = 'figs'
 make_sure_path_exists(pth)

 log = open('log','r')
 with cd(pth):
  # Set up data lists

  w = []
  wb = []
  J = []
  Jb = []
  N = []
  F = []
  t = []
  R = []
  Qt = []

  x = []
  v = []
  
  magtol = 1e-3
  
  string = log.readlines()[3]
  Nmax = [int(s) for s in string.split() if s.isdigit()]
  Nmax = Nmax[0]
  log.seek(0)
  string = log.readlines()[1]
  k = [float(s) for s in string.split()[2:5]]
  log.seek(0)
  Ip = np.diag(map(float, log.readlines()[15].split()))
  log.seek(0)
  Q = np.genfromtxt(islice(log,17,20))
  log.seek(0)
  I = np.dot(Q,Ip)
  nmag = Round_To_n(np.linalg.norm(np.dot(I,np.asarray([0,0,magtol])))/3,0)
  
  mass = float(log.readlines()[11].split()[2])
  log.seek(0)
  Fmag = mass*magtol
  
  lines = log.readlines()[-10000:]
#  lines = log.readlines()[22:100000]
#  lines = log.readlines()[39000:40000]
    
  for line in lines:
   data = line.split("|")
   x.append(map(float,data[1].split()))
   v.append(map(float,data[2].split()))
   w.append(map(float,data[3].split()))
   JJ = map(float,data[4].split())
   J.append(JJ)
   N.append(map(float,data[5].split()))
   F.append(map(float,data[6].split()))
   t.append(map(float,data[7].split()))
   R.append(map(float,data[8].split()))
   RR = np.reshape(map(float,data[8].split()),(3,3))
   JJ = np.dot(np.transpose(np.dot(RR, Q)),JJ)
   ww = np.dot(np.transpose(np.dot(RR,Q)), map(float,data[3].split()))
   wlen = math.sqrt(sum([i*i for i in ww]))
   if wlen != 0:
    ww = ww/wlen
    JJ = JJ/np.linalg.norm(JJ)
   else:
    ww = [0,0,0]
    JJ = [0,0,0]
   wb.append(ww)
   Jb.append(JJ)
   Qt.append(map(float,data[4].split())-np.dot(RR,Q)[:,2])
	
  log.close()
  markevery = round(Nmax/1000)
  mav = int(round(Nmax/100))
  if mav < 2: mav = 2
  if markevery < 1: markevery = 1
  wnorm = [math.sqrt(sum([i*i for i in vec])) for vec in w]
  
#  plot_fig(t,x,markevery,'Position vs. Time','t (s)','x','x.png')
#  plot_fig(t,v,markevery,'Velocity vs. Time','t','v','v.png')
  plot_fig(t,w,markevery,'Angular velocity vs. Time','t (s)','\omega','w.png')
#  plot_fig(t,wb,markevery,'Angular velocity in body frame vs. Time','t','\omega','wb.png')
#  plot_norm(t,wnorm,markevery,r'$|\omega|$ vs. Time','t','\omega','wnorm.png')
  plot_fig(t,J,markevery,'Angular momentum vs. Time','t (s)','J (Nms)','J.png')
#  plot_fig(t,Jb,markevery,'Angular momentum in body frame vs. Time','t','J','Jb.png')
#  plot_fig(t,Qt,markevery,'Evolution of internal alignment','t','(J-Q_{max})','Qt.png')
  plot_fig(t,F,markevery,'Force vs. Time','t (s)','F (N)','F.png')
  plot_fig(t,N,markevery,'Torque vs. Time','t (s)','N (Nm)','N.png')
#  plot_N(t[1:],N[1:],nmag,'Torque vs. Time','t','N','N2.png')
#  plot_N(t[1:],F[1:],Fmag,'Force vs. Time','t','F','F2.png')
#  plot_orbit(x,'Translation','xorbit.png')
  plot_orbit(w,'Spin of angular velocity','worbit.png')
  plot_orbit(wb,'Spin of angular velocity in body frame','wborbit.png')
  plot_orbit(J,'Time evolution of angular momentum','Jorbit.png')
  plot_orbit(Jb,'Angular momentum direction in body frame','Jborbit.png')
  plot_R(R,Q,k,'Time evolution of','.png')
  plot_mav(t,F,mav,'Running average of Force vs. Time','t (s)','F (Nm)','Fav.png')
  plot_mav(t,w,mav,'Running average of angular velocity vs. Time','t (s)','\omega (rad/s)','wav.png')
  plot_mav(t,N,mav,'Running average of torque vs. Time','t (s)','N (Nm)','Nav.png')

