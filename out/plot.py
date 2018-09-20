# Module to plot the results of T-VIEDYN integration (WIP)
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from cd import *
from arrow import *
import StringIO
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from itertools import islice
from io import StringIO
from matplotlib import rc
import codecs
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

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
   
   font = {'weight' : 'bold', 'size' : 40}
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   
   ax.set_aspect('equal')
   ax.plot(x1,x2,x3)
   ax.plot(x1,x2,(mid_z - max_range)*np.ones(len(x3)),alpha=0.1,lw=3,color='k')
   ax.plot(x1,np.ones(len(x2))*(mid_y + max_range),x3,alpha=0.1,lw=3,color='k')
   ax.plot((mid_x - max_range)*np.ones(len(x1)),x2,x3,alpha=0.1,lw=3,color='k')

   ax.set_xlim(mid_x - max_range, mid_x + max_range)
   ax.set_ylim(mid_y - max_range, mid_y + max_range)
   ax.set_zlim(mid_z - max_range, mid_z + max_range)

   ax.set_title(title,fontweight='bold')
   plt.locator_params(nbins=4)
   plt.savefig(figname)

def plot_fig(x, y, markevery, title, xlabel, ylabel, figname, *args, **kwargs):
   ylim = kwargs.get('ylim', None)
   font = {'weight' : 'bold',
			   'size'   : 18}
   fig = plt.figure(figsize=(14,14))
   plt.autoscale(enable=True, axis='x', tight=True)
   plt.tight_layout()
   
   font = {'weight' : 'bold', 'size' : 40}
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')

   mpl.rc('font',**font)
   y1,y2,y3 = zip(*y)
   plt.plot(x,y1,label=r'$'+ylabel+'_{1}$',lw=3.0,markevery=markevery)
   plt.plot(x,y2,label=r'$'+ylabel+'_{2}$',lw=3.0,markevery=markevery)
   plt.plot(x,y3,label=r'$'+ylabel+'_{3}$',lw=3.0,markevery=markevery)
   axes = plt.gca()
   if ylim is not None:
      axes.set_ylim([-ylim,ylim])
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
   plt.tight_layout()
   
   font = {'weight' : 'bold', 'size' : 40}
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   
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

if __name__ == "__main__":
   #This block will be evaluated if this script is called as the main routine
   #and will be ignored if this file is imported from another script.
   inputfile = 'log'
   pth = 'figs'
   skip = 22
   #plt.xkcd()
   
   log = open(inputfile,'r')
   inputf = codecs.open(inputfile, encoding='utf-8').read()
   inputf = inputf.replace('|','')
   lines = np.loadtxt(StringIO(inputf), skiprows=skip)
      
   magtol = 1e-3

   string = log.readlines()[3]
   Nmax = [int(s) for s in string.split() if s.isdigit()]
   Nmax = Nmax[0]
   log.seek(0)
   string = log.readlines()[1]
   k = [float(s) for s in string.split()[2:5]]
   log.seek(0)
   I = np.genfromtxt(islice(log,15,16))
   I = np.diag(I)
   log.seek(0)
   Q = np.genfromtxt(islice(log,17,20))
   log.seek(0)
   # Even though named a, this is the wavelength, sorry
   a = np.genfromtxt(islice(log,10,11))
   a = a[3]*1e-9
   log.seek(0)
   mass = np.genfromtxt(islice(log,11,12))
   mass = mass[2]
   log.close()
   
   markevery = round(Nmax/1000)
   mav = int(round(Nmax/100))
   if mav < 2: mav = 2
   if markevery < 1: markevery = 1
   
   t = lines[:,1]
   x = lines[:,2:5]/a
   w = lines[:,5:8]
   v = lines[:,8:11]/a
   J = 0*w
   N = lines[:,11:14]
   F = lines[:,14:17]/mass
   R = lines[:,17:27]
      
   for i in range(0,w.shape[0]):
      J[i,:] = np.matmul(I,w[i,:])
      w[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),w[i,:]  ))
      N[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),N[i,:]  ))
      J[i,:] = np.matmul(Q,np.matmul(np.reshape(R[i,:],(3,3),order='F'),J[i,:]  ))
#      w[i,:] = 1.1*w[i,:]/np.sqrt(np.sum(np.power(w[i,:],2)))
#      J[i,:] = 1.1*J[i,:]/np.sqrt(np.sum(np.power(J[i,:],2)))

   Jb = np.zeros(J.shape)
   wb = np.zeros(w.shape)

   with cd(pth):
      plot_fig(t,x,markevery,'Position of CM vs. Time','t (s)','x (\lambda)','x.png',ylim=0.5)
      plot_fig(t,w,markevery,'Angular velocity vs. Time','t (s)','\omega','w.png')
      plot_fig(t,v,markevery,'Velocity vs. Time','t (s)','v (\lambda/s)','v.png')
#      plot_fig(t,J,markevery,'Angular momentum vs. Time','t (s)','J (Nms)','J.png')
      plot_fig(t,N,markevery,'Torque vs. Time','t (s)','N (Nm)','N.png')
      plot_fig(t,F,markevery,'Force vs. Time','t (s)','F (N)','F.png')
#      plot_orbit(w,'Spin of angular velocity','worbit.png')
#      plot_orbit(J,'Time evolution of angular momentum','Jorbit.png')
#      plot_mav(t,w,mav,'Running average of angular velocity vs. Time','t (s)','\omega (rad/s)','wav.png')
#      plot_mav(t,N,mav,'Running average of torque vs. Time','t (s)','N (Nm)','Nav.png')

