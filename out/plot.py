import math, codecs, os, errno, sys, getopt, StringIO, numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc, rcParams, pyplot as plt    
from mpl_toolkits.mplot3d import Axes3D
from itertools import islice
from io import StringIO

rcParams.update({'figure.autolayout': True})
file_suffix = ""

class cd:
   """Context manager for changing the current working directory"""
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

def plot_fig(x, y, title, xlabel, ylabel, figname, *args, **kwargs):
   ylim = kwargs.get('ylim', None)
   
   fig = plt.figure(figsize=(10,8))

   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   font = {'weight' : 'bold', 'size' : 36}
   mpl.rc('font',**font)
   
   plt.autoscale(enable=True, axis='x', tight=True)
   plt.tight_layout()
   
   y1,y2,y3 = zip(*y)
   
   label=ylabel.split("(")[0]
   if figname is 'x':
      plt.plot(x,y1,label=r'${x}$',lw=3.0,linestyle='--')
      plt.plot(x,y2,label=r'${y}$',lw=3.0,linestyle=':')
      plt.plot(x,y3,label=r'${z}$',lw=3.0,linestyle='-')
   elif figname is 'v':
      plt.plot(x,y1,label=r'$v_{x}$',lw=3.0, linestyle='--')
      plt.plot(x,y2,label=r'$v_{y}$',lw=3.0, linestyle=':')
      plt.plot(x,y3,label=r'$v_{z}$',lw=3.0, linestyle='-')
   else:
      plt.plot(x,y1,label=r'$'+label+'_{1}$',lw=3.0, linestyle='--')
      plt.plot(x,y2,label=r'$'+label+'_{2}$',lw=3.0, linestyle=':')
      plt.plot(x,y3,label=r'$'+label+'_{3}$',lw=3.0, linestyle='-')
   ax = plt.gca()
   if ylim is not None:
      if ylim.size>1:
         ax.set_ylim([ylim[0],ylim[1]])
      else:
         ax.set_ylim([-ylim,ylim])
   plt.xlabel(r'$'+xlabel+'$',labelpad=10)
   plt.ylabel(r'$'+ylabel+'$',labelpad=10)
   plt.legend(fontsize=20)

   ax.tick_params(direction='in', length=16, width=2, which='major',pad=10)
   ax.tick_params(direction='in', length=10, width=1, which='minor',pad=10)
   ax.tick_params(top=True, right=True,which='both')
   
   plt.title(title,fontweight='bold')
   plt.savefig(figname+file_suffix+'.png')
   plt.close(fig)

if __name__ == "__main__":
   inputfile = 'log'
   skip = 22
   last = 0 # Draw only last n lines, 0 if all
   first = 0 # Draw only first n lines, 0 if all
   
   argv = sys.argv[1:]
   try:
      opts, args = getopt.getopt(argv,"hi:e:f:")
   except getopt.GetoptError:
      print 'python plot.py -h'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h' or opt == '--help':
         print "args are hi:e:f:"
         sys.exit()
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-e"):
         last = int(arg)
      elif opt in ("-f"):
         first = int(arg)
   
   splt = inputfile.split("log")
   file_suffix = file_suffix + splt[-1]
   log = open(inputfile,'r')
   inputf = codecs.open(inputfile, encoding='utf-8').read()
   inputf = inputf.replace('|','')
   lines = np.loadtxt(StringIO(inputf), skiprows=skip)
   if(last!=0):
      lines = lines[last:-1,:]
   if(first!=0):
      lines = lines[0:first,:]
   
   I = np.genfromtxt(islice(log,15,16))
   I = np.diag(I)
   log.close()
   
   t = lines[:,1]
   maxt = np.max(t)
   tstr = ''
   if maxt>1e-6 and maxt<1e-3: 
      tstr = '\mu'
      t = t*1e6
   if maxt>1e-3 and maxt<1.0: 
      tstr = 'm'
      t = t*1e3
   if maxt>1.0: 
      tstr = ''
      
   x = lines[:,2:5]*1e6
   w = lines[:,5:8]
   v = lines[:,8:11]*1e6
   J = np.zeros(w.shape)
   N = lines[:,11:14]
   F = lines[:,14:17]
   R = lines[:,17:27]
      
   for i in range(0,w.shape[0]):
      RR = np.reshape(R[i],(3,3),order='F')
      J[i,:] = np.matmul(I,w[i,:])
      w[i,:] = np.matmul(RR,w[i,:])
      J[i,:] = np.matmul(RR,J[i,:])
      N[i,:] = np.matmul(RR,N[i,:])

   with cd('figs'):
      plot_fig(t,x,r'Position of CM','t (\mathrm{'+tstr+' s})','x (\mathrm{\mu m})','x')
      plot_fig(t,w,'Angular velocity','t \mathrm{('+tstr+' s)}','\omega (\mathrm{rad/s})','w')
      plot_fig(t,J,'Angular momentum','t \mathrm{('+tstr+' s)}','J (\mathrm{Ns})','J')
      plot_fig(t,v,markevery,'Velocity','t \mathrm{('+tstr+' s)}','v (\mathrm{\mu m/s})','v')
      plot_fig(t,N,markevery,'Torque','t \mathrm{('+tstr+' s)}','N (\mathrm{Nm})','N')
      plot_fig(t,F,markevery,'Force','t \mathrm{('+tstr+' s)}','F (\mathrm{N})','F')

