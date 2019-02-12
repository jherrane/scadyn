import h5py, numpy as np
import matplotlib as mpl
from matplotlib import colors as col, pyplot as plt, rc, rcParams
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy.linalg import norm
from scipy.special import eval_genlaguerre as Laguerre
from scipy.misc import factorial
from numpy import pi, cos, sin, arctan2 as atan2, sqrt, abs, exp
rcParams.update({'figure.autolayout': True})

# Parameter declarations
readfile = False
# If read from file
gridname = "grid_xyz.h5"
Ename = "E_field.h5"
labels=["x","y","z"]
# If not
p = 0
l = 3
lmbda = 1.064
w0 = 0.5
lim = 3
zeroXYZ = 2 # Choose the zero plane in python indices, 2: Z = 0

def LG_field(p,l,lmbda,w0,lim):
   rnge = lim*np.linspace(-w0,w0,101)
   X,Y,Z = np.meshgrid(rnge,rnge,rnge)
   
   zR = pi*w0**2/lmbda
   k = 2*pi/lmbda
   ks = sqrt(k/zR)
   wz = w0*sqrt(1+(Z/zR)**2)
   gouy = atan2(Z,zR)
   upl = (-1.)**p*ks*sqrt(factorial(p)/(pi*factorial(p+abs(l))))
   
   Rz = Z**2+zR**2
   
   r = sqrt(X**2+Y**2)
   R = sqrt(2)*sqrt(X**2+Y**2)/wz
   phi = atan2(Y,X)
   
   Lpl = Laguerre(p,abs(l),R**2)
   E = upl*sqrt(2)/(ks*wz)*exp(-0.5*R**2)*R**(abs(l))*Lpl*exp(-1j*l*phi)*exp(1j*(2*pi+abs(l)+1)*gouy-1j*(k*r**2/(2*Rz))+1j*k*Z)
   return X,Y,Z,E
   
def make_segmented_cmap(): 
   white = '#ffffff'
   black = '#000000'
   red = '#ff0000'
   blue = '#0000ff'
   anglemap = col.LinearSegmentedColormap.from_list(
     'anglemap', [black, red, white, blue, black], N=256, gamma=1)
   return anglemap
            
def read_h5(gridname, Ename):
   grid = h5py.File(gridname,"r")
   E = h5py.File(Ename,"r")
   grid = np.asarray(grid['A_r'][:])
   E_r = np.asarray(E['A_r'][:])
   E_i = np.asarray(E['A_i'][:])
   E = E_r +1j*E_i
   E2 = np.power(np.abs(E),2)
   
   # Determine which projection is used and set as x_i, find final grid
   x_i = np.nonzero(grid[1,:])[0]
   x = grid[:,x_i[0]]
   y = grid[:,x_i[1]]
   n = int(np.sqrt(len(x)))
   x = x.reshape(n,n)
   y = y.reshape(n,n)
   
   # Determine the E-component parallel to x,y
   Ezero = 2#np.where(grid[1,:] == 0.0)[0]
   Ez = np.sum(E2,1)
   norm = np.max(Ez)
   Ez = Ez/norm
   Ez = Ez.reshape(n,n)

   Ephase = np.angle(E[:,Ezero])
   Ephase = Ephase.reshape(n,n)
   
   i_label = [x_i[0], x_i[1], Ezero]

   return x, y, Ez, Ephase, i_label

def afield(X,Y,Z,E):
   if(zeroXYZ==0):
      x = Y[X==0]
      y = Z[X==0]
      E = E[X==0]
      i_label = [1, 2, 0]
   elif(zeroXYZ==1):
      x = X[Y==0]
      y = Z[Y==0]
      E = E[Y==0]
      i_label = [0, 2, 1]
   else:
      x = X[Z==0]
      y = Y[Z==0]
      E = E[Z==0]
      i_label = [0, 1, 2]
   n = int(np.sqrt(len(x)))
   x = x.reshape(n,n)
   y = y.reshape(n,n)
   # Determine the E-component parallel to x,y
   E2 = np.power(np.abs(E),2)
   norm = np.max(np.sum(E2))
   Ez = E2/norm
   Ez = Ez.reshape(n,n)

   Ephase = np.angle(E)
   Ephase = Ephase.reshape(n,n)

   return x, y, Ez, Ephase, i_label

def make_figure(x,y,Ez,Ephase,i_label):
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   font = {'weight' : 'bold', 'size' : 48}
   
   mpl.rc('font',**font)
   
   cmap = plt.get_cmap('inferno')

   # Create figure for incident field
   fig = plt.figure(figsize=(12,10),frameon=False)
   ax = fig.gca()
   im = ax.pcolormesh(x, y, Ez, cmap=cmap, shading='gouraud')
   fig.colorbar(im, ax=ax)
   
   ax.set_title(r'Intensity')
   ax.set_xlabel(r'$'+labels[i_label[0]]+' (\lambda)$',labelpad=5) 
   ax.set_ylabel(r'$'+labels[i_label[1]]+' (\lambda)$',labelpad=5) 
   
   plt.savefig('E_field_'+labels[i_label[0]]+'_'+labels[i_label[1]], bbox_inches='tight', pad_inches = 0.05)
   
   cmap = make_segmented_cmap()
   
   fig = plt.figure(figsize=(12,10),frameon=False)
   ax = fig.gca()
   
   if(np.amin(Ephase)>-np.pi+0.01 and np.amax(Ephase)<np.pi-0.01):
      cmap = plt.get_cmap('binary')
      im = ax.pcolormesh(x, y, Ephase, cmap=cmap, shading='gouraud')
      cbar = fig.colorbar(im, ax=ax,ticks=[-np.pi/2, 0, np.pi/2])
      cbar.ax.set_yticklabels([r'$-\pi/2$', r'$0$', r'$\pi/2$'])
   else:
      im = ax.pcolormesh(x, y, Ephase, cmap=cmap, shading='gouraud')
      cbar = fig.colorbar(im, ax=ax,ticks=[-np.pi+0.01, -np.pi/2, 0, np.pi/2, np.pi-0.01])
      cbar.ax.set_yticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
   
   ax.set_title(r'Phase')
   ax.set_xlabel(r'$'+labels[i_label[0]]+' (\lambda)$',labelpad=5) 
   ax.set_ylabel(r'$'+labels[i_label[1]]+' (\lambda)$',labelpad=5) 
   
   plt.savefig('E_phase_'+labels[i_label[0]]+'_'+labels[i_label[1]], bbox_inches='tight', pad_inches = 0.05)
   return
 
if __name__ == "__main__":
   if(readfile):
      x, y, Ez, Ephase, i_label = read_h5(gridname, Ename)
   else:
      X,Y,Z,E = LG_field(p,l,lmbda,w0,lim)
      x, y, Ez, Ephase, i_label = afield(X,Y,Z,E)
   
   make_figure(x, y, Ez, Ephase, i_label)

