import h5py
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
from matplotlib import rc, rcParams
rcParams.update({'figure.autolayout': True})
            
def read_mesh(meshname):
   mesh = h5py.File(meshname,"r")
   coord = np.asarray(mesh['coord'][:])
   etopol = np.asarray(mesh['etopol'][:])
   faces = np.zeros([len(etopol),3],dtype=int)
   etopol = etopol - 1
   for i in range(0,len(etopol)):
      [v0,v1,v2,v3] = coord[etopol[i,:]]
      verts = np.transpose([v0,v1,v2,v3])
      idx = np.argsort([norm(v0), norm(v1), norm(v2), norm(v3)])
      outerface = verts[:,idx[1:4]]
      faces[i,:] = etopol[i,idx[1:4]]
   return coord,etopol,faces
 
if __name__ == "__main__":
   meshname = "shape.h5"
   colors = ['r','k','b']

   # Read vertices (coord) and face data (etopol) from mesh hdf5 file
   coord, etopol, faces = read_mesh(meshname)

   # Create figure and 3D axes
   fig = plt.figure(figsize=(10, 10),frameon=False)
   ax = fig.add_subplot(111, projection='3d')
   
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   font = {'weight' : 'bold', 'size' : 48}
   mpl.rc('font',**font)
   
   # Set nice axis limits
   x,y,z = zip(*coord)
   mx = max(max(x),-min(x),max(y),-min(y), max(z),-min(z))
   ax.set_xlim(-mx, mx)
   ax.set_ylim(-mx, mx)
   ax.set_zlim(-mx, mx)

   points = coord
   
   mesh = Poly3DCollection(points[faces], facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.66)
   ax.add_collection3d(mesh)
   
   ax.set_xlabel('x',labelpad=10) 
   ax.set_ylabel('y',labelpad=10) 
   ax.set_zlabel('z',labelpad=10) 

   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)

   plt.savefig('shape')
