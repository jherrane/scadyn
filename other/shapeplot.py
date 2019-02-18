import h5py, numpy as np, matplotlib as mpl, sys, getopt
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d
from stl import mesh
            
def boundary_faces(T):
   T1 = np.array([T[:,0], T[:,1],T[:,2]]) 
   T2 = np.array([T[:,0], T[:,1],T[:,3]])
   T3 = np.array([T[:,0], T[:,2],T[:,3]]) 
   T4 = np.array([T[:,1], T[:,2],T[:,3]])

   T  = np.concatenate((T1,T2,T3,T4),axis=1)
   T = np.sort(T,axis=0)

   unique_cols, inverse = np.unique(T,axis=1, return_inverse=True)
   counts = np.bincount(inverse)==1
   F = unique_cols[:,counts] 
   
   return F.transpose()
   
def read_mesh(meshname):
   meshfile = h5py.File(meshname,"r")
   V = np.asarray(meshfile['node'][:])
   T = np.asarray(meshfile['elem'][:])-1
   F = boundary_faces(T)

   msh =  mesh.Mesh(np.zeros(F.shape[0], dtype=mesh.Mesh.dtype))
   for i, face in enumerate(F):
      for j in range(3):
         msh.vectors[i][j] = V[face[j],:]

   return msh

def args(argv):
   meshname = "mesh"
   try:
      opts, args = getopt.getopt(argv,"i:")
   except getopt.GetoptError:
      print('shapeplot.py -i <meshname>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-i'):
         meshname = arg
   return meshname

if __name__ == "__main__":
   meshname = args(sys.argv[1:])
   mesh = read_mesh(meshname+".h5")

   fig = plt.figure(figsize=(6, 6),frameon=False)
   ax = mplot3d.Axes3D(fig)
   
   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.vectors, facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.66))
   
   scale = mesh.points.flatten(-1)
   ax.auto_scale_xyz(scale, scale, scale)
   
#   plt.axis('off')   
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)

#   mymesh.save(meshname+'.stl')
   plt.savefig(meshname)
