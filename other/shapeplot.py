import h5py, numpy as np, matplotlib as mpl, sys, getopt
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d
   
def read_mesh(meshname):
    meshfile = h5py.File(meshname,"r")
    V = np.asarray(meshfile['coord'][:])
    F = np.asarray(meshfile['etopol'][:]).astype(int)-1

    facevectors = np.zeros((F.shape[0],3,3))
    for i, face in enumerate(F):
        for j in range(3):
            facevectors[i][j] = V[face[j],:]

    return V, F, facevectors

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
   V, F, facevectors = read_mesh(meshname+".h5")

   fig = plt.figure(figsize=(6, 6),frameon=False)
   ax = mplot3d.Axes3D(fig)
   
   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(facevectors, facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.8))
   
   scale = V.flatten('F')
   ax.auto_scale_xyz(scale, scale, scale)
   
#   plt.axis('off')   
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)
   plt.show()
#   mymesh.save(meshname+'.stl')
#   plt.savefig(meshname)
