import numpy as np, matplotlib as mpl, h5py, sys, pymesh, getopt
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d
from scipy.special import spherical_in, factorial, lpmn
from numpy.random import normal, seed
from numpy import linalg as la

sigma = 0.125  # Standard deviation
ell   = 0.35   # Correlation length
a     = 1.00   # Ellipsoid semiaxis a
b     = 0.80   # Ellipsoid semiaxis b
c     = 0.50   # Ellipsoid semiaxis c
ref   = 3      # Mesh refinement level
gid   = 1      # G-ellipsoid id


ref_ind = 1.686 + 0.0312j

eps_r = np.real(ref_ind**2)
eps_i = np.imag(ref_ind**2)

h     = a*(c/a)**2
sigma = sigma/h
beta  = np.sqrt(np.log(sigma**2+1.0))

# Gaussian correlation
def corr(x,ell):
   if(x/ell<10.):
      return np.exp(-0.5*x**2/ell**2)
   return 0.0

# Generate Gaussian random deviates.
def rand_gauss(cv,n):
   D = np.zeros((n,n))
   
   # Note: cv is assumed to be positive-definite
   eigval, eigvec = la.eigh(cv) 
   
   # In ascending order and drop everything ill-conditioned
   e = eigval[::-1]
   v = np.fliplr(eigvec)
   v[:,e<0.] = 0.0
   e[e<0.] = 0.0
   
   for j in range(n):
      for i in range(n):
         D[i,j] = np.sqrt(cv[i,i]*e[j])*v[i,j]
   
   rn = normal(size=n)
   return np.dot(D,rn)

def deform_surf(mesh):
   n = mesh.num_vertices
   cv = np.diag(np.ones(mesh.num_vertices))
   
   for i in range(n):
      for j in range(n):
         d = la.norm(mesh.vertices[i,:]-mesh.vertices[j,:])
         cv[i,j] = corr(d, ell)
         cv[j,i] = cv[i,j]
   
   h1 = rand_gauss(cv,n)  
   hn = h*np.exp(beta*h1-0.5*beta**2)
        
   return hn

def genellip():
   sphere = pymesh.generate_icosphere(1.0,[0.,0.,0.],ref)
   nodes = sphere.vertices
   M = np.diag([a,b,c])
   nodes = np.dot(M,nodes.T).T
   ellipsoid = pymesh.form_mesh(nodes,sphere.elements)
   return ellipsoid

def deform_mesh(mesh):
   # Compute the deformation of surface normals
   hn = deform_surf(mesh)
   
   mesh.add_attribute("vertex_normal")
   nn = mesh.get_vertex_attribute("vertex_normal")
   X = np.zeros((mesh.num_vertices,3))
   # Node coordinates:
   for i in range(mesh.num_vertices):
      X[i,:] = mesh.vertices[i,:] + (hn[i]-h)*nn[i,:]
         
   return X

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

def draw_mesh(meshname, mesh, refinement):
   # Generate the tetrahedral 3d mesh for the particle. Note that optimal 
   # refinement can be hard to automatize with tetgen, so quartet is used!
   tetramesh = pymesh.tetrahedralize(mesh,refinement,engine='quartet')     
   print('Number of tetras: ' + str(tetramesh.num_elements))
   
   param_r = eps_r*np.ones(tetramesh.voxels.shape[0])
   param_i = eps_i*np.ones(tetramesh.voxels.shape[0])
   
   with h5py.File(meshname+".h5","w") as f:
      dset1 = f.create_dataset("coord", tetramesh.vertices.shape, dtype='double' )
      dset1[...] = tetramesh.vertices
      dset2 = f.create_dataset("etopol", tetramesh.voxels.shape, dtype='int32')
      dset2[...] = tetramesh.voxels+1
      dset3 = f.create_dataset("param_r", param_r.shape, dtype='double')
      dset3[...] = param_r
      dset4 = f.create_dataset("param_i", param_i.shape, dtype='double')
      dset4[...] = param_i
   fig = plt.figure(figsize=(6, 6),frameon=False)
   ax = mplot3d.Axes3D(fig)
   
   # Collect face data as vectors for plotting
   F = boundary_faces(tetramesh.elements)
   
   facevectors = np.zeros((F.shape[0],3,3))
   for i, face in enumerate(F):
      for j in range(3):
         facevectors[i][j] = tetramesh.vertices[face[j],:]
   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(facevectors, facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.66))
   
   scale = tetramesh.vertices.flatten(-1)
   ax.auto_scale_xyz(scale, scale, scale)
   
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)

   plt.savefig(meshname)
   return tetramesh

def args(argv):
   meshname = "mesh"
   gid   = 1      # G-ellipsoid id
   gotGid = False
   refinement = 0.2
   try:
      opts, args = getopt.getopt(argv,"i:o:f:r:")
   except getopt.GetoptError:
      print('generate_gellip.py -i <G-ID> -o <meshname> -r <refinement_level>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-i'):
         try:
            gid = int(arg)
            gotGid = True
         except ValueError:
            print('Argument of -i is an integer, so try again!')
            sys.exit(2)
      if opt in ('-r'):
         try:
            refinement = float(arg)
         except ValueError:
            print('Argument of -r is a float, so try again!')
            sys.exit(2)
      if opt in ('-o'):
         meshname = arg
   if gotGid:
      meshname = meshname + str(gid)
   return meshname, gid, refinement

if __name__ == "__main__":
   meshname, gid, refinement = args(sys.argv[1:])
   seed(gid)
   ellipsoid = genellip()  

   # Discrete triangle representation for a sample G-sphere.
   node = deform_mesh(ellipsoid)
   gellip = pymesh.form_mesh(node,ellipsoid.elements)
   
   tetramesh = draw_mesh(meshname, gellip, refinement)
#   pymesh.save_mesh(meshname+".mesh",tetramesh)
      
