import numpy as np, sys, pymesh, h5py, matplotlib as mpl, getopt
from scipy.special import spherical_in, factorial, lpmn
from numpy.linalg import norm as vlen
from numpy.random import seed, normal
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d

# Parameter declarations
cflg = 1        # Correlation function C_1 = power law, C_2 = Gauss
sigma = 0.285   # Relative std of radial distance
nu = 4.5        # Power law index for C_1 correlation
gamma = 35.0    # Input angle for C_2 correlation
lmin = 2        # Minimum degree in C_1, C_2
lmax = 10

refr_ind = 1.4588 + 0.001j

refr_r = np.real(refr_ind)
refr_i = np.imag(refr_ind)

# For vectorizing we want to, well, use vectors. Here we handle indexing between
# arrays and vectors (both ways).
def c_ind(n, m=None):
   if(m==None):
      nn = np.floor((np.sqrt(np.asarray(n)*8+1)-1)/2)
      nn = nn.astype(int)
      return nn, n-(nn**2+nn)/2
   nn = np.asarray(n)*(np.asarray(n)+1)/2+np.asarray(m)
   return nn
   
cind = np.arange(lmax*(lmax+1)/2, dtype=int)
l, m = c_ind(cind)

# Define power law and gaussian autocorrelation functions
def power_lcoef(l, nu):
   a_l = np.zeros(l.size)
   for i,value in enumerate(l):
      if(value!=0):
         a_l[i] = 1.0/value**nu
   a_l[0:lmin+1]=0
   norm = np.sum(1.0/np.arange(lmin,lmax+1)**nu)
   return a_l/norm

def mgaussian_lcoef(l, ell):
   a_l = 2.0*(l+1)*spherical_in(l,1.0/ell**2)*np.exp(-1.0/ell**2)
   a_l[0:lmin+1]=0.
   norm = np.sum(a_l)
   return a_l/norm
   
# Determine the standard deviation from the autocorrelation
def coef_std(a_l,beta,l):
   std = beta*np.sqrt(a_l)
   for i in cind:
      j = l[i]*(l[i]+1)/2
      std[i] = std[j]*np.sqrt(2.*factorial(l[i]-m[i])/factorial(l[i]+m[i]))
   return std

# The spherical harmonics representation of the G-sample-sphere from the 
# statistics
def sample_gaussian_sphere_coef(std, l):
   a_lm = normal()*std
   b_lm = np.zeros(std.size)
   for mm in range(1,lmax):
      for ll in range(max(mm,lmin),lmax):
         i = c_ind(ll,mm)
         a_lm[i] = normal()*std[i]
         b_lm[i] = normal()*std[i]         
   return a_lm, b_lm
   
# Generate the Gaussian deformation using the spherical harmonics expansion
def deform(sphere, a_lm, b_lm, beta):
   vertices = 0.0*sphere.nodes
   # Node coordinates:
   for i in range(sphere.num_nodes):
      z = sphere.nodes[i,2]
      phi = np.arctan2(sphere.nodes[i,1],sphere.nodes[i,0])
      r = r_gsphere(a_lm, b_lm, z, phi, beta)
      nu = np.sqrt(1.0-z**2)
      vertices[i,:] = [r*nu*np.cos(phi),r*nu*np.sin(phi),r*z]
   
   return vertices

# Find the radius functions of the g-sphere
def r_gsphere(a_lm, b_lm, z, phi, beta):
   if(lmax==0):
      return a_lm[0,0]
   legp, _ = lpmn(lmax, lmax, z)
   CPHI = np.cos(np.arange(0,lmax+1)*phi)
   SPHI = np.sin(np.arange(0,lmax+1)*phi)
   s = 0.
   for ll in range(lmin,lmax):
      ii = c_ind(ll,0)
      s = s + legp[ll,0]*a_lm[ii]
   for mm in range(1,lmax):
      for ll in range(max(mm,lmin),lmax):    
         ii = c_ind(ll,mm)
         s = s + legp[mm,ll]*(a_lm[ii]*CPHI[mm]+b_lm[ii]*SPHI[mm])
   return np.exp(s-0.5*beta**2)

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
   
   param_r = refr_r*np.ones(tetramesh.voxels.shape[0])
   param_i = refr_i*np.ones(tetramesh.voxels.shape[0])
   
   with h5py.File(meshname+".h5","w") as f:
      dset1 = f.create_dataset("node", tetramesh.vertices.shape, dtype='double' )
      dset1[...] = tetramesh.vertices
      dset2 = f.create_dataset("elem", tetramesh.voxels.shape, dtype='int32')
      dset2[...] = tetramesh.voxels+1
      dset3 = f.create_dataset("refr_r", param_r.shape, dtype='double')
      dset3[...] = param_r
      dset4 = f.create_dataset("refr_i", param_i.shape, dtype='double')
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
   gid   = 1      
   gotGid = False
   refinement = 0.15
   try:
      opts, args = getopt.getopt(argv,"i:o:f:r:")
   except getopt.GetoptError:
      print('generate_gsphere.py -i <G-ID> -o <meshname> -r <refinement_level>')
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
   meshname, gid, refinement  = args(sys.argv[1:])

   gamma = gamma*np.pi/180
   ell = 2.0*np.sin(0.5*gamma)
   beta = np.sqrt(np.log(sigma**2+1.0))

   # Autocorrelation
   if(cflg == 1):
      a_l = power_lcoef(l, nu)
   if(cflg == 2):
      a_l = mgaussian_lcoef(l, ell)
   
   seed(gid)
   std = coef_std(a_l, beta, l)

   # Generate a sample Gaussian sphere with id gid, then move to discretization and output
   sphere = pymesh.generate_icosphere(1.0,[0.,0.,0.],3)
   a_lm, b_lm = sample_gaussian_sphere_coef(std, l)
   new_vertices = deform(sphere, a_lm, b_lm, beta)
   gsphere = pymesh.form_mesh(new_vertices,sphere.elements)
   
   pymesh.save_mesh(meshname+".ply", gsphere)
   tetramesh = draw_mesh(meshname, gsphere, refinement)
   pymesh.save_mesh(meshname+".mesh",tetramesh)
