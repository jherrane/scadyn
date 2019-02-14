import numpy as np, sys, pymesh, h5py, matplotlib as mpl
from scipy.special import spherical_in, factorial, lpmn
from numpy.linalg import norm as vlen
from numpy.random import seed, normal
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d

# Parameter declarations
cflg = 1       # Correlation function C_1 = power law, C_2 = Gauss
sigma = 0.285    # Relative std of radial distance
nu = 4.5       # Power law index for C_1 correlation
gamma = 35.0    # Input angle for C_2 correlation
lmin = 2       # Minimum degree in C_1, C_2
lmax = 10
gid = 1        # G-sphere id

eps_r = 1.95   # Real and imaginary parts of permittivity
eps_i = 0.786  

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

if __name__ == "__main__":
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
   
   fig = plt.figure(figsize=(6, 6),frameon=False)
   ax = mplot3d.Axes3D(fig)
   
   # Collect face data as vectors for plotting
   facevectors = np.zeros((gsphere.faces.shape[0],3,3))
   for i, face in enumerate(gsphere.faces):
      for j in range(3):
         facevectors[i][j] = gsphere.vertices[face[j],:]
   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(facevectors, facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.66))
   
   scale = gsphere.vertices.flatten(-1)
   ax.auto_scale_xyz(scale, scale, scale)
   
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)

   plt.savefig("mesh")
   
   # Generate the tetrahedral 3d mesh for the particle. Note that optimal 
   # tet_volume and mesh refinement are highly situational, hard to automatize!
   tetgen = pymesh.tetgen()
   tetgen.points = gsphere.vertices
   tetgen.triangles = gsphere.faces
   tetgen.max_tet_volume = 0.03
   tetgen.verbosity = 1
   tetgen.run()
   tetramesh = tetgen.mesh
   
   param_r = eps_r*np.ones(tetramesh.voxels.shape[0])
   param_i = eps_i*np.ones(tetramesh.voxels.shape[0])
   
   with h5py.File("tetramesh.h5","w") as f:
      dset1 = f.create_dataset("coord", tetramesh.vertices.shape, dtype='double' )
      dset1[...] = tetramesh.vertices
      dset2 = f.create_dataset("etopol", tetramesh.voxels.shape, dtype='int32')
      dset2[...] = tetramesh.voxels+1
      dset3 = f.create_dataset("param_r", param_r.shape, dtype='double')
      dset3[...] = param_r
      dset4 = f.create_dataset("param_i", param_i.shape, dtype='double')
      dset4[...] = param_i
      
      
      
