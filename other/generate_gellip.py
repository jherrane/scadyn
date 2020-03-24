# Copyright 2019 Joonas Herranen (joonas.herranen@iki.fi)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Based on description of Gaussian random ellipsoids by Muinonen, K. & Pieniluoma, T. Light scattering by Gaussian random ellipsoid particles: First results with discrete-dipole approximation Journal of Quantitative Spectroscopy & Radiative Transfer, 2011, 112, 1747-1752 (https://doi.org/10.1016/j.jqsrt.2011.02.013).

# The code generates both a surface triangle mesh by deforming an ellipsoidal mesh according to the Gaussian random ellipsoid rules and from that generates and outputs a tetrahedral volume mesh. The main use of this code is to generate inputs for a volume integral equation solution of electromagnetic scattering, and for that reason the .h5 output includes the refractive indices for each tetrahedra (this code generates homogeneous shapes).

# The most important dependency is PyMesh (https://pymesh.readthedocs.io/en/latest/), which can be used to manipulate and save both triangular and tetrahedral meshes. For extreme shapes, the default choices for tetrahedralization and refinement levels should probably be changed.

import numpy as np, matplotlib as mpl, h5py, sys, pymesh, getopt
from matplotlib import rc, rcParams, pyplot as plt
from mpl_toolkits import mplot3d
from scipy.special import spherical_in, factorial, lpmn
from numpy.random import normal, seed
import meshio
from numpy import linalg as la
from numpy.linalg import norm
from scipy.spatial import ConvexHull as ch

sigma = 0.125  # Standard deviation
ell   = 0.35   # Correlation length
a     = 1.00   # Ellipsoid semiaxis a
b     = 0.80   # Ellipsoid semiaxis b
c     = 0.50   # Ellipsoid semiaxis c
n     = 1000    # Surface mesh vertex count
gid   = 1      # G-ellipsoid id
aeff = 5.2
lmbda = 0.500

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

def sphere_points(n=1000):
   x = np.zeros([3,n])
   offset = 2./n
   increment = np.pi * (3. - np.sqrt(5.));

   for i in range(n):
      yi = ((i * offset) - 1) + (offset / 2);
      r = np.sqrt(1 - yi**2)
      phi = ((i + 1) % n) * increment
      
      x[:,i] = [np.cos(phi) * r, yi, np.sin(phi) * r]
   return x.T
   
def surface_fix(V,F):
   n = np.zeros([3,np.shape(F)[0]])
   centroids = np.zeros([3,np.shape(F)[0]])
   for i, face in enumerate(F):
      n[:,i] = np.cross(V[face[1],:]-V[face[0],:],V[face[2],:]-V[face[0],:])/2.
      centroids[:,i] = (V[face[0],:] + V[face[1],:] + V[face[2],:])/3
      if np.dot(n[:,i], centroids[:,i]) < 0:
         f0 = face[0]; f1 = face[1]; f2 = face[2]
         F[i,0] = f2
         F[i,1] = f1
         F[i,2] = f0
   return F

def genellip():
   nodes = sphere_points(n)

   elements = ch(nodes).simplices
   elements = surface_fix(nodes, elements)
   M = np.diag([a,b,c])
   nodes = np.dot(M,nodes.T).T
   ellipsoid = pymesh.form_mesh(nodes,elements)
#   ellipsoid, info = pymesh.split_long_edges(ellipsoid, 0.1)
   return ellipsoid

def deform_mesh(mesh):
   # Compute the deformation of surface normals
   hn = deform_surf(mesh)
   
   mesh.add_attribute("vertex_normal")
   nn = mesh.get_vertex_attribute("vertex_normal")
   print('normals calculated')
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
   return

def save_h5(meshname, mesh):
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
   return
   
def get_tetra_vol(mesh):
   F = mesh.elements
   V = 0.
   for i, face in enumerate(F):
      p0 = mesh.vertices[face[0],:]
      p1 = mesh.vertices[face[1],:]
      p2 = mesh.vertices[face[2],:]
      p3 = mesh.vertices[face[3],:]
      V += np.abs(np.dot(p0-p3,np.cross(p1-p3,p2-p3)))/6.
   return V
   
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

def write_msh(fname, mesh):
   f = open(fname,"w")
   f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
   f.close
   f = open(fname,'a')
   f.write(str(mesh.num_vertices)+'\n')
   for i, V in enumerate(mesh.vertices):
      f.write(str(i+1)+ ' '+ str(V[0])+ ' '+ str(V[1])+ ' '+ str(V[2])+'\n')
   f.write("$EndNodes\n$Elements\n")
   f.write(str(mesh.num_elements)+'\n')
   for i, V in enumerate(mesh.elements):
      f.write(str(i+1)+ ' 2 2 1 3 '+ str(V[0]+1)+ ' '+ str(V[1]+1)+ ' '+ str(V[2]+1)+'\n')
   f.write("$EndElements")
   f.close
   return
   
def truncation_order(lmbda,aeff):
   ka = 2.*np.pi/lmbda*aeff
   print()
   if ka<=1.:
      return 4.
   else:
      return np.floor(ka+4*ka**(1./3.))

def mesh_vectors(V,F):
   """
   Creates a vector set of the mesh data for nice plotting.
   """
   msh = np.zeros((np.shape(F)[0],3,3))
   for i, face in enumerate(F):
      for j in range(3):
         msh[i][j] = V[face[j],:]
   return msh

def plot_mesh(points, tris):
   fig = plt.figure(figsize=(8,8))
   ax = mplot3d.Axes3D(fig)
   
   meshvectors = mesh_vectors(points, tris)
   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(meshvectors, facecolor=[0.5,0.5,0.5], lw=0.5, edgecolor=[0,0,0], alpha=.8, antialiaseds=True))  
   scale = points.flatten('F')
   ax.auto_scale_xyz(scale, scale, scale) 

   plt.show()
   return
   
def fix_mesh(mesh, detail="normal"):
    bbox_min, bbox_max = mesh.bbox;
    diag_len = norm(bbox_max - bbox_min);
    if detail == "normal":
        target_len = diag_len * 2e-2;
    elif detail == "high":
        target_len = diag_len * 10e-3;
    elif detail == "low":
        target_len = diag_len * 3e-2;
    print("Target resolution: {} mm".format(target_len));

    count = 0;
    mesh, __ = pymesh.remove_degenerated_triangles(mesh, 100);
    mesh, __ = pymesh.split_long_edges(mesh, target_len);
    num_vertices = mesh.num_vertices;
    while True:
        mesh, __ = pymesh.collapse_short_edges(mesh, 1e-6);
        mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                preserve_feature=True);
        mesh, __ = pymesh.remove_obtuse_triangles(mesh, 150.0, 100);
        if mesh.num_vertices == num_vertices:
            break;

        num_vertices = mesh.num_vertices;
        print("#v: {}".format(num_vertices));
        count += 1;
        if count > 10: break;

    mesh = pymesh.resolve_self_intersection(mesh);
    mesh, __ = pymesh.remove_duplicated_faces(mesh);
    mesh = pymesh.compute_outer_hull(mesh);
    mesh, __ = pymesh.remove_duplicated_faces(mesh);
    mesh, __ = pymesh.remove_obtuse_triangles(mesh, 179.0, 5);
    mesh, __ = pymesh.remove_isolated_vertices(mesh);

    return mesh;
      
if __name__ == "__main__":
   meshname, gid, refinement = args(sys.argv[1:])
   seed(gid)
   ellipsoid = genellip()
   print('genellip() done')

   # Discrete triangle representation for a sample G-sphere.
   node = deform_mesh(ellipsoid)
   print('deform_mesh() done')
   gellip = pymesh.form_mesh(node,ellipsoid.elements)
   gellip = fix_mesh(gellip,detail="high")
   print('fix_mesh() done')
#   plot_mesh(gellip.vertices, gellip.elements)
   
   # Generate the tetrahedral 3d mesh for the particle. Note that optimal 
   # refinement can be hard to automatize with tetgen, so quartet is used!
   tetramesh = pymesh.tetrahedralize(gellip,refinement,engine='quartet')   
   F2 = boundary_faces(tetramesh.elements)
#   print('Number of tetras: ' + str(tetramesh.num_elements))
   V = get_tetra_vol(tetramesh)
   gellip_s = pymesh.form_mesh(gellip.vertices*aeff/(3.*V/4./np.pi)**(1./3.),gellip.elements)

   write_msh(meshname+"_s.msh", gellip_s)   
   print('Use l_max = ', truncation_order(lmbda,aeff))
   draw_mesh(meshname, gellip, refinement)
   save_h5(meshname, tetramesh)
   pymesh.save_mesh(meshname+".mesh",tetramesh)
      
