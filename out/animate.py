import h5py, numpy as np, matplotlib as mpl, moviepy.editor as mpy, codecs, pymesh
from matplotlib import rc, pyplot as plt
from mpl_toolkits import mplot3d
from stl import mesh
from moviepy.video.io.bindings import mplfig_to_npimage
from itertools import islice
from io import StringIO
import time

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def boundary_faces(V, T):
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

def decimate_mesh(tol, V, F):
   mesh = pymesh.form_mesh(V,F)
   vertices, faces, info = pymesh.collapse_short_edges_raw(V, F, tol,preserve_feature=True)
   return vertices, faces

def read_mesh(meshname):
   meshfile = h5py.File(meshname,"r")
   V = np.asarray(meshfile['coord'][:])
   T = np.asarray(meshfile['etopol'][:])-1
   F = boundary_faces(V, T)
   V,F = decimate_mesh(0.1,V,F)

   msh =  mesh.Mesh(np.zeros(F.shape[0], dtype=mesh.Mesh.dtype))
   for i, face in enumerate(F):
      for j in range(3):
         msh.vectors[i][j] = V[face[j],:]

   return msh.vectors

def read_log(lg, duration, fps):
   log = open(lg,'r')
   a = np.genfromtxt(islice(log,8,9))
   a = a[2]
   log.seek(0)
   Q = np.genfromtxt(islice(log,17,20))
   log.close()
   log.close()

   inputf = codecs.open(lg, encoding='utf-8').read()
   inputf = inputf.replace('|','')
   lines = np.loadtxt(StringIO(inputf), skiprows=22)
#   lines = lines[:10000,:]
   # First solution gives smooth motion with variable time, second constant clock
#   lines = lines[0::len(lines)/(duration*fps),:]
   t = lines[:,1]
   t0 = t[0]
   tf = t[-1]
   times = np.linspace(t0,tf,duration*fps)
   indices = np.zeros(np.size(times)).astype(int)
   for i, time in enumerate(times):
      idx = find_nearest(t,time)
      times[i] = t[idx]
      indices[i] = idx
      
   lines = lines[indices,:]
   
   t = lines[:,1]
   t0 = t[0]
   tf = t[-1]
   
   xx = lines[:,2:5]/a
   w = lines[:,5:8]
   R = lines[:,17:27]
   
   x,y,z = zip(*xx)
   mx =  np.max(np.abs(xx))
   
   return R, xx, t, mx, w, Q

def get_frame(i):
   i = int(i*fps)
   plt.clf()
   ax = fig.add_subplot(111, projection='3d')
   ax.set_xlim(-mx, mx)
   ax.set_ylim(-mx, mx)
   ax.set_zlim(-mx, mx)
   
   if(t[-1]<1.0):
      plt.title('t = '+str(t[i])+' s', fontsize=36)
   else:
      plt.title('t = '+time.strftime('%H:%M:%S', time.gmtime(t[i])), fontsize=36)

   RR = np.reshape(R[i],(3,3),order='F')
   points = P.dot(np.dot(Q,(RR.T))) + x[i]

   ax.add_collection3d(mplot3d.art3d.Poly3DCollection(points, facecolor=[0.5,0.5,0.5], lw=0.5,edgecolor=[0,0,0], alpha=0.66))
  
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)
   plt.setp( ax.get_zticklabels(), visible=False)

   return mplfig_to_npimage(fig)
 
if __name__ == "__main__":
   meshname = 'mesh'
   logname = 'log'
   
   duration = 30
   fps = 60
   
   P = read_mesh(meshname+'.h5')
   fig = plt.figure(figsize=(6,6),frameon=False)
   R, x, t, mx, w, Q = read_log(logname, duration, fps)

   animation =mpy.VideoClip(get_frame, duration=duration)
   animation.write_videofile(meshname+".mp4", fps=fps)
