import numpy as np, h5py, sys, getopt

# Change tetrahedral mesh elements layer by layer
def layers(elem, depths):
#  input:
#      elem: Tetrahedral element list
#      depths: Layer thicknesses from the surface inwards
#  example:

#      elem2 = layers(elem,[3,2]);
#          Creates new element list elem2, where the 5th columns
#          have been changed so that outermost layer consists of
#          the surface tetrahedra and the next 3-1 neighbors of
#          those tetrahedra. Next layer is the surface of the core
#          and their neighbors.
   L = len(depths)
   ine = elem
   newelem = np.array([[0,0,0,0,0]])
   for i in range(L):
      depth = depths[i]
      for j in range(depth):
         oute,outi = outlayer(ine)
         ine = np.delete(ine,outi,axis=0)
         newelem = np.append(newelem, oute, axis=0)
      newelem[:,4] = newelem[:,4]+1
   newelem = np.append(newelem,ine,axis=0)
   newelem = np.delete(newelem,0,axis=0)
   return newelem

def outlayer(elem):
   f = elem[:,range(4)]
   edges = np.concatenate((f[:,[0,1,2]],
            f[:,[1,0,3]],
            f[:,[0,2,3]],
            f[:,[1,3,2]]),axis=0)
   edgesort = np.sort(edges,1)
   _, i = np.unique(edgesort,axis=0,return_index=True)
   _, j = np.unique(edgesort,axis=0,return_inverse=True)
   vec = np.bincount(j)

   qx = np.argwhere(vec==1)
   
   ii = np.mod(i,elem.shape[0])
   
   outind = np.unique(ii[qx])
   oute = elem[np.unique(ii[qx]),:]
   
   return oute, outind

def read_mesh(meshname):
   meshfile = h5py.File(meshname,"r")
   V = np.asarray(meshfile['coord'][:])
   T = np.asarray(meshfile['etopol'][:],dtype=int)

   return V, T

def args(argv):
   meshname = "mesh"
   outmesh = "layermesh"
   try:
      opts, args = getopt.getopt(argv,"i:o:")
   except getopt.GetoptError:
      print('layered.py -i <meshname> -o <layername>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-i'):
         meshname = arg
      if opt in ('-o'):
         outmesh = arg
   return meshname, outmesh

if __name__ == "__main__":
   meshname, outmesh = args(sys.argv[1:])
   V,T = read_mesh(meshname+".h5")
   if T.shape[1]<5:
      TT = np.ones([T.shape[0],T.shape[1]+1],dtype=int)
      TT[:,:-1] = T
      T = TT
      
   elem = layers(T,[3,8])
   refr_r = [1.686,1.6861,2.]
   refr_i = [0.0312,0.0312,0.003]
   param_r = np.zeros(elem.shape[0])
   param_i = np.zeros(elem.shape[0])
   mx = np.max(elem[:,4])
   
   for i in range(mx):
      ind = elem[:,4]==i+1
      param_i[ind] = refr_i[i]
      param_r[ind] = refr_r[i]
   
   with h5py.File(outmesh+".h5","w") as f:
      dset1 = f.create_dataset("node", V.shape, dtype='double' )
      dset1[...] = V
      dset2 = f.create_dataset("elem", elem.shape, dtype='int32')
      dset2[...] = elem
      dset3 = f.create_dataset("refr_r", param_r.shape, dtype='double')
      dset3[...] = param_r
      dset4 = f.create_dataset("refr_i", param_i.shape, dtype='double')
      dset4[...] = param_i
      
