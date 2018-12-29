import numpy as np
import sys, getopt
 
if __name__ == "__main__":   
   filename = 'trajectory'
   fin = 'path'
   
   argv = sys.argv[1:]
   try:
      opts, args = getopt.getopt(argv,"ho:i:",["out=", "in="])
   except getopt.GetoptError:
      print 'wpath.py -h'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h' or opt == '--help':
         print 'wpath.py -i --in <inputfile>'
         print 'wpath.py -o --out <outputfile>'
         sys.exit()
      elif opt in ("-o", "--out"):
         filename = arg
      elif opt in ("-i", "--in"):
         fin = arg

   data = open(fin,'r')

   s = data.readlines()
   Nxi = int(s[0])
   Nw = int(s[1])
   Npoints = int(s[2])
   data.close()

   data = np.loadtxt(fin, skiprows=3)

   path_xi = data[:,0]
   path_w = data[:,1]

   path_xi = path_xi.reshape((Nxi,Nw,Npoints))
   path_xi = np.cos(path_xi)
   path_w = path_w.reshape((Nxi, Nw, Npoints))

   xis = []
   for i in range(0,Nxi):
      for j in range(0,Nw):
         xis = np.append(xis,[np.round(path_xi[i,j,Npoints-1],2)])
   
   unique, counts = np.unique(xis,return_counts=True)
   unique = unique[counts>2]
   counts = counts[counts>2]
   counts = 1.0*counts/sum(counts)
   
   np.savetxt(filename,np.transpose([unique,counts]),fmt='%.2f %.4f')
