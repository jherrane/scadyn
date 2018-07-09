import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
from matplotlib.patches import Rectangle

def add_arrow(line, position=None, direction='right', size=15, color=None):
   """
   add an arrow to a line.

   line:       Line2D object
   position:   x-position of the arrow. If None, mean of xdata is taken
   direction:  'left' or 'right'
   size:       size of the arrow in fontsize points
   color:      if None, line color is taken.
   """
   if color is None:
      color = line.get_color()

   xdata = line.get_xdata()
   ydata = line.get_ydata()

   if position is None:
      position = xdata.mean()
   # find closest index
   start_ind = np.argmin(np.absolute(xdata - position))
   if direction == 'right':
      end_ind = start_ind + 1
   else:
      end_ind = start_ind - 1

   line.axes.annotate('',
      xytext=(xdata[start_ind], ydata[start_ind]),
      xy=(xdata[end_ind], ydata[end_ind]),
      arrowprops=dict(arrowstyle="->", color=color),
      size=size
   )
 
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
         
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')

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

   fig = plt.figure(figsize=(10,10))

   for i in range(0,Nxi-1):
      line1 = plt.plot(path_xi[i,0,:],path_w[i,0,:],lw=2)[0]
      add_arrow(line1,size=24)
      line2 = plt.plot(path_xi[i,Nw-1,:],path_w[i,Nw-1,:],lw=2)[0]
      add_arrow(line2,size=24)

   someX, someY = 0,0
   currentAxis = plt.gca()
   currentAxis.add_patch(Rectangle((someX - 1, someY - 1.5), 2, 3, edgecolor="grey", facecolor="grey", zorder=4))

   plt.xlabel(r'\cos \xi', fontsize = 40)
   plt.ylabel(r'\omega/\omega_{T}', fontsize = 40)
   plt.yticks(fontsize=32)
   plt.xticks(fontsize=32)
   ax = fig.add_subplot(111)
   ax.set_xlim((-1, 1))
   plt.savefig(filename+'.png')


