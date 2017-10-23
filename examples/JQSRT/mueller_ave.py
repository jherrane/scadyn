# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, getopt
import errno
import numpy as np
import math
import os

class cd:
 """Context manager for changing the current working directory"""
 def __init__(self, newPath):
  self.newPath = os.path.expanduser(newPath)

 def __enter__(self):
  self.savedPath = os.getcwd()
  os.chdir(self.newPath)

 def __exit__(self, etype, value, traceback):
  os.chdir(self.savedPath)

def make_sure_path_exists(path):
 try:
  os.makedirs(path)
 except OSError as exception:
  if exception.errno != errno.EEXIST:
   raise

if __name__ == "__main__":
	pth = 'mueller'
	mueller_in = 'mueller'
	mueller_out = 'mueller_all'
	
	argv = sys.argv[1:]
	try:
		opts, args = getopt.getopt(argv,"hp:o:i:",["path=","out=", "in="])
	except getopt.GetoptError:
		print 'mueller_all.py -args <>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h' or opt == '--help':
			print 'mueller_all.py -args <>'
			sys.exit()
		elif opt in ("-p", "--path"):
			pth = arg
		elif opt in ("-o", "--out"):
			mueller_out = arg
		elif opt in ("-i", "--in"):
			mueller_in = arg
   
	make_sure_path_exists(pth)
	output_exists = os.path.isfile(mueller_out)
	
	input_exists = os.path.isfile(pth+'/'+mueller_in)
	
	if input_exists:
		pass
	else:
		sys.exit("No input file available, quitting... ("+pth+'/'+mueller_in+')')

	if not(output_exists):
		with open(pth+'/'+mueller_in) as f:
			lines = f.readlines()
			with open(mueller_out, "w") as f1:
				f1.writelines(lines)
	else:
		data1 = np.zeros((180,18))
		data2 = np.zeros((180,18))
		data3 = np.zeros((180,18))
		with open(pth+'/'+mueller_in) as f:
			lines = f.readlines()[-180:]
			i = 0
			for line in lines:
				floatline = [float(x) for x in line.split()]
				data1[i][:] = floatline
				i += 1
		with open(mueller_out) as f:
			lines = f.readlines()[-180:]
			lines2 = open(mueller_out).read().splitlines()
			i = 0
			for line in lines:
				floatline = [float(x) for x in line.split()]
				data2[i][:] = floatline
				data3[i][0:2] = data1[i][0:2]
				data3[i][2:] = data1[i][2:]+data2[i][2:]
				string = ["%.4e" % x for x in data3[i][:]]
				lines2[i+1] = "  ".join(string)
				i += 1
			open(mueller_out,'w').write('\n'.join(lines2))
		
		
		
