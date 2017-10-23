# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, getopt
import errno
import numpy as np
import math
import os
from itertools import islice

def make_sure_path_exists(path):
 try:
  os.makedirs(path)
 except OSError as exception:
  if exception.errno != errno.EEXIST:
   raise

if __name__ == "__main__":
	pth = 'out'
	data_in = 'log'
	data_out = 'alignment'
	numlines = 24000
	
	argv = sys.argv[1:]
	try:
		opts, args = getopt.getopt(argv,"hp:o:i:",["path=","out=", "in="])
	except getopt.GetoptError:
		print 'alignment.py -args <>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h' or opt == '--help':
			print 'alignment.py -args <>'
			sys.exit()
		elif opt in ("-p", "--path"):
			pth = arg
		elif opt in ("-o", "--out"):
			data_out = arg
		elif opt in ("-i", "--in"):
			data_in = arg
   
	make_sure_path_exists(pth)
	output_exists = os.path.isfile(data_out)
	
	input_exists = os.path.isfile(pth+'/'+data_in)
	
	if input_exists:
		pass
	else:
		sys.exit("No input file available, quitting...")

	if not(output_exists):
		with open(data_out, "w") as f:
			f.write('acos(wdotq)  which \n')
			
	with open(pth+'/'+data_in) as log:
		Q = np.genfromtxt(islice(log,17,20))
		log.seek(0)
		
		lines = log.readlines()[-numlines:]
		avew = [0.0,0.0,0.0]
		for line in lines:
			data = line.split("|")
			R = np.reshape(map(float,data[8].split()),(3,3))
			QRT = np.dot(Q,np.transpose(R))
			w = map(float,data[3].split())
			if np.linalg.norm(w) != 0:
				w = w/np.linalg.norm(w)
			else:
				w = [0,0,0]
			avew = avew + np.dot(np.transpose(QRT),w)
		
		avew = avew/numlines
		wq1 = np.absolute(np.dot(avew,[1,0,0]))
		wq3 = np.absolute(np.dot(avew,[0,0,1]))
		ave = np.maximum(wq1,wq3)
		if wq1 > wq3:
			b = 1
		else:
			b = 3
			
	with open(data_out,"a") as f:
		f.write(str(ave)+'  ' + str(b) +'\n')
		
