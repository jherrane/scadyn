import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

def cart2sph(x,y,z):
	azimuth = np.arctan2(y,x)
	elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return azimuth, elevation, r

def sph2cart(theta,phi,r):
	x = r * np.sin(phi) * np.cos(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(phi)
	return x, y, z

data=pd.read_csv("test.out",
           delim_whitespace=True,
           skipinitialspace=True)

x = np.asarray(data["x"].tolist())
y = np.asarray(data["y"].tolist())
z = np.asarray(data["z"].tolist())

vals = data["N.a_3"].tolist()
vals = vals/max(np.abs(vals))
vals = np.asarray(vals)

th,ph,r = cart2sph(x,y,z)
theta,phi = np.meshgrid(th,ph)

x_vals, x_idx = np.unique(x, return_inverse=True)
y_vals, y_idx = np.unique(y, return_inverse=True)
vals_array = np.empty(x_vals.shape + y_vals.shape)
vals_array.fill(np.nan) # or whatever your desired missing data flag is
vals_array[x_idx, y_idx] = vals
zz = vals_array.T

X,Y,Z = sph2cart(theta,phi,np.abs(zz))

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
    linewidth=0, antialiased=False, alpha=0.5)

plt.show()
