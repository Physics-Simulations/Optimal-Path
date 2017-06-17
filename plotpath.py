"""
	Utility script to print the paths found in a 3D graph with the potential.
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

COLORS = ['g', 'm', 'b', 'r', 'y']

def plot_path(paths, P, r0, rf, xlim, ylim, labels):
	"""
		This will show a graph with the paths and the potential.

		paths : an array of paths containing the points 

		P : potential as a function of x and y

		r0 : initial point

		rf : final point

		xlim : dimensions of the x axis as tuple (xmin, xmax)

		ylim : dimensions of the y axis as tuple (ymin, ymax)

		labels : a description for each path
	"""
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	X = np.arange(xlim[0], xlim[1], 0.1)
	Y = np.arange(ylim[0], ylim[1], 0.1)
	X, Y = np.meshgrid(X, Y)
	Z = P(X, Y)

	# Plot the surface.
	surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
	                       linewidth=0, antialiased=False)

	for k in xrange(len(paths)):
		ax.plot(paths[k][:,0], paths[k][:,1], P(paths[k][:,0], paths[k][:,1]), 
			ls='-', lw=4.0, color=COLORS[k % 5], label=labels[k])

	ax.plot([r0[0]], [r0[1]], [P(r0[0], r0[1])], marker='^', color='w', ls='None', ms=10)
	ax.plot([rf[0]], [rf[1]], [P(rf[0], rf[1])], marker='d', color='k', ls='None', ms=10)

	# Customize the z axis.
	#ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_zlabel('$U(x,y)$')
	ax.legend()

	# Add a color bar which maps values to colors.
	fig.colorbar(surf, shrink=0.5, aspect=8)
	plt.subplots_adjust(left=0, right=1.0, top=1.0, bottom=0)
	plt.show()