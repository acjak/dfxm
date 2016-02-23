# !/bin/python
"""blah."""

from lib.getedfdata import *
from lib.dfxm import *
from lib.gauss import *
import numpy as np

# import matplotlib
import matplotlib.pylab as plt
# import seaborn as sns

# import scipy.ndimage
# import itertools
#
# import warnings
#
# from mpi4py import MPI
# import time


def makeGvector(bv, xylin, ang):
	gv = np.ones((3, len(xylin), len(xylin)))  # G-vector

	xcpos = len(gv[0, :, 0])/2
	ycpos = len(gv[0, 0, :])/2

	prefix = 2*math.pi/bv
	# prefix = 2*math.sin(ang)

	for i in range(len(xylin)):
		for j in range(len(xylin)):
			if i != xcpos and j != ycpos:
				gv[0, i, j] = prefix*bv/(2*math.pi)*(xylin[j])/((xylin[i])**2 + (xylin[j])**2)
				gv[1, i, j] = prefix*bv/(2*math.pi)*(xylin[i])/((xylin[i])**2 + (xylin[j])**2)
				gv[2, i, j] = prefix

	return gv


def makeGvectorDiamond(bv_111, bv_110, xylin, ang, phi, w):
	gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

	xcpos = len(gv[0, :, 0])/2
	ycpos = len(gv[0, 0, :])/2

	prefix = 2*math.pi/bv_111
	# prefix = 2*math.sin(ang)

	for i in range(len(xylin)):
		for j in range(len(xylin)):
			if i != xcpos or j != ycpos:
				eps_xz = (bv_110/(2*math.pi))*(xylin[j])/((xylin[i])**2 + (xylin[j])**2)
				eps_yz = (-bv_110/(2*math.pi))*(xylin[i])/((xylin[i])**2 + (xylin[j])**2)
				gv[0, i, j] = prefix*(eps_xz*math.cos(w)-eps_yz*math.sin(phi)*math.sin(w))
				gv[1, i, j] = prefix*(eps_yz*math.cos(phi))
				gv[2, i, i] = prefix*(eps_xz*math.cos(w)*math.sin(phi)+eps_xz*math.sin(w))

	return gv


def makeGvectorDiamondPantleon(bv_111, bv_110, xylin, ang, phi, alpha, beta):
	gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

	xcpos = len(gv[0, :, 0])/2
	ycpos = len(gv[0, 0, :])/2

	prefix = 2*math.pi/bv_111
	# prefix = 1
	# prefix = 2*math.sin(ang)

	sa = math.sin(alpha)
	ca = math.cos(alpha)
	print bv_111, bv_110

	for i in range(len(xylin)):
		for j in range(len(xylin)):
			if i != xcpos or j != ycpos:
				xc = -xylin[i]*sa
				eps_xz = (-bv_110/(2*math.pi))*(xc)/((xylin[j])**2 + (xc)**2)
				eps_yz = (bv_110/(2*math.pi))*(xylin[j])/((xylin[j])**2 + (xc)**2)
				gv[0, i, j] = prefix*(eps_xz)  # prefix*(eps_yz*cw-sw)
				gv[1, i, j] = prefix*(eps_yz*ca-sa)  # prefix*(eps_xz)
				gv[2, i, j] = prefix*(eps_yz*ca+sa)

	return gv


def makeGvectorDiamondSonja(bv_111, bv_110, xylin, ang, alpha, beta):
	gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

	xcpos = len(gv[0, :, 0])/2
	ycpos = len(gv[0, 0, :])/2

	sa = math.sin(math.radians(alpha))
	ca = math.cos(math.radians(alpha))
	sb = math.sin(math.radians(beta))
	cb = math.cos(math.radians(beta))

	z = 0

	prefix = (bv_110/bv_111)*ca

	for i in range(len(xylin)):
		for j in range(len(xylin)):
			if i != xcpos or j != ycpos:

				r2p = (ca*(xylin[i]*cb+xylin[j]*sb)-sa*z)**2 + (-xylin[i]*sb+xylin[j]*cb)**2
				gv[0, i, j] = prefix*(1/r2p)*(-sa*sb*z + ca*xylin[j])
				gv[1, i, j] = prefix*(1/r2p)*(sa*cb*z - ca*xylin[i])
				gv[2, i, j] = prefix*(1/r2p)*(sa*(sb*xylin[i] + cb*xylin[j]))+2*math.pi/bv_111

	return gv


def makeQvector(wl, omega, xi):
	q = [(2*math.pi/wl)*2*np.sin(xi/2)*np.sin(omega-xi/2),
		(2*math.pi/wl)*2*np.sin(xi/2)*np.sin(omega-xi/2),
		(2*math.pi/wl)*2*np.sin(xi/2)*np.cos(omega-xi/2)]

	return q


def fullArray(gv, wl, xylin, omega0, alpha, beta):
	# linearray = np.zeros((len(data.alphavals),len(data.betavals)))

	a, b, c = data.getMetaValues()

	dist = np.zeros((3))

	sigmax = 1E-5
	sigmay = 5E-4
	sigmaz = 1E-4

	om_start = 10
	xi_start = 15

	norm_max = 2762.119
	# Projection coordinates
	# x0, y0, x1, y1
	y0 = 0
	pc = [20, y0, 20, y0+len(xylin)]

	q0 = makeQvector(wl, np.radians(omega0), np.radians(2*omega0))
	ql0 = math.sqrt(q0[0]**2 + q0[1]**2 + q0[2]**2)

	with open(data.infFile,  "a") as myfile:
		myfile.write("sigmax: %s \n" % (str(sigmax)))
		myfile.write("sigmay: %s \n" % (str(sigmay)))
		myfile.write("sigmaz: %s \n" % (str(sigmaz)))
		myfile.write("\n")
		myfile.write("om_start: %s \n" % (str(om_start)))
		myfile.write("xi_start: %s \n" % (str(xi_start)))
		myfile.write("\n\n\n")

	k = 0

	fig = plt.figure(figsize=(18, 4))
	gs1 = matplotlib.gridspec.GridSpec(1, 10)
	gs1.update(wspace=.2,  hspace=0.2)

	# sns.set_style("white")

	for i in range(1):  # len(a)/2):
		for j in range(10):
			# Empty array for model calculations.
			matrix_part = np.zeros((len(xylin), len(xylin)))

			# Omega and xi defined from the peak position in degrees.
			omega = b[j+om_start]-data.beta0+0.0002  # +0.00085
			xi = a[i+xi_start]-data.alpha0

			print data.alpha0, data.beta0, omega0
			print omega, xi

			# Import image from data set.
			ind = data.getIndex(float(a[i+xi_start]), float(b[j+om_start]), float(beta))[0]
			img = data.getImage(ind, False)

			if j == 2:
				im = img

			# Creat subplot in plot array.
			axarr = plt.subplot(gs1[k])
			axarr.locator_params(nbins=5, axis='x')
			# axarr.set_title('%.3f %.3f' % (omega, xi),fontsize=8)
			axarr.set_title('%.4f %.4f' % (b[j+om_start], a[i+xi_start]), fontsize=8)

			# Get an average across the dislocation line.
			line = np.sum(img, axis=1)/len(img[0, :])
			# line = np.sum(img[0:2,0], axis=1)/len(img[0, :])
			# line = img[:, 0]

			# Make Q vector from the same angle as data.
			q = makeQvector(wl, np.radians(omega0+omega), np.radians(2*omega0+xi))

			# Length of q.
			ql = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)
			dq = ql0-ql

			# sigma_a = 1.3e-5
			sigma_a = 7.698e-5
			# sigma_a = 1.e-3
			# sigma_a = 1.

			sigmax = 2*ql*sigma_a  # 2*ql*sigma_a
			sigmay = ql*1e-5  # ql*1e-5
			sigmaz = ql*(1/math.tan(np.radians(omega0+omega)))*sigma_a

			sigma = [sigmax, sigmay, sigmaz]

			matrix_part[:, :] = 1.
			for m in range(len(q)):
				# dist[m] = abs(gv[m, :, :] - q[m])
				matrix_part[:, :] *= np.exp(-(abs(gv[m, :, :] - q[m])**2)/(2*sigma[m]**2))

			matrix_part[:, :] *= 2000.

			if dq != 0:
				matrix_part[:, :] += 2.e-12*1./((dq-5.5e-8)**2)

			proj_line = tools.getProjection(img, pc[0], pc[1], pc[2], pc[3], len(xylin))
			mean_line = np.sum(img, axis=1)

			axarr.plot(xylin, mean_line)
			axarr.plot([10, 10], [0, 1.1*np.max(mean_line)])

			# Make a line from the model and plot it on the same subplot.
			model_line = np.sum(matrix_part[:, :], axis=1)[::-1]

			maxpoint = np.argmax(model_line)
			maxval = max(model_line)
			for f in range(len(model_line)):
				distpixel = (maxpoint - f)
				fraction = abs((distpixel*pixelsize)/4e9)
				# model_line[f] = 0.3*maxval*np.exp(-(float(abs(distpixel))/(2*10.**2)))
				print f, float(abs(distpixel)), maxval*np.exp(-(float(abs(distpixel))/(2*2.**2)))
				# if fraction != 0:
				# 	theta = math.atan(fraction)
				# 	qpixel = (4*math.pi/wl)*math.sin(theta)
				# 	print math.degrees(theta), qpixel, fraction, distpixel
				# 	model_line[f] += (1/(qpixel**(2)))/10000000000

			axarr.plot(xylin, model_line)
			axarr.set_xticklabels(map(str, map(int,axarr.get_xticks()/1000.0)))
			axarr.set_xlabel('um')

			# Loop over all pixels in the cross-sectional model image.
			for e in range(len(xylin)):
				for d in range(len(xylin)):

					# Find the distance from Q to G vector in x, y and z.
					for m in range(len(q)):
						# dist[m] = abs(gv[m, e, d] - q[m])
						# matrix_part[e, d] *= math.exp(-(dist[m]**2)/(2*sigma[m]**2))
						if e == 5 and d == 3:
							print m, dist[m], gv[m, e, d], q[m], math.exp(-(dist[m]**2)/(2*sigma[m]**2))
							if m == 2:
								print ""

					# Use the distance to fold with a Gaussian.
					# matrix_part[e, d] = 1.
					# matrix_part[e, d] *= math.exp(-(dist[0]**2)/(2*sigmax**2))
					# matrix_part[e, d] *= math.exp(-(dist[1]**2)/(2*sigmay**2))
					# matrix_part[e, d] *= math.exp(-(dist[2]**2)/(2*sigmaz**2))
					# matrix_part[e, d] *= 30.
					# if dq != 0:
					# 	matrix_part[e, d] += 1.e-14*1./((dq-5.5e-8)**2)

			# Plot data.
			# axarr.plot(xylin, line)
			# axarr.plot([xylin[np.argmax(line)], xylin[np.argmax(line)]], [0, np.max(line)])
			# axarr.plot([40, 40], [0, np.max(line)])

			# proj_line = tools.getProjection(img, 20, 80, 80, 20, 100)
			# proj_array = np.zeros((len(xylin)))
			# for f in range(3):
			# 	proj_array += tools.getProjection(img, pc[0], pc[1], pc[2], pc[3], len(xylin))
			#
			# proj_line = proj_array/3

			# fig2, ax2 = plt.subplots()
			# ax2.imshow(img, cmap="Greens", interpolation=None)
			# fig2.savefig(data.directory + '/disl_im_%g.png' % (j))
			# fig2.clf()

			# axarr.plot(xylin, np.max(mean_line)*model_line/max(model_line))

			# axarr.set_ylim(0, 600)
			# axarr.imshow(img, cmap='Greens')

			# axarr.plot([np.argmax(proj_line), np.argmax(proj_line)], [0, np.max(proj_line)])
			# axarr.plot([250, 250], [0, np.max(line)])

			# axarr.imshow(matrix_part, cmap='Greens')
			# axarr.set_ylim(0,2000)
			# axarr.xaxis.set_major_formatter(plt.NullFormatter())
			# axarr.yaxis.set_major_formatter(plt.NullFormatter())

			k += 1

	fig.savefig(data.directory + '/comp_screwdisl.pdf')
	fig3, ax3 = plt.subplots()
	ax3.imshow(im, cmap='Greens', interpolation=None)
	ax3.autoscale(False)
	ax3.plot([pc[0], pc[2]], [pc[1], pc[3]])
	fig3.savefig(data.directory + '/image_with_line.png')
	plt.show()


path = '/data/hxrm/Dislocation_may_2015/dislocations/strain'

path = '/data/hxrm/Dislocation_november_2015/diamond/ff_strain'
bg_path = '/data/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff2_'
# filename2 = 'ff2_'
sampletitle = 'topo_strain'
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'strain_tt'

test_switch = True

# poi = [750, 750]
# size = [800, 200]

# poi = [300, 650]  # [300, 650]
# size = [52, 430]
poi = [277, 417]
# size = [2, 20]
size = [50, 40]

pixelsize = 180.  # 180.

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)

try:
	directory = data.directory
except AttributeError:
	directory = 0
tools = DFXM(path, data.data_files, directory, roi, datatype, data.dirhash, data.meta, test_switch)
# data.printMeta()
meta = data.getMetaArray()
# hist, datx, daty = data.makeMeanGrid()
a, b, c = data.getMetaValues()

data.setTest(True)
data.adjustOffset(False)
# alignToCenter(meta, a, b, hist)
# overlapImg(meta, a, b, hist, [10, 26], [11, 26])
# compareTheta(meta, eta, cols)

# omega0 = 10.94822
# ang = np.radians(omega0)
#
# bv = 0.192  # Burger's vector

wl = 0.07293  # Wavelength

bv_111 = 0.35668/math.sqrt(3)   # Burger's vector
bv_110 = 0.35668/math.sqrt(2)   # Burger's vector
diffrx = 150.  # +0.
alpha_disl = -32.5  # -32.5
omega0 = math.degrees(math.asin(wl/(2*bv_111)))
ang = math.radians(omega0)

ln = size[1]*pixelsize

xylin = np.ogrid[-ln/2:ln/2:pixelsize]

gv = makeGvectorDiamondSonja(bv_111, bv_110, xylin, ang, alpha_disl, diffrx)
# print xylin
# print len(xylin)
fullArray(gv, wl, xylin, omega0, alpha_disl, diffrx)

# bins, length, xr = diagonalPlot(meta, a, b, hist, 1, 0.53, datatype)
