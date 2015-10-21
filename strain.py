# !/bin/python
"""blah."""

from lib.getedfdata import *
from lib.gauss import *
import numpy as np

# import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

import scipy.ndimage
import itertools

import warnings

from mpi4py import MPI
import time


def compareTheta(data_or, eta, cols):
	first = np.where(data_or[:, 0] == eta)[0][12]
	second = np.where(data_or[:, 0] == eta)[0][34]

	fig, axarr = plt.subplots(nrows=2, ncols=cols, figsize=(16, 8))
	index = np.zeros((length, 2), dtype=np.int16)
	for i in range(cols):
		index[i, 0] = first + i
		index[i, 0] = second - i

	data.makePlotArray(index, 'comparetheta_%g' % eta)

	plt.show()


def diagonalPlot(meta, alpha, beta, hist, bins, xpos, datatype, rank):
	hist_array = np.zeros((len(hist[:, 0]), len(hist[0, :]), 1))

	# hist_array[:, :, 0] = hist

	startpos = [3, 6]
	startpos2 = [6, 34]
	length = 13
	show_index = 2
	index = np.zeros((length, 1), dtype=np.int16)

	i1 = np.zeros((2, 2, length))
	xr = []

	for i in range(length):
		i1[0, 0, i] = startpos[0]+1*i
		i1[0, 1, i] = startpos[1]+2.5*i
		hist[i1[0, 0, i], i1[0, 1, i]] = 3000
		i1[1, 0, i] = startpos2[0]+2*i
		i1[1, 1, i] = startpos2[1]+10*i
		# hist[i1[1, 0, i], i1[1, 1, i]] = 3000
		# print float(alpha[i1]), float(beta[i2])
		index[i, 0] = data.getIndex(float(alpha[i1[0, 0, i]]), float(beta[i1[0, 1, i]]))[0]
		# index[i, 1] = data.getIndex(float(alpha[i1[1, 0, i]]), float(beta[i1[1, 1, i]]))[0]
		xr.append(beta[i1[0, 1, i]] - 10.992-7.86E-4)

	data.makePlotArray(index, bins, xpos, 'diagplot_%g_%g' % (startpos[0], startpos[1]))

	# fig2, ax2 = plt.subplots()
	# for i in range(len(data.imgarray[:, 0, 0])):
	# 	ax2.imshow(data.imgarray[i, :, :])
	# 	fig2.savefig('output/diag_oldfilter_compare_%g.png' % i)
		# fig2.clf()

	if rank == 0:
		sns.set_style("white")
		hist_array[:, :, 0] = hist
		fig,  ax = data.makeHistogram(hist_array, alpha, beta, 'diagplot_hist_%g_%g' % (startpos[0], startpos[1]))
		# data.showArea(i1[0, 0, show_index], i1[0, 1, show_index])
		# ax.scatter(i1[0, 1, show_index]+0.5, i1[0, 0, show_index]+0.5, s=100, color='magenta')
	return bins, length, xr


def horizontalPlot(meta, alpha, beta, hist, bins, xpos, datatype, rank):
	hist_array = np.zeros((len(hist[:, 0]), len(hist[0, :]), 1))

	startpos = [10, 0]
	# startpos2 = [6, 34]
	length = 40
	show_index = 2
	index = np.zeros((length, 1), dtype=np.int16)

	i1 = np.zeros((2, 2, length))
	xr = []

	for i in range(length):
		i1[0, 0, i] = startpos[0]
		i1[0, 1, i] = startpos[1]+1*i
		hist[i1[0, 0, i], i1[0, 1, i]] = 3000
		# i1[1, 0, i] = startpos2[0]+2*i
		# i1[1, 1, i] = startpos2[1]+10*i
		# hist[i1[1, 0, i], i1[1, 1, i]] = 3000
		# print float(alpha[i1]), float(beta[i2])
		index[i, 0] = data.getIndex(float(alpha[i1[0, 0, i]]), float(beta[i1[0, 1, i]]))[0]
		# index[i, 1] = data.getIndex(float(alpha[i1[1, 0, i]]), float(beta[i1[1, 1, i]]))[0]
		xr.append(beta[i1[0, 1, i]] - 10.992-7.86E-4)

	data.makePlotArray(index, bins, xpos, 'horizplot_%g_%g' % (startpos[0], startpos[1]))

	# fig2, ax2 = plt.subplots()
	# for i in range(len(data.imgArray[:, 0, 0])):
	# 	ax2.imshow(data.imgArray[i, :, :])
	# 	fig2.savefig('output/horizontal/horiz_compare_%g.png' % i)

	if rank == 0:
		hist_array[:, :, 0] = hist
		fig,  ax = data.makeHistogram(hist_array, alpha, beta, 'horizplot_%g_%g' % (startpos[0], startpos[1]))
		data.showArea(i1[0, 0, 1], i1[0, 1, 1])
		ax.scatter(i1[0, 1, show_index]+0.5, i1[0, 0, show_index]+0.5, s=100, color='magenta')
	return bins, length, xr



def overlapImg(meta, a, b, hist, pos1, pos2):
	pos1 = [10, 20]
	index = data.getIndex(float(data.alphavals[pos1[0]]), float(data.betavals[pos1[1]]))

	img = data.getImage(index[0], False)

	# pos2 = [pos1[0]-1, pos1[1]]

	index = data.getIndex(float(data.alphavals[pos2[0]]), float(data.betavals[pos2[1]]))
	img1 = data.getImage(index[0], False)

	bestsum = 1E8
	guessrange = range(-50, 50)
	for i in guessrange:
		for j in guessrange:

			if i != 0 and j != 0:
				newimg = img[i:, j:] - img1[:-i, :-j]

			if i == 0 and j != 0:
				newimg = img[:, j:] - img1[:, :-j]

			if i != 0 and j == 0:
				newimg = img[i:, :] - img1[:-i, :]

			if i == 0 and j == 0:
				newimg = img - img1

			newsum = abs(np.sum(newimg))
			if newsum < bestsum:
				bestsum = newsum
				bestindex = [i, j]

	newimg = img[bestindex[0]:, bestindex[1]:] - img1[:-bestindex[0], :-bestindex[1]]
	print bestsum,  bestindex

	return bestindex

def makeGrid(meta, a, b):

	sns.set_style("white")
	sns.set_context("paper")

	# plt.figure(figsize = (14, 14))
	# gs1 = matplotlib.gridspec.GridSpec(4, 4)
	# gs1.update(wspace=0.025,  hspace=0.03)

	# fig, ax = plt.subplot()

	avals = a[5:15]
	bvals = b[15:16]  # b[15:25:3]
	print len(a)
	print avals
	print bvals
	ab_vals = list(itertools.product(avals, bvals))

	print len(ab_vals)

	for i in range(len(ab_vals)):
		axarr = plt.subplot()  # (gs1[i])
		# print ab_vals[i][1]
		# print ab_vals[i][0]
		index = data.getIndex(float(ab_vals[i][0]), float(ab_vals[i][1]))
		# print index
		img = data.getImage(index[0], False)
		# ta = np.ones((len(img[:, 0]), len(img[0, :]), 4),  dtype=np.uint8)*0

		# ta[:, :, 3] = 255
		# ta[:, :, 1] = img# 255*img1/np.max(img1)

		axarr.imshow(img)
		axarr.autoscale(True)
		axarr.set_title('%.4f %.4f' % (meta[index[0], 1], meta[index[0], 0]))
		plt.savefig('output/%s_%s.pdf' % ('straintest2', i))
		# plt.show()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
	start = time.time()


path = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/strain'

filename = 'strainmap_tt_2'
sampletitle = filename
bg_filename = 'bg1_5s_'

datatype = 'strain_tt'


#poi = [750, 750]
# size = [800, 200]

poi = [1250, 1150]
# size = [100, 100]
size = [600, 300]
# poi = [1000, 1000]
# size = [1000, 1000]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

eta = -0.01
cols = 8

data = GetEdfData(path, filename, bg_filename, roi, datatype)
# data.printMeta()
meta = data.getMetaArray()
hist, datx, daty = data.makeMeanGrid()
a, b = data.getMetaValues()

data.setTest(True)
data.adjustOffset(True)
# alignToCenter(meta, a, b, hist)
# overlapImg(meta, a, b, hist, [10, 26], [11, 26])
# compareTheta(meta, eta, cols)

bins, length, xr = diagonalPlot(meta, a, b, hist, 1, 0.53, datatype, rank)
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	strainpic = data.makeStrainArrayMPI(data.imgarray, bins, length, xr)


# bins, length, xr = horizontalPlot(meta, a, b, hist, 1, 0.53, datatype, rank)
# with warnings.catch_warnings():
# 	warnings.simplefilter("ignore")
# 	gaussarray = data.makeGaussArrayMPI(data.imgarray, bins, length, xr)

# strainpic = makeStrainPlot(bins, length, xr)
if rank == 0:
	end = time.time()
	print "Time:", end-start
	data.plotStrain(strainpic)
	# data.plotPPOS(gaussarray, length)
	# data.plotStrainHeatmap(data.imgarray)
	plt.show()
# HorizCompare(meta, a, b, hist, datatype)
# makeGrid(meta, a, b)
# verticalPlot(meta, a, b, hist, 1, 0.3, datatype)
