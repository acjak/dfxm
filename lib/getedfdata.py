# !/bin/python
"""Class for loading DFXM data sets.

The class can be loaded from another Python file. This gives access to all metadata in the data set as well as direct access to an image by giving either coordinates or an index.

To use the class in another Python file:

	from lib.getedfdata import *

where getedfdata.py is in a sub directory called 'lib'.

An example of using the class can be found in the test function at the bottom of this file.
"""
import os
import math
import numpy as np
import scipy
import EdfFile
from os import listdir
from os.path import isfile,  join
import cv2
from time import localtime,  strftime
import matplotlib
import matplotlib.pylab as plt
import seaborn as sns
from mpi4py import MPI
import time


class GetEdfData(object):

	"""docstring for GetEdfData."""

	def __init__(self,  path, filename, bg_filename, roi, datatype):
		super(GetEdfData,  self).__init__()

		self.comm = MPI.COMM_WORLD
		self.rank = self.comm.Get_rank()
		self.size = self.comm.Get_size()

		self.datatype = datatype
		self.sampletitle = filename
		self.path = path
		self.roi = roi
		self.getFilelists(filename, bg_filename)
		self.getBGarray()
		self.getMetaData()
		self.makeROIAdjustmentArray()
		if self.rank == 0:
			self.makeOutputFolder()
			self.printInfoToFile()

	def makeOutputFolder(self):
		timestamp = strftime("%d%m%y-%H%M",  localtime())
		self.directory = 'output/' + timestamp
		print self.directory

		if not os.path.exists(self.directory):
			os.makedirs(self.directory)

	def printInfoToFile(self):
		self.infFile = self.directory + "/inf.txt"
		with open(self.infFile,  "a") as myfile:
			myfile.write("Datatype: %s \nSampletitle: %s \nPath: %s \n" % (self.datatype, self.sampletitle, self.path))
			myfile.write("ROI: %s \n" % (str(self.roi)))

			myfile.write("\n\n\n")

	def setTest(self, testcase):
		self.test = testcase

	def adjustOffset(self, case):
		self.adjustoffset = case

	def getFilelists(self,  filename,  bg_filename):
		onlyfiles = [f for f in listdir(self.path) if isfile(join(self.path,  f))]

		self.data_files = []
		self.bg_files = []

		for k in onlyfiles[:]:
			if k[:len(filename)] == filename:
				self.data_files.append(k)

			if k[:len(bg_filename)] == bg_filename:
				self.bg_files.append(k)

	def getROI(self):
		return self.roi

	def setROI(self,  roi):
		self.roi = roi

	def makeROIAdjustmentArray(self):
		self.adj_array = np.zeros((len(self.alphavals),  4))
		for i in range(len(self.alphavals)):
			offset = (i-10)*2  # *4
			# print i,  offset
			self.adj_array[i,  0] = self.roi[0]
			self.adj_array[i,  1] = self.roi[1]
			self.adj_array[i,  2] = self.roi[2]+offset
			self.adj_array[i,  3] = self.roi[3]+offset

	def adjustROI(self,  xoff,  yoff):
		for i in range(len(self.alphavals)):
			self.adj_array[i,  0] = self.adj_array[i,  0]+xoff
			self.adj_array[i,  1] = self.adj_array[i,  1]+xoff
			self.adj_array[i,  2] = self.adj_array[i,  2]+yoff
			self.adj_array[i,  3] = self.adj_array[i,  3]+yoff

	def getBGarray(self):
		bg_file_with_path = self.path + '/' + self.bg_files[0]
		bg_class = EdfFile.EdfFile(bg_file_with_path)
		bg_img = bg_class.GetData(0).astype(np.int64)[self.roi[2]:self.roi[3],  self.roi[0]:self.roi[1]]

		self.bg_combined = np.zeros(np.shape(bg_img))

		if self.rank == 0:
			print "Reading background files (ROI)..."

		for i in range(len(self.bg_files)):
			bg_file_with_path = self.path + '/' + self.bg_files[i]
			bg_class = EdfFile.EdfFile(bg_file_with_path)
			self.bg_combined += bg_class.GetData(0).astype(np.int64)[self.roi[2]:self.roi[3],  self.roi[0]:self.roi[1]]

		self.bg_combined /= len(self.bg_files)

		bg_img_full = bg_class.GetData(0).astype(np.int64)
		self.bg_combined_full = np.zeros(np.shape(bg_img_full))

		if self.rank == 0:
			print "Reading backg files (Full)..."

		for i in range(len(self.bg_files)):
			bg_file_with_path = self.path + '/' + self.bg_files[i]
			bg_class = EdfFile.EdfFile(bg_file_with_path)
			self.bg_combined_full += bg_class.GetData(0).astype(np.int64)

		self.bg_combined_full /= len(self.bg_files)

	def getMeanData(self):
		meandatafile = 'datamean_%s.txt' % self.sampletitle

		if os.path.isfile(meandatafile) == True:
			self.data_mean = np.loadtxt(meandatafile)
		else:
			self.data_mean = np.zeros((len(self.data_files)))
			if self.rank == 0:
				print "Reading %g data files..." % len(self.data_files)
			for i in range(len(self.data_files)):
				file_with_path = self.path + '/' + self.data_files[i]
				img = EdfFile.EdfFile(file_with_path)
				self.data_mean[i] = img.GetData(0).astype(np.int64)[self.roi[2]:self.roi[3],  self.roi[0]:self.roi[1]].mean()-self.bg_combined.mean()
			np.savetxt('datamean_%s.txt' % self.sampletitle,  self.data_mean)

	def getMean(self):
		self.getMeanData()
		return self.data_mean

	def getMetaData(self):
		self.meta = np.zeros((len(self.data_files),  4))

		def getHeader(filenumber):
			file_with_path = self.path + '/' + self.data_files[i]
			img = EdfFile.EdfFile(file_with_path)
			header = img.GetHeader(0)

			mot_array = header['motor_mne'].split(' ')
			motpos_array = header['motor_pos'].split(' ')

			return mot_array, motpos_array

		if self.rank == 0:
			print "Reading meta data..."

		for i in range(len(self.data_files)):
			mot_array, motpos_array = getHeader(i)

			if self.datatype == 'topotomo':
				self.meta[i,  0] = round(float(motpos_array[mot_array.index('diffrx')]),  8)

			if self.datatype == 'strain_eta':
				self.meta[i,  0] = round(float(motpos_array[mot_array.index('obpitch')]),  8)

			if self.datatype == 'strain_tt':
				self.meta[i,  0] = round(float(motpos_array[mot_array.index('obyaw')]),  8)

			sn = float(self.data_files[i][-8:-4])
			theta = (11.01-10.97)/40
			self.meta[i,  1] = round(10.9700+theta*sn+theta/2,  8)
			self.meta[i,  2] = sn

		alphavals = sorted(list(set(self.meta[:,  0])))
		betavals = sorted(list(set(self.meta[:,  1])))
		self.alphavals = np.zeros((len(alphavals)))
		self.betavals = np.zeros((len(betavals)))
		for i in range(len(alphavals)):
			self.alphavals[i] = float(alphavals[i])
		for i in range(len(betavals)):
			self.betavals[i] = float(betavals[i])
		if self.rank == 0:
			print "Meta data from %s files read." % str(len(self.data_files))

	def removeHotPixels(self,  data_new):
		width = 15
		hp = [[130,  784], [750,  1390]]

		for i in range(len(hp)):
			if hp[i][0] <= len(data_new[:, 0])-width and hp[i][1] <= len(data_new[0, :])-width:
				data_new[hp[i][1]-width:hp[i][1]+width, hp[i][0]-width:hp[i][0]+width] = 0
		return data_new

	def rebin(self, a, bs):
		shape = (a.shape[0]/bs, a.shape[1]/bs)
		sh = shape[0], a.shape[0]//shape[0], shape[1], a.shape[1]//shape[1]
		return a.reshape(sh).sum(-1).sum(1)

	def gaus(self, x, a, x0, sigma):
		return a*np.exp(-(x-x0)**2/(2*sigma**2))

	def fitGaussian(self, x, y):
		from scipy.optimize import curve_fit
		# n = len(x)                          # the number of data
		# mean = sum(x*y)/n                   # note this correction
		# sigma = sum(y*(x-mean)**2)/n        # note this correction
		# sigma = -3.E-3
		try:
			popt, pcov = curve_fit(self.gaus, x, y, p0=[max(y), x[np.argmax(y)], -3.E-3], maxfev=100000)
			return popt, pcov
		except RuntimeError:
			print "Error - curve_fit failed"

	def getIndex(self, alpha, beta):
		if alpha != -1 and beta == -1:
			index = np.where(self.meta[:, 0] == alpha)
		if alpha == -1 and beta != -1:
			index = np.where(self.meta[:, 1] == beta)
		if alpha != -1 and beta != -1:
			i1 = np.where(self.meta[:, 0] == alpha)
			i2 = np.where(self.meta[:, 1] == beta)
			index = list(set(i1[0]).intersection(i2[0]))
		return index

	def getImage(self, index, full):
		tmpfile = 'output/tmp/img_' + str(index) + '.npy'
		file_with_path = self.path + '/' + self.data_files[index]
		if self.rank == 0:
			print file_with_path
		# print index

		# if os.path.isfile(tmpfile):
		# 	im = np.load(tmpfile)
		if True:
			img = EdfFile.EdfFile(file_with_path)
			if self.adjustoffset is True:
				alpha = self.meta[index, 0]
				# beta = self.meta[index, 1]
				# print alpha
				a_index = np.where(self.alphavals == alpha)
				# b_index = np.where(self.betavals == beta)
				# print b_index[0]
				roi = self.adj_array[a_index[0]][0]
				# print roi
			else:
				print roi
				roi = self.roi
			if full is True:
				im = img.GetData(0).astype(np.int64)-self.bg_combined_full
			else:
				im = img.GetData(0).astype(np.int64)[roi[2]:roi[3], roi[0]:roi[1]]-self.bg_combined
			im = self.cleanImage(im)
			np.save(tmpfile,im)
		return im

	def cleanImage(self, img):

		# kernel = np.ones((5, 5), np.float32)/25
		# img = cv2.filter2D(img, -1, kernel)

		img = self.rfilter(img,18,3)
		img[img < 0] = 0

		return img

	def smooth(self, a, n):
		"""Do a moving average."""
		ret = np.cumsum(a, dtype=float)
		ret[n:] = ret[n:] - ret[:-n]
		return ret / n

	def rfilter(self,img,nstp,slen):
		mask = np.ones(np.shape(img))
		stp = 180./nstp

		stack = np.zeros((len(img[:,0]),len(img[0,:]),nstp))

		for n in range(nstp):
			rot = n*stp
			imgr = scipy.ndimage.interpolation.rotate(img, rot)
			mskr = scipy.ndimage.interpolation.rotate(mask, rot)

			for j in range(len(mskr[:,0])):
				ids = np.nonzero(mskr[j,:])
				if bool(ids[0].any()):
					imgr[j,ids[0]] = self.smooth(imgr[j,ids[0]],slen)

			imgb = scipy.ndimage.interpolation.rotate(imgr, -rot)

			idx0 = [np.floor((len(imgb[:,0])-len(img[:,0]))/2),np.floor((len(imgb[0,:])-len(img[0,:]))/2)]
			idx1 = [np.floor((len(imgb[:,0])+len(img[:,0]))/2),np.floor((len(imgb[0,:])+len(img[0,:]))/2)]

			imgb = imgb[idx0[0]:idx1[0],idx0[1]:idx1[1]]

			stack[:,:,n] = imgb
		return np.amin(stack,2)


	def printMeta(self):
		print "Alpha values:\n",  self.alphavals
		print "Beta values:\n",  self.betavals

	def getMetaValues(self):
		return self.alphavals, self.betavals

	def getMetaArray(self):
		return self.meta

	def getBG(self):
		return self.bg_combined

	def makeMultiColor(self, data, data2):
		img = np.zeros((len(data[:, 0]), len(data[0, :]), 4),  dtype=np.uint8)
		img[:, :, 3] = 255
		img[:, :, 1] = 255*data/np.max(data)
		img[:, :, 2] = 255*data2/np.max(data2)
		img[:, :, 0] = 255*data2/np.max(data2)
		return img

	def makeMeanGrid(self):
		tt_vals = sorted(list(set(self.meta[:, 0])))
		theta_vals = sorted(list(set(self.meta[:, 1])))
		tt_step_size = (max(tt_vals)-min(tt_vals))/len(tt_vals)
		theta_step_size = (max(theta_vals)-min(theta_vals))/len(theta_vals)
		tt_vals_max = max(tt_vals)
		theta_vals_max = max(theta_vals)

		if theta_step_size == 0:
			theta_step_size = 1
			theta_vals_max = min(theta_vals)+1
		if tt_step_size == 0:
			tt_step_size = 1
			tt_vals_max = min(tt_vals)+1

		t2t_grid_x,  t2t_grid_y = np.ogrid[min(theta_vals):theta_vals_max:theta_step_size, min(tt_vals):tt_vals_max:tt_step_size]

		hist = np.zeros((len(tt_vals), len(theta_vals)))
		mean = self.getMean()
		for i in range(len(self.meta[:, 0])):
			hist[tt_vals.index(self.meta[i, 0]), theta_vals.index(self.meta[i, 1])] = mean[i]

		return hist,  t2t_grid_x,  t2t_grid_y

	def makeHistogram(self, hist, alpha, beta, savefilename):
		import matplotlib.pylab as plt
		import seaborn as sns
		sns.set_style("white")
		sns.set_context("paper")

		if len(hist[0, 0, :]) == 1:
			fig, ax = plt.subplots(ncols=len(hist[0, 0, :]), figsize=(6 * len(hist[0, 0, :]), 6))
			ax.pcolor(hist[:, :, 0], norm=matplotlib.colors.LogNorm(), cmap='Greens')
			ax.set_xticks(np.arange(hist[:, :, 0].shape[1])+0.5,  minor=False)
			ax.set_xticklabels(beta, rotation=70)
			ax.set_yticks(np.arange(hist[:, :, 0].shape[0])+0.5,  minor=False)
			ax.set_yticklabels(alpha)
			ax.set_xlabel('Theta')
			if self.datatype == 'strain_tt':
				ax.set_ylabel('2Theta')
			if self.datatype == 'strain_eta':
				ax.set_ylabel('Eta')
			if self.datatype == 'topotomo':
				ax.set_ylabel('Phi')
		else:
			fig, ax = plt.subplots(ncols=len(hist[0, 0, :]), figsize=(6 * len(hist[0, 0, :]), 6))
			for i in range(len(hist[0, 0, :])):
				ax[i].pcolor(hist[:, :, i], norm=matplotlib.colors.LogNorm(), cmap='Greens')
				ax[i].set_xticks(np.arange(hist[:, :, i].shape[1])+0.5,  minor=False)
				ax[i].set_xticklabels(beta, rotation=70)
				ax[i].set_yticks(np.arange(hist[:, :, i].shape[0])+0.5,  minor=False)
				ax[i].set_yticklabels(alpha)
				ax[i].set_xlabel('Theta')
				if self.datatype == 'strain_tt':
					ax[i].set_ylabel('2Theta')
				if self.datatype == 'strain_eta':
					ax[i].set_ylabel('Eta')
				if self.datatype == 'topotomo':
					ax[i].set_ylabel('Phi')
		plt.axis('tight')
		plt.tight_layout()
		if self.rank == 0:
			fig.savefig(self.directory + '/%s_%s.pdf' % (savefilename, self.datatype))
		return fig, ax

	def makePlotArray(self, index, bins, xpos, savefilename):
		sns.set_style("white")
		sns.set_context("paper")

		# diff = np.zeros((len(index[:, 0])))

		xoff = np.zeros((len(index[:, 0])))
		yoff = np.zeros((len(index[:, 0])))

		xoff[2] = -15
		xoff[8] = -20
		xoff[9] = -75
		xoff[10] = -80
		xoff[11] = -45
		xoff[12] = -50

		img = self.getImage(index[0, 0], False)
		img = self.rebin(img, bins)
		npix = len(img[:, 0])*len(img[0, :])
		self.pixarr = np.zeros((len(index[:, 0]), npix))
		self.imgarray = np.zeros((len(index[:, 0]), len(img[:, 0]), len(img[0, :])))

		if len(index[0, :]) == 1:
			fig_array, axarr = plt.subplots(nrows=2, ncols=len(index[:, 0]), figsize=(16, 8))
			for i in range(len(index[:, 0])):
				self.adjustROI(xoff[i], yoff[i])
				if self.rank == 0:
					print i
				print self.roi
				img1 = self.getImage(index[i, 0], False)
				img1 = self.rebin(img1, bins)

				self.imgarray[i, :, :] = img1

				self.makeROIAdjustmentArray()

				if self.rank == 0:
					x = int(len(img1[0, :])*xpos)
					ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

					ta[:, :, 3] = 255
					# ta[:, :, 0] = 255*img1/np.max(img1)
					ta[:, :, 1] = 255*img1/np.max(img1)
					# ta[:, :, 1] = 255*img2/np.max(img2)

					axarr[0, i].imshow(ta)
					axarr[0, i].xaxis.set_major_formatter(plt.NullFormatter())
					axarr[0, i].yaxis.set_major_formatter(plt.NullFormatter())
					axarr[0, i].set_title('%.4f %.4f' % (10.992-self.meta[index[i, 0], 1], self.meta[index[i, 0], 0]))
					axarr[0, i].autoscale(False)
					axarr[0, i].plot([x, x], [0, len(img1[:, 0])], color='magenta')

					axarr[1, i].plot(range(len(img1[:, x])), img1[:, x], 'green', alpha=0.5)
					labels = axarr[1, i].get_xticks().tolist()
					axarr[1, i].set_xticklabels(labels, rotation=-90)

		if len(index[0, :]) != 1:
			fig_array, axarr = plt.subplots(nrows=2, ncols=len(index[:, 0]), figsize=(20, 6))
			for i in range(len(index[:, 0])):
				img1 = self.getImage(index[i, 0], False)
				img2 = self.getImage(index[i, 1], False)
				img1 = self.rebin(img1, bins)
				img2 = self.rebin(img2, bins)
				x = int(len(img1[0, :])*xpos)
				ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

				ta[:, :, 3] = 255
				ta[:, :, 0] = 255*img1/np.max(img1)
				ta[:, :, 2] = 255*img1/np.max(img1)
				ta[:, :, 1] = 255*img2/np.max(img2)

				axarr[0, i].imshow(ta)
				axarr[0, i].xaxis.set_major_formatter(plt.NullFormatter())
				axarr[0, i].yaxis.set_major_formatter(plt.NullFormatter())
				axarr[0, i].set_title('%.4f %.4f' % (10.992-self.meta[index[i, 0], 1], self.meta[index[i, 0], 0]))
				axarr[0, i].autoscale(False)
				axarr[0, i].plot([x, x], [0, len(img1[:, 0])], color='white')

				axarr[1, i].plot(img1[:, x], 'magenta', alpha=0.5)
				axarr[1, i].plot(img2[:, x], 'green', alpha=0.5)
		# 		max1, max2 = self.findPeaks(img1[:, x], img2[:, x])
		# 		axarr[1, i].plot([max1, max1], [0, max(img1[:, x])], color='Magenta')
		# 		axarr[1, i].plot([max2, max2], [0, max(img2[:, x])], color='Green')
		# 		# x1 = range(len(img1[:, x]))
		# 		# y1 = popt[0]*np.exp(-(x1-popt[1])**2/(2*popt[2])**2)
		# 		# axarr[1, i].plot(x1, y1, color='Blue')

		# 		axarr[0, i].plot([x-5, x+5], [max1, max1], color='Magenta')
		# 		axarr[0, i].plot([x-5, x+5], [max2, max2], color='Green')
		# 		diff[i] = max2-max1
		# 		# labels = [item.get_text() for item in axarr.get_xticklabels()]
		# 		# axarr[1, i].set_xticklabels(labels, rotation=70)

		# fig,  ax = plt.subplots(figsize=(8, 8))
		# ax.get_xaxis().get_major_formatter().set_useOffset(False)
		# off = diff/2
		# off = off[:-1]
		# offall = np.append(off, -off[::-1])
		# x1 = 10.992 - self.meta[index[:, 0], 1]
		# x1 = x1[:-1]
		# x = np.append(-x1, x1[::-1])
		# print x,  offall
		# # x = 10.99 - self.meta[index[:, 0], 1]
		# ax.set_ylabel('Lattice rotation [degrees]')
		# ax.set_xlabel('Dislocation distance [um]')
		# # ax.set_xticks(offall)
		# ax.scatter(offall*0.08, x, s=50,  marker='o')
		plt.axis('tight')
		# plt.tight_layout()
		# if self.test == False:
		# fig.savefig('output/xdiff_%s_%s_x%g.pdf' % (savefilename, self.datatype, xpos))
		if self.rank == 0:
			fig_array.savefig(self.directory + '/%s_%s_x%g.pdf' % (savefilename, self.datatype, xpos))

	def makeLargePlotArray(self, index, bins, xpos, savefilename):
		sns.set_style("white")
		sns.set_context("paper")

		plt.figure(figsize=(20, 16))
		gs1 = matplotlib.gridspec.GridSpec(8,  8)
		gs1.update(wspace=0.025,  hspace=0.03)

		if len(index[0, :]) == 1:
			for i in range(len(index[:, 0])):
				axarr = plt.subplot(gs1[i])
				print index[i, 0]
				img1 = self.getImage(index[i, 0], False)
				img1 = self.rebin(img1, bins)

				x = int(len(img1[0, :])*xpos)
				ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

				ta[:, :, 3] = 255
				ta[:, :, 1] = img1  # 255*img1/np.max(img1)

				axarr[0, i].imshow(ta)
				axarr[0, i].xaxis.set_major_formatter(plt.NullFormatter())
				axarr[0, i].yaxis.set_major_formatter(plt.NullFormatter())
				axarr[0, i].set_title('%.4f %.4f' % (self.meta[index[i, 0], 1], self.meta[index[i, 0], 0]))
				axarr[0, i].autoscale(False)
				axarr[0, i].plot([x, x], [0, len(img1[:, 0])], color='magenta')

				axarr[1, i].plot(range(len(img1[:, x])), img1[:, x], 'green', alpha=0.5)
				labels = axarr[1, i].get_xticks().tolist()
				axarr[1, i].set_xticklabels(labels, rotation=-90)

		if len(index[0, :]) != 1:
			for i in range(len(index[:, 0])):
				axarr = plt.subplot(gs1[i])
				self.roi[2] = 350+i*15-20  # 1130-i*17-100
				self.roi[3] = 350+i*15+20  # 1130-i*17+100
				img1 = self.getImage(index[i, 0], False)

				img2 = self.getImage(index[i, 1], False)
				img1 = self.rebin(img1, bins)
				img2 = self.rebin(img2, bins)
				x = int(len(img1[0, :])*xpos)
				ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

				ta[:, :, 3] = 255
				ta[:, :, 0] = 255*img1/np.max(img1)
				ta[:, :, 2] = 255*img1/np.max(img1)
				ta[:, :, 1] = 255*img2/np.max(img2)

				axarr.imshow(ta)
				axarr.xaxis.set_major_formatter(plt.NullFormatter())
				axarr.yaxis.set_major_formatter(plt.NullFormatter())
				axarr.set_title('%.4f %.4f' % (self.meta[index[i, 0], 1], self.meta[index[i, 0], 0]))
				axarr.autoscale(False)
				axarr.plot([x, x], [0, len(img1[:, 0])], color='blue')

				axarr = plt.subplot(gs1[i+4*8])
				axarr.plot(img1[:, len(img1[0, :])/2], 'magenta', alpha=0.5)
				axarr.plot(img2[:, len(img2[0, :])/2], 'green', alpha=0.5)
				# labels = [item.get_text() for item in axarr.get_xticklabels()]
				# axarr[1, i].set_xticklabels(labels, rotation=70)

		# plt.axis('tight')
		# plt.tight_layout()
		if self.rank == 0:
			plt.savefig(self.directory + '/%s_%s_x%g.pdf' % (savefilename, self.datatype, xpos))

	def pixelArray(self, img, index, i):
		npix = len(img[:, 0])*len(img[0, :])
		self.pixarr = np.zeros((len(index[:, 0]), npix))

		img2 = np.reshape(img, (npix))

		print np.shape(self.pixarr)
		print np.shape(img2)

	def findPeaks(self, y1, y2):
		# max2min = 20
		max2 = np.argmax(y2)  # [max2min:45])
		max1 = np.argmax(y1[:max2])
		# print max1
		# print y2[max1+10:]
		max2 = np.argmax(y2[max1+2:])+max1+2  # [max2min:45])
		# x = range(len(y1))

		# popt, pcov = self.fitGaussian(x, y2)
		# print max1
		# print popt
		# print pcov

		return max1, max2  # ,popt# max1+max2min, max2+max2min, max2min

	def showArea(self, i1, i2):
		sns.set_style("white")
		sns.set_context("paper")

		index = self.getIndex(float(self.alphavals[i1]), float(self.betavals[i2]))
		img = self.getImage(index[0], True)
		pic = np.zeros((len(img[:, 0]), len(img[0, :]), 4),  dtype=np.uint8)
		pic[:, :, 3] = 255
		pic[:, :, 1] = 255*img/np.max(img)
		# pic[:, :, 2] = 255*img/np.max(img)
		# pic[:, :, 3] = 255
		# pic[:, :, 1] = img# 255*img/np.max(img)

		fig,  ax = plt.subplots(figsize=(8, 8))
		ax.plot([self.roi[0], self.roi[0]], [self.roi[2], self.roi[3]], color='magenta')
		ax.plot([self.roi[1], self.roi[1]], [self.roi[2], self.roi[3]], color='magenta')
		ax.plot([self.roi[0], self.roi[1]], [self.roi[2], self.roi[2]], color='magenta')
		ax.plot([self.roi[0], self.roi[1]], [self.roi[3], self.roi[3]], color='magenta')
		ax.set_title('%g; %g' % (float(self.alphavals[i1]), float(self.betavals[i2])))
		ax.imshow(pic, interpolation='None')
		if self.rank == 0:
			fig.savefig(self.directory + '/fullimg.png')
		# fig.colorbar(img1)
		return fig, ax

	def getProjection(self, data, x0, y0, x1, y1):
		num = 500
		x,  y = np.linspace(x0,  x1,  num),  np.linspace(y0,  y1,  num)
		# zi = scipy.ndimage.map_coordinates(data,  np.vstack((x, y))) # THIS DOESN'T WORK CORRECTLY
		zi = scipy.ndimage.map_coordinates(np.transpose(data),  np.vstack((x, y)))  # THIS SEEMS TO WORK CORRECTLY
		length = math.sqrt((x1-x0)**2+(y1-y0)**2)
		return zi, length

	def makeStrainArrayMPI(self, alldata, bins, length, xr):
		print self.rank, self.size

		def strainRange(data_part, xr):
			strainpic = np.zeros((np.shape(data_part[0, :, :])))

			for i in range(len(data_part[0, :, 0])):
				if self.rank == 0 and 10*float(i)/len(data_part[0, :, 0]) % 1 == 0.0:
					done = 100*(float(i)/len(data_part[0, :, 0]))
					print "Calculation is %g perc. complete..." % done
				for j in range(len(data_part[0, 0, :])):
					try:
						popt, pcov = self.fitGaussian(xr, data_part[:, i, j])
						strain = popt[1]/(10.992)
						if strain >= -0.0002 and strain <= 0.0002:
							strainpic[i, j] = strain
						# else:
						# 	print "Gaussian peak weird."
							# strainpic[i,j] = 0  # np.nan
					except TypeError:
						# strainpic[i,j] = 0
						print i, j, "Gaussian could not be fitted."

			strainpic[0, 0] = self.rank
			return strainpic

		ypix = (self.roi[1]-self.roi[0])/bins

		# Chose part of data set for a specific CPU (rank).
		local_n = ypix/self.size
		istart = self.rank*local_n
		istop = (self.rank+1)*local_n
		local_data = alldata[:, :, istart:istop]

		# Calculate strain on part of data set.
		strainpic_part = strainRange(local_data, xr)

		# CPU 0 (rank 0) combines data parts from other CPUs.
		if self.rank == 0:
			# Make empty arrays to fill in data from other cores.
			recv_buffer = np.zeros((np.shape(alldata[:, :, istart:istop])))
			strainpic = np.zeros((np.shape(alldata[0, :, :])))

			datarank = strainpic_part[0, 0]
			strainpic_part[0, 0] = 0
			strainpic[:, istart:istop] = strainpic_part
			for i in range(1,  self.size):
				self.comm.Recv(recv_buffer,  MPI.ANY_SOURCE)
				datarank = recv_buffer[0][0, 0]
				recv_buffer[0][0, 0] = 0
				strainpic[:, datarank*local_n:(datarank+1)*local_n] = recv_buffer[0]
		else:
			# all other process send their result
			self.comm.Send(strainpic_part)

		# root process prints results
		if self.comm.rank == 0:
			return strainpic

	def makeGaussArrayMPI(self, alldata, bins, length, xr):
		print self.rank, self.size
		if self.rank == 0:
			start = time.time()

		def gaussRange(data_part, xr):
			theta = np.arange(-0.0025*(length/2), 0.0025*(length/2)+.001, .001)
			gaussarray = np.zeros((len(data_part[0, :, 0]), len(data_part[0, 0, :]), 3))

			for i in range(len(data_part[0, :, 0])):

				if self.rank == 0 and 10*float(i)/len(data_part[0, :, 0]) % 1 == 0.0:
					done = 100*(float(i)/len(data_part[0, :, 0]))
					print "Calculation is %g perc. complete..." % done
				for j in range(len(data_part[0, 0, :])):

					popt, pcov = self.fitGaussian(xr, data_part[:, i, j])
					ygauss = self.gaus(theta, popt[0], popt[1], popt[2])
					intpixel = sum([abs(x)*y for x, y in zip(theta, ygauss)])
					gaussarray[i, j, 0] = intpixel
					gaussarray[i, j, 1] = popt[1]
					gaussarray[i, j, 2] = popt[2]

			gaussarray[0, 0] = self.rank
			return gaussarray

		ypix = (self.roi[1]-self.roi[0])/bins
		# xpix = (self.roi[3]-self.roi[2])/bins

		# Chose part of data set for a specific CPU (rank).
		local_n = ypix/self.size
		istart = self.rank*local_n
		istop = (self.rank+1)*local_n
		local_data = alldata[:, :, istart:istop]

		if self.rank == 0:
			end = time.time()
			print "Init time: ", end-start
		# Calculate gaussian on part of data set.
		gaussarray_part = gaussRange(local_data, xr)

		if self.rank == 0:
			# Make empty arrays to fill in data from other cores.
			reclenx = len(alldata[:, 0, 0])
			recleny = len(alldata[0, :, 0])
			reclenz = len(alldata[0, 0, istart:istop])
			reclenr = 3
			recv_buffer = np.zeros((reclenx, recleny, reclenz, reclenr))
			gaussarray = np.zeros((len(alldata[0, :, 0]), len(alldata[0, 0, :]), 3))

			datarank = gaussarray_part[0, 0]
			gaussarray_part[0, 0] = 0
			gaussarray[:, istart:istop, :] = gaussarray_part
			for i in range(1,  self.size):
				self.comm.Recv(recv_buffer,  MPI.ANY_SOURCE)
				datarank = recv_buffer[0][0, 0, 0]
				recv_buffer[0][0, 0, 0] = recv_buffer[0][0, 1, 0]
				recv_buffer[0][0, 0, 1] = recv_buffer[0][0, 1, 1]
				recv_buffer[0][0, 0, 2] = recv_buffer[0][0, 1, 2]
				gaussarray[:, datarank*local_n:(datarank+1)*local_n, :] = recv_buffer[0]

		else:
			# all other process send their result
			self.comm.Send(gaussarray_part)

		# root process prints results
		if self.comm.rank == 0:
			return gaussarray

	def plotPPOS(self, gaussarray, length):
		sns.set_context("talk")

		figstrain, axppos = plt.subplots(3, 1)

		intpic = np.zeros((np.shape(gaussarray[:, :, 0])))
		ppospic = np.zeros((np.shape(gaussarray[:, :, 0])))
		fwhmpic = np.zeros((np.shape(gaussarray[:, :, 0])))

		ppospic[:, :] = gaussarray[:, :, 1]
		fwhmpic[:, :] = 2*math.sqrt(2*math.log(2))*gaussarray[:, :, 2]
		intpic[:, :] = gaussarray[:, :, 0]

		im0 = axppos[0].imshow(intpic, cmap='jet', interpolation='None')
		im1 = axppos[1].imshow(fwhmpic, cmap='jet', interpolation='None')
		im2 = axppos[2].imshow(ppospic, cmap='BrBG', interpolation='None')

		axppos[0].set_title('INT')
		axppos[1].set_title('FWHM')
		axppos[2].set_title('PPOS')

		def fmt(x,  pos):
			a,  b = '{:.2e}'.format(x).split('e')
			b = int(b)
			return r'${} \times 10^{{{}}}$'.format(a,  b)

		figstrain.subplots_adjust(right=0.8)
		cbar_ax0 = figstrain.add_axes([0.85,  0.65,  0.02,  0.2])
		cbar_ax1 = figstrain.add_axes([0.85,  0.4,  0.02,  0.2])
		cbar_ax2 = figstrain.add_axes([0.85,  0.1,  0.02,  0.2])
		figstrain.colorbar(im0, cax=cbar_ax0)  # , format=ticker.FuncFormatter(fmt))
		clb1 = figstrain.colorbar(im1, cax=cbar_ax1)  # , format=ticker.FuncFormatter(fmt))
		clb2 = figstrain.colorbar(im2, cax=cbar_ax2)  # , format=ticker.FuncFormatter(fmt))
		clb1.set_clim(-0.0135, -0.0045)
		clb2.set_clim(-0.0035, 0.0035)

		# linestart = [116, 124]
		# linestop = [193, 30]
		# clb.set_clim(-0.00013, 0.00013)
		# axppos[0].autoscale(False)
		# axppos[0].plot([linestart[0], linestop[0]], [linestart[1], linestop[1]])

		# z, length = self.getProjection(strainpic, linestart[0], linestart[1], linestop[0], linestop[1])
		# f3, ax3 = plt.subplots()
		# print length
		# try:
		# 	linerange = np.arange(0, 90*length/1000, 90*length/499/1000.)
		# 	print len(linerange), len(z)
		# 	ax3.plot(linerange, z)
		# except ValueError:
		# 	linerange = np.arange(0, 90*length/1000, 90*length/500/1000.)
		# 	print len(linerange), len(z)
		# 	ax3.plot(linerange, z)
		# ax3.set_ylabel(r'Strain [$\Delta\theta/\theta$]')
		# ax3.set_xlabel(r'[$\mu m$]')

		if self.rank == 0:
			figstrain.savefig(self.directory + '/int-ppos-map_%s_%s-%s-%s-%s.pdf' % (self.datatype, self.roi[0], self.roi[1], self.roi[2], self.roi[3], ))
		return figstrain, axppos

	def plotStrainHeatmap(self, alldata):
		img = np.zeros((len(alldata[0, :, 0]), len(alldata[0, 0, :]), 4),  dtype=np.uint8)

		# for i in range(len(alldata[:, 0, 0])/2+1):
		# 	rawdata = alldata[i, :, :]
		# 	print i, np.max(rawdata)
		# 	img[:, :, 3] = 255  # -20*i
		# 	img[:, :, 1] = (255-i*10)*rawdata/np.max(rawdata)
		# 	img[:, :, 2] = 0  # 255*rawdata/np.max(rawdata)
		# 	img[:, :, 0] = 0  # 255*rawdata/np.max(rawdata)
		#
		# for i in range(len(alldata[:, 0, 0])/2):
		# 	print i+len(alldata[:, 0, 0])/2+1, np.max(rawdata)
		# 	rawdata = alldata[i+len(alldata[:, 0, 0])/2+1, :, :]
		# 	img[:, :, 3] = 255  # -20*i
		# 	img[:, :, 1] = 0
		# 	img[:, :, 2] = (180+i*10)*rawdata/np.max(rawdata)
		# 	img[:, :, 0] = (180+i*10)*rawdata/np.max(rawdata)

		img = alldata[10, :, :]

		print np.mean(img), np.max(img)
		print np.mean(alldata[2, :, :]), np.max(alldata[2, :, :])

		linestart = [999, 0]
		linestop = [0, 999]

		z, length = self.getProjection(img, linestart[0], linestart[1], linestop[0], linestop[1])

		fig, ax = plt.subplots()
		fig2, ax2 = plt.subplots()

		ax2.plot(z)

		ax.imshow(img)

	def adjustGradient(self):
		xpix = (self.roi[3]-self.roi[2])
		ypix = (self.roi[1]-self.roi[0])
		length = int(math.sqrt(2*ypix**2))
		gradient = np.ones((length, length))
		for i in range(length):
			gradient[i, :] *= -0.0001 + 2*i*0.0001/length

		gradient = scipy.ndimage.interpolation.rotate(gradient, 125)
		# im2 = axstrain[1].imshow(gradient)
		l = len(gradient)
		gradient = gradient[l/2-xpix/2:l/2+xpix/2, l/2-ypix/2:l/2+ypix/2]
		return gradient

	def plotStrain(self, strainpic):
		import matplotlib.ticker as ticker
		sns.set_context("talk")
		figstrain, axstrain = plt.subplots(2, 1)
		# im = axstrain[0].imshow(strainpic, cmap="BrBG")
		# im2 = axstrain[1].imshow(self.imgarray[2, :, :], cmap="Greens")

		gradient = self.adjustGradient()

		# im = axstrain[0].imshow(strainpic, cmap="BrBG")
		im = axstrain[0].imshow(strainpic+gradient, cmap="BrBG")

		im2 = axstrain[1].imshow(self.imgarray[2, :, :], cmap="Greens")

		axstrain[0].set_title("%g %g %g %g" % (self.roi[0], self.roi[1], self.roi[2], self.roi[3]))
		axstrain[0].set_title(r'$\epsilon_{220}$')

		def fmt(x,  pos):
			a,  b = '{:.2e}'.format(x).split('e')
			b = int(b)
			return r'${} \times 10^{{{}}}$'.format(a,  b)

		figstrain.subplots_adjust(right=0.8)
		cbar_ax1 = figstrain.add_axes([0.85,  0.55,  0.02,  0.35])
		cbar_ax2 = figstrain.add_axes([0.85,  0.1,  0.02,  0.35])
		clb = figstrain.colorbar(im, cax=cbar_ax1, format=ticker.FuncFormatter(fmt))
		figstrain.colorbar(im2, cax=cbar_ax2)

		linestart = [275, 110]
		linestop = [275, 175]
		clb.set_clim(-0.00013, 0.00013)
		axstrain[0].autoscale(False)
		axstrain[0].plot([linestart[0], linestop[0]], [linestart[1], linestop[1]])

		z, length = self.getProjection(strainpic, linestart[0], linestart[1], linestop[0], linestop[1])
		f3, ax3 = plt.subplots()
		linerange = np.linspace(0, 90*length/1000, len(z))
		ax3.plot(linerange, z)
		ax3.set_ylabel(r'Strain [$\Delta\theta/\theta$]')
		ax3.set_xlabel(r'[$\mu m$]')

		if self.rank == 0:
			np.save(self.directory + '/strainmap_array.txt', strainpic+gradient)
			f3.savefig(self.directory + '/strainmap_line.pdf')
			figstrain.savefig(self.directory + '/strainmap.pdf')
			f4, ax4 = plt.subplots()
			strain = np.reshape(strainpic, len(strainpic[:, 0])*len(strainpic[0, :]))
			# strain[strain<-0.0005] = 0
			# strain[strain>0.0005] = 0
			sns.distplot(strain, kde=False,  rug=False)
			ax4.set_xlim(np.min(strain)-abs(0.1*np.min(strain)), np.max(strain)+0.1*np.max(strain))
			# ax4.set_xlim(-0.0004,0.0004)
			ax4.set_xlabel(r'$\theta$ offset [$^o$]')
			ax4.set_title('Strain distribution')
			f4.savefig(self.directory + '/straindistribution.pdf')
		return figstrain, axstrain

if __name__ == '__main__':
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	if rank == 0:
		start = time.time()


	path = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/strain'

	filename = 'strainmap_tt_2'
	sampletitle = filename
	bg_filename = 'bg1_5s_'

	datatype = 'strain_tt'

	poi = [1250, 1150]
	size = [600, 300]

	roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

	# Initialize the class as 'dfxm'.
	dfxm = GetEdfData(path, filename, bg_filename, roi, datatype)

	dfxm.setTest(True)
	dfxm.adjustOffset(True)

	# To get alpha and beta values present in the data set:
	a, b = dfxm.getMetaValues()

	# To get an array of all pictures. First column is the 2Theta value, 2nd row is theta value.
	meta = dfxm.getMetaArray()
	print a, b
	# Get the index of
	index = dfxm.getIndex(22.046, 10.9805)
	img = dfxm.getImage(index[0], False)

	# index = test2.getIndex(1., 10.9885)
	img2 = dfxm.getImage(index[0])

	# Subtracting bg
	img2 -= bg
	# Binning 2x2
	img2 = dfxm.rebin(img2, 2)

	# Subtracting bg
	img -= bg
	# Binning 2x2
	img = dfxm.rebin(img, 2)
	# img = test2.makeMultiColor(img, img2)

	plt.imshow(img)
	plt.colorbar()
	plt.show()
