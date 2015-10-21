import sys
import math
import numpy as np
import EdfFile
from os import listdir
from os.path import isfile, join
import matplotlib
import matplotlib.pylab as plt

import scipy.signal as sp
from scipy import ndimage
import seaborn as sns

from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

import cv2


def getFilelists(mypath,filename,bg_filename,phi_val):
	onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]

	data_files = []
	bg_files = []

	for k in onlyfiles[:]:
		if phi_val != -1:
			if k[:len(filename)] == filename and float(k[-13:-9]) == phi_val:
				data_files.append(k) 
		else:
			if k[:len(filename)] == filename:
				data_files.append(k) 
		

	for k in onlyfiles:
		if k[:len(bg_filename)] == bg_filename:
			bg_files.append(k)

	return data_files,bg_files


def getImgSize(roi_xmin,roi_xmax,roi_ymin,roi_ymax):
	imgsize_x = roi_xmax-roi_xmin
	imgsize_y = roi_ymax-roi_ymin
	roi = [roi_xmin,roi_xmax,roi_ymin,roi_ymax]
	return imgsize_x,imgsize_y,roi

def getBGarray(mypath,bg_files,roi):
	bg_file_with_path = mypath + '/' + bg_files[0]
	bg_class = EdfFile.EdfFile(bg_file_with_path)
	bg_img = bg_class.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]#

	bg_combined = np.zeros(np.shape(bg_img))




	print "Reading background files..."

	for i in range(len(bg_files)):
		bg_file_with_path = mypath + '/' + bg_files[i]
		bg_class = EdfFile.EdfFile(bg_file_with_path)
		bg_combined += bg_class.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]

	bg_combined /= len(bg_files)

	return bg_combined


def getAllData(mypath,data_files,roi,phi_val):
	file_with_path = mypath + '/' + data_files[0]
	data_class = EdfFile.EdfFile(file_with_path)
	data_bg = data_class.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]

	data_combined = data_bg-data_bg

	low_threshold = 35
	up_threshold = 200

	data_all = np.zeros((len(data_combined[:,0]),len(data_combined[0,:]),len(data_files)))
	data_mean = np.zeros((len(datafiles)))
	data_or = np.zeros((len(data_files),4))

	a = []
	print "Reading data files..."
	for i in range(len(data_files)):

		if float(data_files[i][-8:-4]) > -1 and float(data_files[i][-13:-9]) > phi_val :
			file_with_path = mypath + '/' + data_files[i]
			img = EdfFile.EdfFile(file_with_path)
			header = img.GetHeader(0)
			data1 = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined
			#data_mean[i] = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]].mean()-bg_combined.mean()

			#low_threshold = 35
			#up_threshold = 200
			#data1[data1 < low_threshold] = 0
			#data1[data1 > up_threshold] = up_threshold

			data_all[:,:,i] = data1

			mot_array = header['motor_mne'].split(' ')
			motpos_array = header['motor_pos'].split(' ')

			
			data_or[i,0] = float(motpos_array[mot_array.index('obpitch')])
			sn = float(data_files[i][-8:-4])
			theta = (11.01-10.97)/40
			print i, sn, data_or[i,0]
			data_or[i,1] = 10.97+theta*sn+theta/2
			data_or[i,2] = sn

	return data_all,data_or,data_mean



def getIndex(data_or):
	print "Getting index..."
	return len(data_all[:,0,0]),len(data_all[0,:,0]),len(list(set(data_or[:,0]))),len(list(set(data_or[:,1]))),len(list(set(data_or[:,2])))

def makePlotArray(data_all,data_or,value):
	tt_vals = sorted(list(set(data_or[:,0])))
	theta_vals = sorted(list(set(data_or[:,1])))
	nrows = int(math.ceil(len(tt_vals)))
	print len(theta_vals), len(tt_vals), math.ceil(len(tt_vals))
	print 'nrows: ' + str(nrows)
	empty_plots = (nrows*4)%len(tt_vals)
	#nrows = int(math.ceil(len(np.where(data_or[:,0] == value)[0]) / 5.))
	#empty_plots = (nrows*4)%len(np.where(data_or[:,0] == value)[0])
	
	plt.figure(figsize = (16,8),dpi=150)
	gs1 = matplotlib.gridspec.GridSpec(nrows, len(theta_vals))
	gs1.update(wspace=0.025, hspace=0.03)

	print gs1

	j = 0

	#for subsm in np.where(data_or[:,0] == value)[0]:
	#for tt in range(len(tt_vals)):
	#	for theta in range(len(theta_vals)):
			#if data_or[subsm,1] < 10.99-0.007 or data_or[subsm,1] > 10.99+0.007:
	for l in range(len(data_or[:,0])):
		#if l%2 == 0:
#		if data_or[l,0] < 22.116 and data_or[l,0] > 22.051:
#			if data_or[l,1] < 10.99:#-0.003 or data_or[l,1] > 10.99+0.01:
		with sns.axes_style('dark'):
			#print l, data_or[l,1], data_or[l,0]
			axs = plt.subplot(gs1[l])
			axs.imshow(data_all[:,:,l],cmap='Blues')
		axs.xaxis.set_major_formatter(plt.NullFormatter())
		axs.yaxis.set_major_formatter(plt.NullFormatter())
		axs.set_aspect('equal')
			#axs.set_title('%.3f' % (data_or[l,0]),fontsize=5)
		j += 1

				#print subsm
	#plt.show()
	plt.savefig('strainmap_eta_100threshold.png',bbox_inches='tight')

def compareTheta(data_all,data_or,eta,cols):
	first = np.where(data_or[:,0] == eta)[0][10]
	second = np.where(data_or[:,0] == eta)[0][36]
	#[print x for np.amax(x) in data_all[:,:,np.where(data_or[:,0] == eta)[0]]]
	#print np.where(data_or[:,0] == eta)
	#first = 7
	#second = 39
	#plt.figure(figsize = (16,8),dpi=150)

	fig,axarr = plt.subplots(nrows=2,ncols=cols)

	for i in range(cols):
		#print data_or[i,0], eta
		#if data_or[i,0] == eta:
		ta = np.ones((len(data_all[:,0,0]),len(data_all[0,:,0]),4), dtype=np.uint8)*0
		ta[:,:,2] = 255*data_all[:,:,first+i]/np.max(data_all[:,:,first+i])
		ta[:,:,3] = 255*data_all[:,:,first+i]/np.max(data_all[:,:,first+i])
		axarr[0,i].imshow(ta,'Greens')
		#axarr[i].imshow(data_all[:,:,first+i],'Greens')
		ta = np.ones((len(data_all[:,0,0]),len(data_all[0,:,0]),4), dtype=np.uint8)*0
		ta[:,:,0] = 255*data_all[:,:,second-i]/np.max(data_all[:,:,second-i])
		ta[:,:,3] = 255*data_all[:,:,second-i]/np.max(data_all[:,:,second-i])
		#ta[np.where(ta[:,:,3] > 100)] = 255
		axarr[0,i].imshow(ta, alpha=0.8)#,interpolation='bilinear')
		#axarr[i].imshow(data_all[:,:,second-i],'Reds',alpha=0.5,interpolation='bilinear')
		axarr[0,i].xaxis.set_major_formatter(plt.NullFormatter())
		axarr[0,i].yaxis.set_major_formatter(plt.NullFormatter())
		axarr[0,i].set_title('%.3f' % (data_or[first+i,1]),fontsize=18)
		#axarr[0,i].plot([100,100],[10,190])
		#axarr[0,i].set_axis_off()

		axarr[1,i].plot(data_all[:,100,first+i],'b',alpha=0.5)

		axarr[1,i].plot(data_all[:,100,second-i],'r',alpha=0.5)

	plt.show()
		#plt.savefig('strainmap_comparetheta/strainmap_eta_%g_comparetheta.png' % eta,bbox_inches='tight')

def compareTwoTheta(data_all,data_or,eta,theta,cols):
	print eta[0]
	fig,axarr = plt.subplots(nrows=2,ncols=cols,figsize=(15,5))

	for i in range(cols):
		print i, np.where(data_or[:,0] == eta[i])[0][7+theta], np.where(data_or[:,0] == eta[i])[0][39-theta]
		
		ta = np.ones((len(data_all[:,0,0]),len(data_all[0,:,0]),4), dtype=np.uint8)*0
		ta[:,:,2] = 255*data_all[:,:,np.where(data_or[:,0] == eta[i])[0][7+theta]]/np.max(data_all[:,:,np.where(data_or[:,0] == eta[i])[0][7+theta]])
		ta[:,:,3] = 255*data_all[:,:,np.where(data_or[:,0] == eta[i])[0][7+theta]]/np.max(data_all[:,:,np.where(data_or[:,0] == eta[i])[0][7+theta]])
		axarr[0,i].imshow(ta)#,'Greens')
		
		ta = np.ones((len(data_all[:,0,0]),len(data_all[0,:,0]),4), dtype=np.uint8)*0
		ta[:,:,0] = 255*data_all[:,:,np.where(data_or[:,0] == eta[i])[0][39-theta]]/np.max(data_all[:,:,np.where(data_or[:,0] == eta[i])[0][39-theta]])
		ta[:,:,3] = 255*data_all[:,:,np.where(data_or[:,0] == eta[i])[0][39-theta]]/np.max(data_all[:,:,np.where(data_or[:,0] == eta[i])[0][39-theta]])
		
		axarr[0,i].imshow(ta, alpha=0.8)#,interpolation='bilinear')
		
		axarr[0,i].xaxis.set_major_formatter(plt.NullFormatter())
		axarr[0,i].yaxis.set_major_formatter(plt.NullFormatter())
		axarr[0,i].set_title('%.3f' % eta[i],fontsize=8)
		
		axarr[1,i].plot(data_all[:,100,np.where(data_or[:,0] == eta[i])[0][7+theta]],'b',alpha=0.5)

		axarr[1,i].plot(data_all[:,100,np.where(data_or[:,0] == eta[i])[0][39-theta]],'r',alpha=0.5)

	#plt.show()
	the = data_or[np.where(data_or[:,0] == eta[i])[0][7+theta],1]
	#axarr.set_aspect('equal')
	#plt.axis('tight')
	#plt.tight_layout()
	plt.savefig('strainmap_comparetheta/strainmap_theta_%g_comparetwotheta.png' % the,bbox_inches='tight')



print "Filename stuff..."
mypath = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/strain'

filename = 'strainmap_eta_'
sampletitle = filename
bg_filename = 'bg1_5s_'

phi_val= -1

value = 1

roi_xmin = 1100
roi_xmax = 1300
roi_ymin = 700
roi_ymax = 900

imgsize_x,imgsize_y,roi = getImgSize(roi_xmin,roi_xmax,roi_ymin,roi_ymax)
datafiles,bg_files = getFilelists(mypath,filename,bg_filename,phi_val)

bg_combined = getBGarray(mypath,bg_files,roi)
print 'Background files read.'
data_all,data_or,data_mean = getAllData(mypath,datafiles[:700],roi,phi_val)
#np.savetxt('datamean_map2.txt',data_mean)
#data_mean = np.loadtxt('datamean_%s.txt' % sampletitle)

print 'Data files read.'

eta_vals = sorted(list(set(data_or[:,0])))
theta_vals = sorted(list(set(data_or[:,1])))

#eta = eta_vals[6] 

theta = 4
cols = 8
eta_start = 6
eta = eta_vals[eta_start:eta_start+cols] 
#compareTheta(data_all,data_or,eta,cols)
compareTwoTheta(data_all,data_or,eta,3,cols)
compareTwoTheta(data_all,data_or,eta,4,cols)
compareTwoTheta(data_all,data_or,eta,5,cols)
compareTwoTheta(data_all,data_or,eta,6,cols)
compareTwoTheta(data_all,data_or,eta,7,cols)
compareTwoTheta(data_all,data_or,eta,8,cols)
compareTwoTheta(data_all,data_or,eta,9,cols)
compareTwoTheta(data_all,data_or,eta,10,cols)
compareTwoTheta(data_all,data_or,eta,11,cols)
compareTwoTheta(data_all,data_or,eta,12,cols)


