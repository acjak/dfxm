import sys
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


def getAllData(mypath,data_files,data_or,roi,alpha,beta):
	file_with_path = mypath + '/' + data_files[0]
	data_class = EdfFile.EdfFile(file_with_path)
	data_bg = data_class.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]

	data_combined = data_bg-data_bg

	threshold = 0
	upper_thres = 300

	data_all = np.zeros((len(data_combined[:,0]),len(data_combined[0,:]),len(data_files)))
	data_mean = np.zeros((len(datafiles)))


	a = []
	print "Reading data files..."
	for i in range(len(data_files)):

		#if float(data_files[i][-8:-4]) > -1 and float(data_files[i][-13:-9]) > phi_val :
		if data_or[i,0] == alpha:
			file_with_path = mypath + '/' + data_files[i]
			img = EdfFile.EdfFile(file_with_path)
			header = img.GetHeader(0)
			#data_all[:,:,i] = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]
			#data_mean[i] = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]].mean()-bg_combined.mean()
			sn = float(data_files[i][-8:-4])
			print i, sn, data_or[i,0]

	return data_all

def getMeanData(mypath,data_files,data_or,roi):
	# file_with_path = mypath + '/' + data_files[0]
	# data_class = EdfFile.EdfFile(file_with_path)
	# data_bg = data_class.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]

	# data_combined = data_bg-data_bg

	# threshold = 0
	# upper_thres = 300

	#data_all = np.zeros((len(data_combined[:,0]),len(data_combined[0,:]),len(data_files)))
	data_mean = np.zeros((len(datafiles)))


	a = []
	print "Reading data files..."
	for i in range(len(data_files)):
		file_with_path = mypath + '/' + data_files[i]
		img = EdfFile.EdfFile(file_with_path)
		header = img.GetHeader(0)
		#data_all[:,:,i] = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]
		data_mean[i] = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]].mean()-bg_combined.mean()
		sn = float(data_files[i][-8:-4])
		print i, sn, data_or[i,0]

	return data_mean

def getMetaData(mypath,data_files):
	# file_with_path = mypath + '/' + data_files[0]
	# data_class = EdfFile.EdfFile(file_with_path)

	data_or = np.zeros((len(data_files),4))

	print "Reading data files..."
	for i in range(len(data_files)):

		file_with_path = mypath + '/' + data_files[i]
		img = EdfFile.EdfFile(file_with_path)
		header = img.GetHeader(0)

		mot_array = header['motor_mne'].split(' ')
		motpos_array = header['motor_pos'].split(' ')

		data_or[i,0] = int(round(float(motpos_array[mot_array.index('diffrx')]),0))
		sn = float(data_files[i][-8:-4])
		theta = (11.01-10.97)/40
		print i, sn, data_or[i,0]
		data_or[i,1] = 10.97+theta*sn+theta/2
		data_or[i,2] = sn

	return data_or


def getIndex(data_or):
	print "Getting index..."
	return len(data_all[:,0,0]),len(data_all[0,:,0]),len(list(set(data_or[:,0]))),len(list(set(data_or[:,1]))),len(list(set(data_or[:,2])))

def reshapeData(data_all):
	print "Reshaping data from [x,y,u*v*z] to [u,v,x*y*z]..."
	reshaped = np.reshape(np.transpose(np.reshape(data_all,[xnum,ynum,unum,vnum,znum]),(2,3,0,1,4)),[unum,vnum,xnum*ynum*znum])
	print "Last index of array is now " + str(len(reshaped[0,0,:])) + " Long."
	return reshaped

def getImage(data_all,data_or,pmo, sath_val):
	img = np.zeros(np.shape(data_all[:,:,2]))
	for i in range(len(data_all[0,0,:])):
		if data_or[i,0] == pmo:#  and data_or[i,2] == sath_val:
			img += data_all[:,:,i]
			# for j in range(len(img[0,:])):
			# 	for k in range(len(img[:,0])):
			# 		if k != 0 and k != len(img[:,0])-1 and j != 0 and j != len(img[0,:])-1:
			# 			val = img[k,j]
			# 			if val >= 3*img[k-1,j] and val >= 3*img[k+1,j] and val >= 3*img[k,j-1] and val >= 3*img[k,j+1]:
			# 				img[k,j] = 0
			#sns.set_style("white")
	img = img/8
	plt.imshow(img,cmap='Greys')
	plt.colorbar()

	sath = data_or[i,2]
	pmo = data_or[i,0]
	#filename = '/Users/andcj/bitsync/Documents/PostDoc/plots_for_henning/data_outliersremoved_small_nobg_thres30_pmo%g_averaged.jpg'%(pmo)
	#print filename
	#plt.savefig(filename)
	#plt.show()
	return img




def makeGrid(data_mean,data_or,phi_val,sampletitle):
	tt_vals = sorted(list(set(data_or[:,0])))
	theta_vals = sorted(list(set(data_or[:,1])))
	tt_step_size = (max(tt_vals)-min(tt_vals))/len(tt_vals)
	theta_step_size = (max(theta_vals)-min(theta_vals))/len(theta_vals)

	theta_vals_round = theta_vals
	tt_vals_round = tt_vals

	if theta_step_size == 0:
		theta_step_size = 1
		theta_vals_max = min(theta_vals)+1
	if tt_step_size == 0:
		tt_step_size = 1
		tt_vals_max = min(tt_vals)+1

	t2t_grid_x, t2t_grid_y = np.ogrid[min(theta_vals):max(theta_vals):theta_step_size,min(tt_vals):max(tt_vals):tt_step_size]
	
	data_histogram = np.zeros((len(tt_vals),len(theta_vals)))
	index_histogram = np.zeros((len(tt_vals),len(theta_vals)))
	print len(data_mean)
	for i in range(len(data_mean[:])):
		data_histogram[tt_vals.index(data_or[i,0]),theta_vals.index(data_or[i,1])] = data_mean[i]
		index_histogram[tt_vals.index(data_or[i,0]),theta_vals.index(data_or[i,1])] = i

	for i in range(len(theta_vals)):
		theta_vals_round[i] = "%.4f" %(theta_vals[i])

	for i in range(len(tt_vals)):
		tt_vals_round[i] = "%.4f" %(tt_vals[i])

	for i in range(len(tt_vals)):
		max_index = np.argmax(data_histogram[i,:])
		max_val = data_histogram[i,max_index]
		#print data_histogram[i,max_index], max_index
		for j in range(len(theta_vals)):
			
			#print data_histogram[i,j]/data_histogram[i,max_index], data_histogram[i,j], data_histogram[i,max_index], max_index
			data_histogram[i,j] /= max_val


	# fig, ax = plt.subplots()
	# heatmap = ax.pcolor(data_histogram,norm=matplotlib.colors.LogNorm(),cmap='Greens')
	# ax.set_xticks(np.arange(data_histogram.shape[1])+0.5, minor=False)
	# ax.set_xticklabels(theta_vals_round,rotation=70)
	# ax.set_yticks(np.arange(data_histogram.shape[0])+0.5, minor=False)
	# ax.set_yticklabels(tt_vals_round)
	# ax.set_xlabel('Theta')
	# ax.set_ylabel('Phi')
	# #abar = fig.colorbar(heatmap)
	#plt.axis('tight')
	#plt.tight_layout()
	#plt.savefig('topotomo_heatmap_%s.png' % sampletitle)
	#plt.show()
	return data_histogram, t2t_grid_x, t2t_grid_y, theta_vals, index_histogram

def makeOffsetGrid(data_histogram,t2t_grid_x,t2t_grid_y,gaussianpeaks,data_or,sampletitle):
	angle_stepsize = t2t_grid_x[1]-t2t_grid_x[0]
	offset = (gaussianpeaks[:,1]-10.99)/angle_stepsize
	offset_integer = np.around(offset, 0)#map(int,offset)

	offset_x0 = offset_integer-offset
	print offset[20]
	print offset_x0[20]
	print offset_integer[20]
	print gaussianpeaks[20,1]
	print 0.0082712/(t2t_grid_x[1]-t2t_grid_x[0])

	gauss = np.zeros((len(t2t_grid_x)))

	for i in range(len(t2t_grid_x)):
		gauss[i] = exp(-(angle_stepsize*(10.99-(t2t_grid_x[0]+i*angle_stepsize))**2/(2*gaussianpeaks[20,2]**2)))
	#gauss = exp(-(angle_stepsize*(gaussianpeaks[:,1]-10.99)**2/(2*gaussianpeaks[:,2]**2)))

	print gauss



	tt_vals = sorted(list(set(data_or[:,0])))
	theta_vals = sorted(list(set(data_or[:,1])))

	theta_vals_round = theta_vals
	tt_vals_round = tt_vals

	data_histogram = np.zeros((len(tt_vals),len(theta_vals)))
	for i in range(len(data_mean[:])):
		data_histogram[tt_vals.index(data_or[i,0]),theta_vals.index(data_or[i-offset_integer[tt_vals.index(data_or[i,0])],1])] = data_mean[i]


	#normalize so every peak = 1
	for i in range(len(tt_vals)):
		max_index = np.argmax(data_histogram[i,:])
		max_val = data_histogram[i,max_index]
		for j in range(len(theta_vals)):
			data_histogram[i,j] /= max_val
			#data_histogram[i,j] *= exp(-(angle_stepsize*(10.99-(t2t_grid_x[0]+j*angle_stepsize))**2/(2*gaussianpeaks[i,2]**2)))



	#fig, ax = plt.subplots()
	#heatmap = ax.pcolor(data_histogram,norm=matplotlib.colors.LogNorm(),cmap='Greens')
	
	for i in range(len(theta_vals)):
		theta_vals_round[i] = "%.4f" %(theta_vals[i])

	for i in range(len(tt_vals)):
		tt_vals_round[i] = "%.4f" %(tt_vals[i])
	# ax.set_xticks(np.arange(data_histogram.shape[1])+0.5, minor=False)
	# ax.set_xticklabels(theta_vals_round,rotation=70)
	# ax.set_yticks(np.arange(data_histogram.shape[0])+0.5, minor=False)
	# ax.set_yticklabels(tt_vals_round)
	# ax.set_xlabel('Theta')
	# ax.set_ylabel('Phi')
	# ax.set_title('Corrected grid,'+sampletitle)
	#abar = fig.colorbar(heatmap)
	#plt.axis('tight')
	#plt.tight_layout()
	#plt.savefig('topotomo_corrected_heatmap_%s.png' % sampletitle)
	#plt.show()
	return data_histogram, offset_integer

def makeLinePlot(data_histogram,t2t_grid_x, t2t_grid_y):
	fig2, bx = plt.subplots()
	lineplotX = bx.plot(t2t_grid_x,data_histogram[10,:])
	x_formatter = ScalarFormatter(useOffset=False)
	bx.xaxis.set_major_formatter(x_formatter)
	#bx.set_xticklabels(,rotation=70)
	bx.set_xlabel('Theta')
	bx.set_ylabel('Mean intensity')
	bx.set_yscale('log')
	plt.savefig('strain_heatmap_lineplot_x.png')

	fig3, cx = plt.subplots()
	lineplotY = cx.plot(t2t_grid_y[0],data_histogram[:,10])
	y_formatter = ScalarFormatter(useOffset=False)
	cx.yaxis.set_major_formatter(y_formatter)
	#bx.set_xticklabels(,rotation=70)
	cx.set_xlabel('2Theta')
	cx.set_ylabel('Mean intensity')
	#cx.set_yscale('lin')
	plt.savefig('strain_heatmap_lineplot_y_lin.png')

def findGaussians(data_histogram,x):
	x = np.reshape(x,len(x))
	gaussianpeaks = np.zeros((len(data_histogram[:,0]),3))

	for i in range(len(data_histogram[:,0])):
		y = data_histogram[i,:]
		popt,pcov = fitGaussian(x,y,i)
		gaussianpeaks[i,:] = popt
		
	return gaussianpeaks

def gaus(x,a,x0,sigma):
		return a*exp(-(x-x0)**2/(2*sigma**2))

def fitGaussian(x,y,i):
	n = len(x)                          #the number of data
	mean = sum(x*y)/n                   #note this correction
	sigma = sum(y*(x-mean)**2)/n        #note this correction
	try:
		popt,pcov = curve_fit(gaus,x,y,p0=[max(y),x[np.argmax(y)],0.005],maxfev=10000)
		return popt,pcov
	except RuntimeError:
		print "Error - curve_fit failed"




def removeHotPixels(data_new):
	width = 15
	hp = [[130,784],\
				 [750,1390],\
				 [779,1728],\
				 [1443,737],\
				 [405,1157],\
				 [1926,1244],\
				 [9,1763],\
				 [935,36],\
				 [81,375],\
				 [2037,1392],\
				 [1311,1454]]
				 

	for i in range(len(hp)):
		if hp[i][0] <= len(data_new[:,0])-width and hp[i][1] <= len(data_new[0,:])-width:
			data_new[hp[i][1]-width:hp[i][1]+width,hp[i][0]-width:hp[i][0]+width] = 0
	return data_new


def plotThetaOffset(gaussianpeaks,sampletitle):
	plt.clf()
	#ax = plt.Figure()
	y = gaussianpeaks[:,1]-10.99
	x = t2t_grid_y[0][:]
	#plt.plot(x, y, 'g-')
	#plt.fill_between(x,y-gaussianpeaks[:,2]/2,y+gaussianpeaks[:,2]/2,alpha=0.3,color='g')
	#plt.savefig('thetaoffset_%s_test.png' % sampletitle)
	#plt.show()

def createImages(data_histogram,t2t_grid_x,t2t_grid_y,data_mean,datafiles,mypath,bg_combined,roi,theta_vals):
	plt.clf()
	a_first = 0
	mean_before = 1

	img_mean = np.zeros((len(t2t_grid_y[0,:])))

	for i in range(len(t2t_grid_y[0,:])):
		a = np.argmax(data_histogram[i])
		
		file_index = int(index_histogram[i,a])
		print a, theta_vals[a],data_mean[file_index]

		low_threshold = 20
		up_threshold = 200
		if a == 29:
			file_with_path = mypath + '/' + datafiles[file_index-10]
			img = EdfFile.EdfFile(file_with_path)
			data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined

			data1 = rebin(data,4)

			data1[data1 < low_threshold] = 0
			data1[data1 > up_threshold] = up_threshold

			#data /= data.mean()
			f = np.fft.fft2(data1)
			fshift = np.fft.fftshift(f)
			magnitude_spectrum = 20*np.log(np.abs(fshift))

			rows, cols = data1.shape
			crow,ccol = rows/2 , cols/2
			bw = 30
			fshift[crow-bw:crow+bw, ccol-bw:ccol+bw] = 50
			f_ishift = np.fft.ifftshift(fshift)
			img_back = np.fft.ifft2(f_ishift)
			img_back = np.abs(img_back)

			plt.clf()
			plt.ion()
			fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2,figsize=(15,15))
			ax1.imshow(data1,cmap='Greens')
			ax2.imshow(magnitude_spectrum,cmap='Greens')
			ax3.imshow(img_back,cmap='Greens')
			#ax2.colorbar(ax1)
			ax1.set_title(str(theta_vals[a]))
			
			plt.show()

			posx = 200
			posy = 200
			maxval = 0
			linetrace = []
			while posx < len(data[:,0]):
				maxval = 0
				for i in [-2,-1]:
					for j in [-2,-1,0,1,2]:
						if i != 0 and j != 0:
							if data1[posx+i,posy+j] > maxval: 
								maxval = data1[posx+i,posy+j]
								px = posx+i
								py = posy+i
								print data1[posx+i,posy+j], px, py
				if px != posx or py != posy:
					posx = px
					poxy = py
					linetrace.append([posx,posy])
					print linetrace[0][0], linetrace[0][1]
				#print linetrace[0]
					
					ax4.plot(linetrace)
					plt.draw()
				if px == -6 and py == 198:
					print linetrace
					break
				#print posx, posy





			
		#plt.savefig('topomap_test2/leftofpeak7_'+str("%04d" % (i))+'.png')

	# max_val = max(img_mean)
	# for i in range(len(t2t_grid_y[0,:])):
	# 	img_mean[i] /= img_mean[i]

	# plt.plot(t2t_grid_y[0,:],img_mean)
	# plt.show()

def binImg(data,bs):
	binnedArray = np.zeros((len(data[:,0])/bs,len(data[0,:])/bs))
	for i in range(len(binnedArray[:,0])):
		for j in range(len(binnedArray[0,:])):
			binnedArray[i,j] = np.sum(data[i*bs:i*bs+bs,j*bs:j*bs+bs])
	return binnedArray

def rebin(a, bs):
	shape = (a.shape[0]/bs,a.shape[1]/bs)
	sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
	return a.reshape(sh).sum(-1).sum(1)

def getCleanImage(data_histogram,t2t_grid_x,t2t_grid_y,mypath,datafiles,bg_combined):
	# def press(event):
	#     print('press', event.key)
	#     sys.stdout.flush()
	#     if event.key=='x':
	#         #visible = ax1.get_visible()
	#         #ax1.set_visible(not visible)
	#         fig2.canvas.draw()

	plt.clf()
	plt.ion()
	fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2,figsize=(15,15))
	#fig2.canvas.mpl_connect('key_press_event', press)
	plt.show()
	for i in range(len(t2t_grid_y[0,:])):
		a = np.argmax(data_histogram[i])
		file_index = int(index_histogram[i,a])

		file_with_path = mypath + '/' + datafiles[file_index-10]
		img = EdfFile.EdfFile(file_with_path)
		data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined

		low_threshold = 0
		up_threshold = 40

		file_with_path = mypath + '/' + datafiles[file_index-10]
		img = EdfFile.EdfFile(file_with_path)
		data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined

		data[data < low_threshold] = 0
		data[data > up_threshold] = up_threshold

		data /= np.amax(data)

		imag = np.array(255*data,dtype=np.uint8)
		

		#ax = fig.add_subplot(2,2,1)
		ax1.imshow(imag,cmap='Greens')
		ax1.set_title('Omega ='+str(t2t_grid_y[0,i]))

		#imag = cv2.bilateralFilter(imag,1,100,255)
		#bx = fig.add_subplot(2,2,2)
		#ax2.imshow(imag,cmap='Greens')
		#edges = cv2.Canny(imag,1,80,apertureSize = 3)
		#imag[imag> up_threshold] = up_threshold

		#cx = fig.add_subplot(2,2,3)
		#ax3.imshow(imag,cmap='Greens')

		#imag = cv2.Canny(imag,1,80,apertureSize = 3)
		
		#dx = fig.add_subplot(2,2,4)
		#ax4.imshow(imag,cmap='Greens')
		#plt.colorbar(ax)
		#plt.savefig('topomap_test3/leftofpeak2_'+str("%04d" % (i))+'.png')
		#plt.tight_layout()

		
    	#if chr=='q':
        #	break 
		#sys.stdin.read(1)
		raw_input("Press Enter to continue...")
		plt.draw()


def makePlotArray(data_or,data_histogram,index_histogram,mypath,datafiles,sampletitle):
	
	plt.figure(figsize = (20,15))
	# gs1 = matplotlib.gridspec.GridSpec(4, 20)
	# gs1.update(wspace=0.025, hspace=0.03)

	low_threshold = -20	
	up_threshold = 2000
	offset = 8

	cols = 20

	for i in range(len(set(data_or[:2600,0]))):
		with sns.axes_style('dark'):
			# axs = plt.subplot(gs1[i])
			a = np.argmax(data_histogram[i])
			file_index = int(index_histogram[i,a])
			print a, file_index#, data_or[file_index,0], data_or[file_index,1],len(set(data_or[:,0]))

			file_with_path = mypath + '/' + datafiles[file_index-offset]
			img = EdfFile.EdfFile(file_with_path)

			data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined
			data = rebin(data,2)

			file_with_path = mypath + '/' + datafiles[file_index+offset]
			img = EdfFile.EdfFile(file_with_path)

			data2 = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined
			data2 = rebin(data2,2)
	
			#data[data < low_threshold] = 0
			#data[data > up_threshold] = up_threshold

			img = np.zeros((len(data[:,0]),len(data[0,:]),4), dtype=np.uint8)



			img[:,:,3] = 255
			img[:,:,1] = 255*data/np.max(data)
			img[:,:,2] = 255*data2/np.max(data)
			img[:,:,0] = 255*data2/np.max(data)
			# print np.max(data2), np.max(data)
			# print np.mean(data2), np.mean(data)

			# axs.set_title('%.3f' % data_or[file_index,0],fontsize=8)
			# axs.imshow(img)
			plt.figure(figsize = (20,15))
			plt.title('%.3f' % data_or[file_index,0],fontsize=8)
			plt.imshow(img)
			#plt.savefig('topotomo_disl_location/%s_%g.png' % (sampletitle,i),bbox_inches='tight')
			plt.clf()
			plt.show()

		# axs.xaxis.set_major_formatter(plt.NullFormatter())
		# axs.yaxis.set_major_formatter(plt.NullFormatter())
		# axs.set_aspect('equal')


				#print subsm
	#plt.show()
	#plt.savefig('topotomo_single_disl_movement/single_disl_1_roi.png',bbox_inches='tight')

def makePlotArrayZoom(data_or,data_histogram,index_histogram,mypath,datafiles):
	
	plt.figure(figsize = (20,12))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025, hspace=0.03)

	#maxfig = plt.figure()

	max_array = np.zeros((32,2))

	low_threshold = -20	
	up_threshold = 2000
	offset = 4
	xpos = 75
	img_height = 200
	cols = 20
	j = 0
	for i in range(len(set(data_or[:,0]))):
		with sns.axes_style('dark'):
			#axs = plt.subplot(gs1[i])
			#a = np.argmax(data_histogram[i])
			#file_index = int(index_histogram[i,a])
			#print a, file_index#, data_or[file_index,0], data_or[file_index,1],len(set(data_or[:,0]))
			a = np.argmax(data_histogram[i])
			file_index = int(index_histogram[i,a])
			if data_or[file_index,0] > -47. and data_or[file_index,0] < 19. :

				axs = plt.subplot(gs1[j])

				roi[2] = 415+j*17
				roi[3] = 415+j*17+100
				

				file_with_path = mypath + '/' + datafiles[file_index-offset]
				img = EdfFile.EdfFile(file_with_path)

				data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined[:100]
				data = rebin(data,2)

				file_with_path = mypath + '/' + datafiles[file_index+offset]
				img = EdfFile.EdfFile(file_with_path)

				data2 = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined[:100]
				data2 = rebin(data2,2)

				
				#data[data < low_threshold] = 0
				#data[data > up_threshold] = up_threshold

				img = np.zeros((len(data[:,0]),len(data[0,:]),4), dtype=np.uint8)



				img[:,:,3] = 255
				img[:,:,1] = 255*data/np.max(data)
				img[:,:,2] = 255*data2/np.max(data2)
				img[:,:,0] = 255*data2/np.max(data2)

				axs.set_title('%.3f' % data_or[file_index,0],fontsize=8)
				axs.imshow(img)
				axs.autoscale(False)
				axs.plot([xpos,xpos],[0,50])
				axs.xaxis.set_major_formatter(plt.NullFormatter())
				axs.yaxis.set_major_formatter(plt.NullFormatter())
				axs.set_aspect('equal')

				axs = plt.subplot(gs1[j+4*8])
				axs.plot(data[:,xpos],color='Green')
				axs.plot(data2[:,xpos],color='Magenta')
				max_array[j,0] = data_or[file_index,0]
				max_array[j,1],max1,max2,max2min = findPeaks(data_or[file_index,0],data[:,xpos],data2[:,xpos])
				axs.plot([max2,max2],[0,1000],color='Magenta')
				axs.plot([max1,max1],[0,1000],color='Green')
				# localmax = (np.diff(np.sign(np.diff(data2[:,xpos]))) < 0).nonzero()[0] + 1 
				# localmax2 = sp.find_peaks_cwt(data2[:,xpos],np.arange(1,20))
				# #axs.plot(localmax,data2[localmax,xpos],'o')
				# axs.plot(localmax2,data2[localmax2,xpos],'o')
				axs.xaxis.set_major_formatter(plt.NullFormatter())
				axs.yaxis.set_major_formatter(plt.NullFormatter())

				axs = plt.subplot(gs1[j])
				axs.plot([xpos-5,xpos+5],[max1,max1],color='Green')
				axs.plot([xpos-5,xpos+5],[max2,max2],color='Magenta')
				j+=1
				


				#print subsm
	#plt.plot(max_array[:,0],max_array[:,1])
	plt.show()
	#plt.savefig('topotomo_single_disl_movement/single_disl_3_tomo.png',bbox_inches='tight')
	plt.clf()
	plt.plot(max_array[:,0],max_array[:,1])
	plt.xlabel('Angle [degrees]')
	plt.ylabel('Strain peak distance')
	#plt.savefig('topotomo_single_disl_movement/single_disl_3_peak_dist.png',bbox_inches='tight')
	plt.show()

def makePlotArrayZoomMap2(data_or,data_histogram,index_histogram,mypath,datafiles):
	
	plt.figure(figsize = (20,12))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025, hspace=0.03)

	#maxfig = plt.figure()

	max_array = np.zeros((32,2))

	low_threshold = -20	
	up_threshold = 2000
	offset = 4
	xpos = 75
	img_height = 200
	cols = 20
	j = 0
	for i in range(len(set(data_or[:,0]))):
		with sns.axes_style('dark'):
			#axs = plt.subplot(gs1[i])
			#a = np.argmax(data_histogram[i])
			#file_index = int(index_histogram[i,a])
			#print a, file_index#, data_or[file_index,0], data_or[file_index,1],len(set(data_or[:,0]))
			a = np.argmax(data_histogram[i])
			file_index = int(index_histogram[i,a])
			if data_or[file_index,0] > 121. and data_or[file_index,0] < 220. :

				axs = plt.subplot(gs1[j])

				roi[2] = 1130-j*17-100
				
				roi[3] = 1130-j*17+100
				

				file_with_path = mypath + '/' + datafiles[file_index-offset]
				img = EdfFile.EdfFile(file_with_path)

				data = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined[:200]
				data = rebin(data,2)

				file_with_path = mypath + '/' + datafiles[file_index+offset]
				img = EdfFile.EdfFile(file_with_path)

				data2 = img.GetData(0).astype(np.int64)[roi[2]:roi[3],roi[0]:roi[1]]-bg_combined[:200]
				data2 = rebin(data2,2)

				
				#data[data < low_threshold] = 0
				#data[data > up_threshold] = up_threshold

				img = np.zeros((len(data[:,0]),len(data[0,:]),4), dtype=np.uint8)



				img[:,:,3] = 255
				img[:,:,1] = 255*data/np.max(data)
				img[:,:,2] = 255*data2/np.max(data2)
				img[:,:,0] = 255*data2/np.max(data2)

				axs.set_title('%.3f' % data_or[file_index,0],fontsize=8)
				axs.imshow(img)
				axs.autoscale(False)
				axs.plot([xpos,xpos],[0,50])
				axs.xaxis.set_major_formatter(plt.NullFormatter())
				axs.yaxis.set_major_formatter(plt.NullFormatter())
				axs.set_aspect('equal')

				# axs = plt.subplot(gs1[j+4*8])
				# axs.plot(data[:,xpos],color='Green')
				# axs.plot(data2[:,xpos],color='Magenta')
				# max_array[j,0] = data_or[file_index,0]
				# max_array[j,1],max1,max2,max2min = findPeaks(data_or[file_index,0],data[:,xpos],data2[:,xpos])
				# axs.plot([max2,max2],[0,1000],color='Magenta')
				# axs.plot([max1,max1],[0,1000],color='Green')
				# # localmax = (np.diff(np.sign(np.diff(data2[:,xpos]))) < 0).nonzero()[0] + 1 
				# # localmax2 = sp.find_peaks_cwt(data2[:,xpos],np.arange(1,20))
				# # #axs.plot(localmax,data2[localmax,xpos],'o')
				# # axs.plot(localmax2,data2[localmax2,xpos],'o')
				# axs.xaxis.set_major_formatter(plt.NullFormatter())
				# axs.yaxis.set_major_formatter(plt.NullFormatter())

				# axs = plt.subplot(gs1[j])
				# axs.plot([xpos-5,xpos+5],[max1,max1],color='Green')
				# axs.plot([xpos-5,xpos+5],[max2,max2],color='Magenta')
				j+=1
				


				#print subsm
	#plt.plot(max_array[:,0],max_array[:,1])
	plt.show()
	#plt.savefig('topotomo_single_disl_movement/single_disl_map2_1_tomo.png',bbox_inches='tight')
	plt.clf()
	plt.plot(max_array[:,0],max_array[:,1])
	plt.xlabel('Angle [degrees]')
	plt.ylabel('Strain peak distance')
	#plt.savefig('topotomo_single_disl_movement/single_disl_map2_1_peak_dist.png',bbox_inches='tight')
	#plt.show()

def findPeaks(x,y1,y2):
	max2min = 20
	max2 = np.argmax(y2[max2min:45])
	max1 = np.argmax(y1[max2min:max2+max2min])

	# print sp.find_peaks_cwt(y2,np.arange(5,40))
	# print (np.diff(np.sign(np.diff(y2))) < 0).nonzero()[0] + 1 
	# test =  np.r_[True, y2[1:] > y2[:-1]] & np.r_[y2[:-1] > y2[1:], True]

	# print y2[test]

	if max2+max2min > max1+max2min:
		dist_max = max2-max1
		print dist_max, max1, max2
	else: 
		dist_max = 'NaN'
	return dist_max,max1+max2min,max2+max2min,max2min

print "Filename stuff..."
mypath = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/mapping'

filename = 'map2'
sampletitle = filename
bg_filename = 'bg1_5s_'

phi_val= -1

roi_xmin = 600
roi_xmax = 1400
roi_ymin = 0
roi_ymax = 2048

# roi_xmin = 700
# roi_xmax = 1100
# roi_ymin = 400
# roi_ymax = 1800

# roi_xmin = 1100
# roi_xmax = 1300
# roi_ymin = 700
# roi_ymax = 900

# roi_xmin = 0
# roi_xmax = 1700
# roi_ymin = 0
# roi_ymax = 1700

imgsize_x,imgsize_y,roi = getImgSize(roi_xmin,roi_xmax,roi_ymin,roi_ymax)
datafiles,bg_files = getFilelists(mypath,filename,bg_filename,phi_val)

bg_combined = getBGarray(mypath,bg_files,roi)
print 'Background files read.'
data_or = getMetaData(mypath,datafiles)
print 'Meta data read.'
#data_mean = getMeanData(mypath,datafiles,data_or,roi)
#np.savetxt('datamean_map2.txt',data_mean)
data_mean = np.loadtxt('datamean_%s.txt' % sampletitle)

#data_all = getAllData(mypath,datafiles,data_or,roi,theta,omega)

print 'Data files read.'
#data_all = removeHotPixels(data_all)

#img_sum = getSingleImage(data_or,theta,omega)

data_histogram,t2t_grid_x, t2t_grid_y, theta_vals, index_histogram = makeGrid(data_mean,data_or,phi_val,sampletitle)



#gaussianpeaks = findGaussians(data_histogram,t2t_grid_x)

#plotThetaOffset(gaussianpeaks,sampletitle)

#getCleanImage(data_histogram,t2t_grid_x,t2t_grid_y,mypath,datafiles,bg_combined)


#data_hist_offset = makeOffsetGrid(data_histogram,t2t_grid_x,t2t_grid_y,gaussianpeaks,data_or,sampletitle)

makePlotArray(data_or,data_histogram,index_histogram,mypath,datafiles,sampletitle)
#makePlotArrayZoom(data_or,data_histogram,index_histogram,mypath,datafiles)
#makePlotArrayZoomMap2(data_or,data_histogram,index_histogram,mypath,datafiles)

#createImages(data_histogram,t2t_grid_x,t2t_grid_y,data_mean,datafiles,mypath,bg_combined,roi,theta_vals)


