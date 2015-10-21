#!/bin/python

from lib.getedfdata import *
import numpy as np

import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

def followDislocation(alpha,beta,bins,xpos):
	print hist[10,:]
	startpos = [0,12]
	length = 32
	offset = 6
	index = np.zeros((length,2),dtype=np.int16)
	for i in range(length):
		#i1 = startpos[0]+2*i
		i1 = i
		i2 = np.argmax(hist[i,:])-offset
		i2_2 = np.argmax(hist[i,:])+offset
		#i2 = startpos[1]+1*i

		hist[i1,i2] = 3000
		data.roi[2] = 350+i*17-100#1130-i*17-100
		data.roi[3] = 350+i*17+100#1130-i*17+100
		index[i,0] = data.getIndex(float(alpha[i1]),float(beta[i2]))[0]
		index[i,1] = data.getIndex(float(alpha[i1]),float(beta[i2_2]))[0]

	data.makeLargePlotArray(index,bins,xpos,'diagplot_%g_%g' % (startpos[0],startpos[1]))
	hist_array[:,:,0] = hist
	fig,ax = data.makeHistogram(hist_array,alpha,beta,'diagplot_hist_%g_%g' % (startpos[0],startpos[1]))
	data.showArea(i1,i2)

path = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/mapping'

filename = 'map1'
sampletitle = filename
bg_filename = 'bg1_5s_'

datatype = 'topotomo'

roi = [1100,1300,250,450]

eta = -0.01
cols = 8

data = GetEdfData(path,filename,bg_filename,roi,datatype)
data.printMeta()
meta = data.getMetaArray()
print "Getting histogram data..."
data.setTest(True)
# data.adjustOffset(True)
hist,datx,daty = data.makeMeanGrid()
a,b = data.getMetaValues()
hist_array = np.zeros((len(hist[:,0]),len(hist[0,:]),1))
hist_array[:,:,0] = hist
print "Making histogram plot..."

data.adjustOffset(False)



followDislocation(a,b,1,0.5)
#diagonalPlot(meta,a,b,hist,1,0.5,datatype)
plt.show()