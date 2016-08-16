import numpy as np
import EdfFile
import tomopy
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt



def getFilelists(mypath,filename):
	onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]

	data_files = []

	for k in onlyfiles[:]:
		if k[:len(filename)] == filename:
			data_files.append(k)

	return data_files

def makeDataArray(data_files, mypath, roi):
	im0 = EdfFile.EdfFile(mypath + '/' + data_files[0])
	d = im0.GetStaticHeader(0)
	print roi
	# arr = np.zeros((len(data_files), int(d['Dim_2']), int(d['Dim_1'])))
	arr = np.zeros((len(data_files), int(roi[3]-roi[2]), int(roi[1]-roi[0])))
	for (i, ar) in enumerate(arr):
		#print "Reading file number:", i
		im = EdfFile.EdfFile(mypath + '/' + data_files[i])
		arr[i, :, :] = im.GetData(0).astype(np.int64)[roi[2]:roi[3], roi[0]:roi[1]]
	return arr

poi = [512, 512]

size = [1000, 1000]


roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

mypath = '/Users/andcj/hxrm_data/diamond_topo/'
filename = 'ff1'

data_files = getFilelists(mypath, filename)

arr = makeDataArray(data_files, mypath, roi)

plt.imshow(arr[:, 0, :], cmap="Greys_r")
plt.show()
