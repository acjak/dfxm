#!/usr/local/bin/python

## To use it on the gauss2 computer:
# Change first line to: #!/users/blissadm/bin/python

import numpy as np
from scipy import ndimage, signal
import EdfFile
import sys
from os import listdir
from os.path import isfile,  join
import matplotlib.pyplot as plt

def getFilelists(filename):
    onlyfiles = [f for f in listdir('./') if isfile(join('./',  f))]

    data_files = []

    for k in onlyfiles:
        if k[:len(filename)] == filename:
            data_files.append(k)

    return data_files

def findBrightSpot(img):
    r = 20
    spotmax = 0
    maxpos = [1000,1000]

    weight = np.ones((5,5))*30000.

    # correlation = ndimage.correlate(weight,img)
    # matches = (correlation-np.mean(correlation)) > 1*np.std(correlation)

    # hotspot = signal.convolve2d(img,weight)
    #
    # print np.argmax(hotspot,0)
    # print np.argmax(hotspot,1)
    #
    # maxpos[0] = [np.max(np.argmax(hotspot,0))]
    # maxpos[1] = [np.max(np.argmax(hotspot,1))]

    for i in range(len(img[:,0])):
        for j in range(len(img[0,:])):
            if i > r and j > r and i < len(img[:,0])-r and j < len(img[0,:])-r:
                spot_brightness = np.mean(img[i-r:i+r,j-r:j+r])
                if spot_brightness >= spotmax:
                    spotmax = spot_brightness
                    maxpos = [i,j]

    print spotmax, maxpos
    # print img[maxpos[0]-r:maxpos[0]+r,maxpos[1]-r:maxpos[1]+r]
    return spotmax, maxpos

data_files = getFilelists(sys.argv[1])

f = EdfFile.EdfFile(data_files[0])
data = f.GetData(0).astype(np.int64)
imgstack = np.zeros((len(data[:,0]),len(data[0,:]), len(data_files)))
n = 0
for fname in data_files:
    f = EdfFile.EdfFile(fname)
    imgstack[:,:,n] = f.GetData(0).astype(np.int64)
    print fname
    n += 1

# meanimg = np.mean(imgstack, 2)
imgmin = np.min(imgstack,2)



plt.imshow(imgmin)
# plt.autoscale(False)
if len(sys.argv) > 2:
    if sys.argv[2] == '1':
        spotmax, maxpos = findBrightSpot(imgmin)
        plt.scatter(maxpos[1], maxpos[0], s=100, color='magenta')
plt.show()
