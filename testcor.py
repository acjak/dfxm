#!/usr/local/bin/python

import numpy as np
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

meanimg = np.mean(imgstack, 2)

plt.imshow(np.min(imgstack,2))
plt.show()
