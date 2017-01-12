# !/bin/python
"""blah."""

import EdfFile
import matplotlib.pyplot as plt
import numpy as np

directory = '/Users/andcj/hxrm_data/nanorod/'
focusdir = '/Users/andcj/hxrm_data/mll_focus/'

focusclass = EdfFile.EdfFile(focusdir + 'ff_mll_3_0051.edf')
imgclass = EdfFile.EdfFile(directory + 'static_2_0001.edf')
ffclass = EdfFile.EdfFile(directory + 'static_2_flatfield_0001.edf')

img_focus = focusclass.GetData(0).astype(np.int64)
img = imgclass.GetData(0).astype(np.int64)
img_fft = np.fft.fft2(img)

print np.shape(np.abs(img_fft))

ff = ffclass.GetData(0).astype(np.int64)

cr = img * ff.mean()/ ff

print cr.max()

cr[cr > 600] = 600
# cr[cr < 0] = 0

# cr_fft = np.fft.fft2(cr)
# fftmap = cr_fft
# print fftmap.mean()
# fftmap[fftmap > 1000000] = 1000000
# print fftmap.mean()
# fftmap = np.fft.fftshift(fftmap)
#
# img[img > 2000] = 2000
# img[img < 0] = 0
#
# fftmap[600:998, 1004:1035] = 0.
# fftmap[1066:1400, 1004:1035] = 0.
#
# fftmap = np.fft.ifftshift(fftmap)
# fftmap = np.fft.ifft2(fftmap)
#
# ff[ff > 2000] = 2000
# ff[ff < 0] = 0
#
# img_focus[img_focus > 200] = 200
# img_focus[img_focus < 0] = 0

plt.imshow(cr, cmap='Greys')

plt.show()
