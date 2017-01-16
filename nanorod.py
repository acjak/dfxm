# !/bin/python
"""blah."""

import EdfFile
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
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

cr = img * ff.mean() / ff

print cr.max()

cr[cr > 600] = 600
img[img > 1000] = 1000
img[img < 0] = 0
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

fig, ax = plt.subplots(1, 2, figsize=(14, 7))

s = 2048

ax[0].imshow(img, cmap='Greys')
ax[0].autoscale(False)
ax[0].plot(
	[s - 100 - 300, s - 100],
	[s - 100, s - 100],
	linewidth=5,
	color='blue')
ax[0].text(s - 420, s - 150, u'5 um', color='blue', fontsize=16)

ax[1].imshow(cr, cmap='Greys')
ax[1].autoscale(False)
ax[1].plot(
	[s - 100 - 300, s - 100],
	[s - 100, s - 100],
	linewidth=5,
	color='blue')
ax[1].text(s - 420, s - 150, u'5 um', color='blue', fontsize=16)

fig.savefig('output/nanorod.png')

plt.show()
