#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import time
from lib.getedfdata import *
from lib.dfxm import *
import numpy as np
import itertools
from mpi4py import MPI
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
from matplotlib.colors import colorConverter
from scipy import ndimage
from mpl_toolkits.axes_grid.inset_locator import inset_axes


def fullimage(directory, roi, a, b, c, title=None):
	ind = data.getIndex(a, b, c)[0]
	img = data.getImage(ind, True)
	img_small = data.getImage(ind, False)

	fig0, ax0 = plt.subplots(1, 2, figsize=(8, 4))
	ax0[0].imshow(img, cmap='Greens')
	ax0[0].autoscale(False)
	ax0[0].plot([roi[0], roi[1]], [roi[2], roi[2]], 'b')
	ax0[0].plot([roi[0], roi[0]], [roi[2], roi[3]], 'b')
	ax0[0].plot([roi[1], roi[1]], [roi[2], roi[3]], 'b')
	ax0[0].plot([roi[0], roi[1]], [roi[3], roi[3]], 'b')

	ax0[1].imshow(img_small, cmap='Greens')

	if title is None:
		fig0.savefig(directory + '/fullimage.png')
	else:
		fig0.savefig('{}/fullimage_{}.png'.format(directory, title))

	return img_small


def actualangles(alp, b, com, theta):
	arang = max(alp) - min(alp)
	brang = max(b) - min(b)
	tt = 2 * theta - 0.004
	# alpha = np.arange(-arang / 2, arang / 2, arang / len(alp))
	alpha = np.linspace(-arang / 2, arang / 2, len(alp)) + tt
	beta = np.linspace(-brang / 2, brang / 2, len(b)) + theta
	# beta = np.arange(-brang / 2, brang / 2, brang / len(b))
	return alpha, beta


def histogram(alpha, beta, omega, theta):
	tools.getMeanData(data.bg_combined)
	hist, t2t_grid_x, t2t_grid_y = tools.makeMeanGrid(tools.data_mean, alpha, beta, omega)

	com = list(ndimage.measurements.center_of_mass(hist))
	# compos = [np.rint(com[0]), np.rint(com[1])]
	atrue, btrue = actualangles(alpha, beta, com, theta)

	# a0 = com[0] * (max(a) - min(a)) / len(a) + min(a)
	# b0 = com[1] * (max(b) - min(b)) / len(b) + min(b)

	fig, ax = tools.makeHistogram(hist, atrue, btrue, 'meangrid', com)

	dynamicfields = hist
	dynamicfields[dynamicfields < hist.mean()] = 0.
	dynamicfields[dynamicfields >= hist.mean()] = 1.

	color2 = colorConverter.to_rgba('red', alpha=0.5)
	cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list(
		'my_cmap2',
		[color2, color2],
		256)

	ax.pcolor(
		dynamicfields,
		norm=matplotlib.colors.LogNorm(),
		cmap=cmap2)

	return fig, ax, atrue, btrue


def add_subplot(fig, ax, ai, bi, ci, apind, ind):
	index = data.getIndex(ai, bi, ci)
	img0 = data.getImage(index[0], False)

	xpos = 0.18 + 0.07 * ind
	subax = fig.add_axes([xpos, 0.35, 0.06, 0.08], axisbg='white')
	subax.tick_params(axis=u'both', which=u'both', length=0, labelsize=8)

	line = img0[:, np.shape(img0)[0] / 2]
	subax.plot(line, color='black')

	# subax.imshow(img0)
	subax.set_xticklabels([])
	subax.set_yticklabels([])

	p0 = subax.get_position()
	invfig = ax.transData.inverted()
	p1 = invfig.transform(fig.transFigure.transform(p0))

	subax.plot(
		[np.shape(img0)[0] / 2 - 2, np.shape(img0)[0] / 2 - 2],
		[0, max(img0[:, np.shape(img0)[0] / 2]) * 1.1],
		'r-')

	ax.plot(
		[p1[0, 0] + (p1[1, 0] - p1[0, 0]) / 2, apind[1] + 0.5],
		[p1[0, 1] - (p1[0, 1] - p1[1, 1]), apind[0] + 0.5],
		'w-', transform=ax.transData, linewidth=2)

	ax.plot(
		[p1[0, 0] + (p1[1, 0] - p1[0, 0]) / 2, apind[1] + 0.5],
		[p1[0, 1] - (p1[0, 1] - p1[1, 1]), apind[0] + 0.5],
		'k-', transform=ax.transData, linewidth=1.5)

	subax2 = fig.add_axes([xpos, 0.25, 0.06, 0.08], axisbg='white')
	subax2.set_xticklabels([])
	subax2.set_yticklabels([])
	subax2.imshow(img0)
	subax2.autoscale(False)
	midp = np.shape(img0)[0] / 2
	subax2.plot([midp, midp], [0, midp * 2], 'w-')

	return line, fig, ax


def add_arrow(fig, ax, ai, bi, ci, ind):
	ax.arrow(xpos + 0.08 / 2, 0.25 + 0.08, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
	return fig, ax


path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/ff_strain'
bg_path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff2_'
filename2 = 'nf_'
sampletitle = filename
bg_filename = 'bg_'

datatype = 'strain_tt'

theta = 10.25

test_switch = True


poi = [600, 500]
size = [50, 50]
s = 40

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]

print roi

data = GetEdfData(
	path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0

tools = DFXM(
	path,
	data.data_files,
	directory,
	roi,
	datatype,
	data.dirhash,
	data.meta,
	test_switch)

a, b, c = data.getMetaValues()

omega = c[30]

emptypos = [a[len(a) / 2 - 3], b[len(b) / 2 - 3], omega]

fullimage(directory, roi, emptypos[0], emptypos[1], emptypos[2])
fig, ax, atrue, btrue = histogram(a, b, omega, theta)

print btrue

poi = [300, 550]
size = [50, 50]

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]

data.roi = roi
data.bg_combined = data.bg_combined_full[roi[2]:roi[3], roi[0]:roi[1]]
anglepos = [a[len(a) / 2], b[len(b) / 2 + 3], omega]

img = fullimage(directory, roi, anglepos[0], anglepos[1], anglepos[2], title='disl')

for i in range(11):
	apind = [len(a) / 2, len(b) / 2 - 5 + i]
	anglepos = [a[apind[0]], b[apind[1]]]
	line, fig, ax = add_subplot(fig, ax, anglepos[0], anglepos[1], omega, apind, i)
	np.savetxt(directory + '/line_{:06.4f}.txt'.format(btrue[apind[1]]), line)

# subax.imshow(img, cmap='Greens')
# subax.set_axis_off()
# subax.imshow(img)


fig.savefig(directory + '/histogram.pdf')
