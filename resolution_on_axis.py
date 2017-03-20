# !/bin/python
"""blah."""
# import EdfFile
from lib.getedfdata import GetEdfData
from lib.dfxm import DFXM
import warnings
import numpy as np
import time
import math
# import scipy.interpolate as inter
from scipy import stats
from mpi4py import MPI
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

print rank, mpisize

if rank == 0:
	start = time.time()

# for i in range(len(sys.argv)):
# 	print i, sys.argv[i]

path = '/u/data/andcj/hxrm/Resolution_march_2016/rollscan_big'
bg_path = '/u/data/andcj/hxrm/Resolution_march_2016/rollscan_top_far'

filename = 'scan1_'
bg_filename = 'scan1_0960.edf'
# filename2 = 'ff2_'
sampletitle = 'rollscan_big'
datatype = 'res_paper'
motors = ['samrz', 'detx', 'diffrz']

test_switch = True

poi = [1024, 1024 - 800]
# size = [2, 20]
size = [40, 100]

pixelsize = 180.  # 180.

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]


runvars = [
	path,
	filename,
	bg_path,
	bg_filename,
	datatype,
	roi,
	test_switch,
	motors]

data = GetEdfData(runvars)

try:
	directory = data.directory
except AttributeError:
	directory = 0


tools = DFXM(path, data.data_files, directory, roi, datatype, data.dirhash, data.meta, test_switch)

data.setTest(True)
data.adjustOffset(False)

meta = data.getMetaArray()
a, b, c = data.getMetaValues()


def setroi(poi, size):
	data.roi = [
		poi[0] - size[0] / 2,
		poi[0] + size[0] / 2,
		poi[1] - size[1] / 2,
		poi[1] + size[1] / 2]


def plotImageArray(sampletitle):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025, hspace=0.03)

	for i in range(len(data.imgarray[:, 0, 0])):
		img = data.imgarray[i, :, :]
		axarr = plt.subplot(gs1[i])

		axarr.imshow(img, cmap="Greens")
		axarr.set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	plt.savefig(data.directory + '/%s_array.pdf' % (sampletitle))


def makeIndexList_ROLL(a, b, c):
	index_list = []

	for i in a:
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(i, b[len(b) / 2], c[0])
		index_list.append(index[0])

	return index_list


def makeIndexList_TT(a, b, c):
	index_list = []

	for i in b:
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(a[len(a) / 2], i, c[0])
		index_list.append(index[0])

	return index_list


def getIntList(index_list):
	center_int = np.zeros((len(index_list)))

	srcur = meta[index_list, 3]
	srcur_max = max(srcur)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')

		for i in range(len(data.imgarray)):
			center_int[i] = srcur_max * np.average(data.imgarray[i]) / srcur[i]
			# center_int[i] = np.mean(data.imgarray[i])
		return center_int


def fitMyShitUp(center_int, xr, div):
	popt, pcov = tools.fitGaussian(xr, center_int, div)
	return popt[0], popt[1], popt[2]


def makeGaussian(xr, ampl, midp, sigma):
	xr_fine = np.arange(min(xr), max(xr), abs(xr[1] - xr[0]) / 10)
	yr = []
	for i in xr_fine:
		yr.append(tools.gaus(i, ampl, midp, sigma))
	return xr_fine, yr


def convertToDegrees(b):
	# detx_specular = data.beta0
	# tt0 = 10.89*2
	sample_det = 5720
	b_degrees = []
	for i, detx in enumerate(b):
		b_degrees.append(math.degrees(math.atan(detx / sample_det)))
		print i, detx, math.degrees(math.atan(detx / sample_det))
	return b_degrees


def dtr(angle_degrees):
	return 2 * np.pi * angle_degrees / 360


def nonsmiley():
	beta = []
	ampl1 = np.zeros((100))
	ampl2 = np.zeros((100))
	sig_alpha = 0.0012 / 2.35
	sig_beta = 0.0313 / 2.35 / math.sin(20 * math.pi / 180)
	for ii in range(100):
		beta.append((ii - 50) / 50. * sig_beta * 3)
		alpha = beta[ii] * (1 - math.cos(10 * math.pi / 180))
		ampl1[ii] = math.exp(-0.5 * (beta[ii] / sig_beta)**2) * math.exp(-0.5 * (alpha / sig_alpha)**2)
		ampl2[ii] = math.exp(-0.5 * (beta[ii] / sig_beta)**2)

	max_ampl1 = max(ampl1)
	ampl1 = ampl1 / max_ampl1
	max_ampl2 = max(ampl2)
	ampl2 = ampl2 / max_ampl2
	return np.array(beta), ampl2


def getFits():
	center_int_tt = getIntList(index_list_tt)
	center_int_roll = getIntList(index_list_roll)

	ampl_tt, midp_tt, sigma_tt = fitMyShitUp(center_int_tt, b_d, 3.E-3)
	ampl_roll, midp_roll, sigma_roll = fitMyShitUp(center_int_roll, a, 3.E-3)

	xr_tt, yr_tt = makeGaussian(b_d, ampl_tt, midp_tt, sigma_tt)
	xr_roll, yr_roll = makeGaussian(a, ampl_roll, midp_roll, sigma_roll)

	fwhm_tt = 2 * math.sqrt(2 * math.log(2)) * sigma_tt
	fwhm_roll = 2 * math.sqrt(2 * math.log(2)) * sigma_roll

	tt_results = [ampl_tt, midp_tt, sigma_tt, xr_tt, yr_tt, fwhm_tt]
	roll_results = [ampl_roll, midp_roll, sigma_roll, xr_roll, yr_roll, fwhm_roll]

	return tt_results, roll_results, center_int_tt, center_int_roll


def makePlots(tt_results, roll_results, center_int_tt, center_int_roll, x):
	ampl_tt, midp_tt, sigma_tt, xr_tt, yr_tt, fwhm_tt = tt_results
	ampl_roll, midp_roll, sigma_roll, xr_roll, yr_roll, fwhm_roll = roll_results

	fig, ax = plt.subplots(figsize=(7, 7))
	# fig2, ax2 = plt.subplots(figsize=(7,7))
	ax2 = ax.twiny()

	legend_tt = '2T Ampl:' + str(ampl_tt) + '\n 2T Midp:' + str(midp_tt) + '\n 2T FWHM:' + str(fwhm_tt)
	legend_roll = 'Roll Ampl:' + str(ampl_roll) + '\n Roll Midp:' + str(midp_roll) + '\n Roll FWHM:' + str(fwhm_roll)

	ln1 = ax.plot(xr_tt, yr_tt, label=legend_tt, color='red')
	ax.plot(b_d, center_int_tt, color='black', linestyle=':')
	ax.set_xlabel('2theta [degrees]')

	ln2 = ax2.plot(xr_roll, yr_roll, label=legend_roll)
	ax2.plot(a, center_int_roll, color='black', linestyle=':')
	ax2.set_xlabel('samrz [degrees]')

	lns = ln1 + ln2
	labs = [l.get_label() for l in lns]
	ax.legend(lns, labs, loc=0, prop={'size': 8})
	ax.ticklabel_format(useOffset=False, axis='x')

	# ax.legend(loc='upper right', prop={'size':8})
	# ax2.legend(loc='upper right', prop={'size':8})

	# np.savetxt(data.directory + '/fit_center.txt', np.vstack((xr_tt, yr_tt)).T)
	# np.savetxt(data.directory + '/center.txt', np.vstack((b_d, center_int_tt)).T)
	fig.savefig(data.directory + '/' + sampletitle + '_tt_{:04}.pdf'.format(x))

	# fig3, ax3 = plt.subplots()
	#
	# ax3.imshow(data.imgarray[30], interpolation=None)
	# fig3.savefig(data.directory + '/' + sampletitle + '_image{:04}.pdf'.format(x))


def makeTTPlot(ampl, midp, fwhm, tt_results):
	ampl_tt, midp_tt, sigma_tt, xr_tt, yr_tt, fwhm_tt = tt_results

	midp_tt = 2 * np.pi * midp_tt / 360
	fwhm_tt = 2 * np.pi * fwhm_tt / 360
	xr_tt = 2 * np.pi * xr_tt / 360
	angle = 2 * np.pi * np.array(b_d) / 360

	fig4, ax4 = plt.subplots(1, 2, figsize=(10, 5))

	legend_tt = 'Midp: {:05.3f} mrad\n FWHM: {:05.3} mrad'.format(midp_tt * 1000, fwhm_tt * 1000)
	ax4[0].plot(angle * 1000, center_int_tt, 'ko', ms=3.0)
	lns = ax4[0].plot(xr_tt * 1000, yr_tt, label=legend_tt, color='red')
	ax4[0].set_xlabel(r'$2\Theta$ [mrad]')
	ax4[0].set_ylabel('Intensity [arb. units]')
	ax4[0].yaxis.set_major_formatter(plt.NullFormatter())
	# xlabels = ax4[0].get_xmajorticklabels()
	# xlabels[1::2] = ' ' * len(xlabels[1::2])
	# ax4[0].set_xticklabels(xlabels)
	for label in ax4[0].xaxis.get_ticklabels()[1::2]:
		label.set_visible(False)
	# lns = ln1 + ln2
	labs = [l.get_label() for l in lns]
	ax4[0].legend(lns, labs, loc=0, prop={'size': 10})
	ax4[0].set_ylim(0, max(yr_tt) * 1.15)
	# ax4[0].ticklabel_format(useOffset=False, axis='x')

	pxsteps = 0.086 * (xsteps - 1000) / 4.3 * 4.

	x, y, z = fitMyShitUp(ampl, pxsteps, 20)
	xr, yr = makeGaussian(pxsteps, x, y, z)

	midp = 1000 * 2 * np.pi * (midp - np.mean(midp)) / 360.
	fwhm = 1000 * 2 * np.pi * np.array(fwhm) / 360.

	ax4[1].plot(pxsteps, fwhm, 'r.')
	ax5 = ax4[1].twinx()
	ln1 = ax5.plot(pxsteps, midp, 'b.', label='MIDP data')
	# ax4[1][2].plot(pxsteps, ampl, '.')
	# ax4[1][2].plot(xr, yr, '-')

	ax4[1].set_ylabel('FWHM [mrad]', color='r')
	ax5.set_ylabel('MIDP [mrad]', color='b')
	ax4[1].tick_params('y', colors='r')
	ax5.tick_params('y', colors='b')
	# ax4[1][2].set_title('AMPL')
	ax4[1].set_ylim(np.mean(fwhm) - np.std(fwhm) * 10, np.mean(fwhm) + np.std(fwhm) * 10)
	ax4[1].set_xlabel('Sample position [um]')
	fig4.tight_layout()

	# slope, intercept, r_value, p_value, std_err = stats.linregress(pxsteps, fwhm)
	# legend_fwhm = 'FWHM slope: {}\n FWHM intercept: {}'.format(slope, intercept)
	# ax5.plot(pxsteps, slope * pxsteps + intercept, 'green', label=legend_fwhm)

	slope, intercept, r_value, p_value, std_err = stats.linregress(pxsteps, midp)
	legend_midp = r'MIDP fit slope: %6.3f $m^{-1}$' % (slope * 1000)
	ln2 = ax5.plot(pxsteps, slope * pxsteps + intercept, 'b', label=legend_midp)
	lns = ln1 + ln2
	labs = [l.get_label() for l in ln2]
	ax5.legend(ln2, labs, loc=0, prop={'size': 10})

	fig4.savefig(data.directory + '/' + 'tt_resolution.pdf')

	poi = [1024, 1024]
	size = [2048, 2048]
	setroi(poi, size)

	fig3, ax3 = plt.subplots()
	img = data.getImage(1505, True)
	ax3.imshow(img, interpolation=None)
	fig3.savefig(data.directory + '/' + sampletitle + '_image.pdf')

	fig6, ax6 = plt.subplots()
	ax6.plot(np.sum(img[900:1100, 1100:1200], 0))
	fig6.savefig(data.directory + '/' + sampletitle + '_profile.pdf')


def makeRockRollPlot(roll_results):
	ampl_roll, midp_roll, sigma_roll, xr_roll, yr_roll, fwhm_roll = roll_results

	roll_angle = 2 * np.pi * np.array(a) / 360
	rollfit_angle = 2 * np.pi * xr_roll / 360

	mp = dtr(midp_roll) * 1000
	# xrw, yrw = makeGaussian(rollfit_angle * 1000, ampl_roll, mp, 1.1 / (2 * math.sqrt(2 * math.log(2))))
	xrw, yrw = nonsmiley()
	print xrw
	print yrw

	fig7, ax7 = plt.subplots(1, 2, figsize=(10, 5))

	legend_roll = 'Smiley term'  # .format(dtr(fwhm_roll) * 1000)
	legend_rollwide = 'No smiley term'
	ax7[1].plot(roll_angle * 1000 - mp, center_int_roll, 'ko', ms=3.0)
	ln1 = ax7[1].plot(rollfit_angle * 1000 - mp, yr_roll, label=legend_roll, color='red')
	# ln2 = ax7[1].plot(xrw - mp, yrw, label=legend_rollwide, color='blue')
	ln2 = ax7[1].plot(dtr(xrw) * 1000, yrw * ampl_roll, label=legend_rollwide, color='blue')
	ax7[1].set_xlabel(r'$\beta$ [mrad]')
	ax7[1].set_xlim(-2., 2.)
	# ax7[1].set_xlim(dtr(midp_roll) * 1000 - 2., dtr(midp_roll) * 1000 + 2.)

	lns = ln1 + ln2

	labs = [l.get_label() for l in lns]
	ax7[1].legend(lns, labs, loc=0, prop={'size': 10})

	for label in ax7[1].xaxis.get_ticklabels()[1::2]:
		label.set_visible(False)

	fitvals = np.loadtxt('resolution_paper/output/gaussoutput.txt')
	rock = np.loadtxt('resolution_paper/output/rocking.txt')
	rockfit = np.loadtxt('resolution_paper/output/gaussfit.txt')

	print fitvals

	rock_angle = 2 * np.pi * rock[:, 0] / 360
	rockfit_angle = 2 * np.pi * rockfit[:, 0] / 360
	xrd, yrd = makeGaussian(rockfit_angle * 1000, fitvals[0], dtr(fitvals[1]) * 1000, 1.232e-2 / (2 * math.sqrt(2 * math.log(2))))

	legend_rock = 'Midp: {:03.3f} mrad\n FWHM: {:06.4} mrad'.format(dtr(fitvals[1]) * 1000, dtr(fitvals[2]) * 1000 * 2.35)
	legend_darwin = 'Darwin width of (111) diamond\nFWHM: {:06.4} mrad'.format(1.232e-2)

	ax7[0].plot(rock_angle * 1000, rock[:, 1], 'ko', ms=3.0)
	ln1 = ax7[0].plot(rockfit_angle * 1000, rockfit[:, 1], label=legend_rock, color='red')
	ln2 = ax7[0].plot(xrd, yrd, label=legend_darwin, color='blue')
	ax7[0].set_xlabel(r'$\theta$ [mrad]')
	lns = ln1 + ln2
	labs = [l.get_label() for l in lns]
	ax7[0].legend(lns, labs, loc=0, prop={'size': 10})
	major_formatter = plt.FormatStrFormatter('%2.2f')
	ax7[0].xaxis.set_major_formatter(major_formatter)
	ax7[0].set_ylim(-2, max(yrd) * 1.35)

	ax7[0].set_ylabel('Intensity [arb. units]')
	ax7[1].set_ylabel('Intensity [arb. units]')

	ax7[1].yaxis.set_major_formatter(plt.NullFormatter())
	ax7[0].yaxis.set_major_formatter(plt.NullFormatter())

	fig7.savefig(data.directory + '/' + 'rockandroll.pdf')


index_list_tt = makeIndexList_TT(a, b, c)
index_list_roll = makeIndexList_ROLL(a, b, c)

b_d = convertToDegrees(b)
print b[len(b) / 2], a[len(a) / 2]

xsteps = np.linspace(300., 1700., 30)

ampl = []
midp = []
fwhm = []

for i, x in enumerate(xsteps):
	poi = [x, 1024]
	setroi(poi, size)

	tt_results, roll_results, center_int_tt, center_int_roll = getFits()

	ampl.append(tt_results[0])
	midp.append(tt_results[1])
	fwhm.append(tt_results[5])

	# makePlots(tt_results, roll_results, center_int_tt, center_int_roll, x)

poi = [1024, 1024]
setroi(poi, size)

tt_results, roll_results, center_int_tt, center_int_roll = getFits()

if rank == 0:
	makeTTPlot(ampl, midp, fwhm, tt_results)
	makeRockRollPlot(roll_results)

# fig2.savefig(data.directory + '/' + sampletitle + '_roll_onaxis.pdf')

# plotImageArray(sampletitle)
# tools.plotPPOS(gaussarray)
# tools.plotPPOSBig(gaussarray)
# saveGaussArray(sampletitle, gaussarray)
# tools.plotStrain(strainpic, data.imgarray)
