import matplotlib
import matplotlib.pylab as plt
import math
import numpy as np

def makePlots(sampletitle, amplpic, midppic, fwhmpic):
    fig0, ax0 = plt.subplots(1, 1, dpi=500, figsize=[7,7])
    # plt.tight_layout()
    fig1, ax1 = plt.subplots(1, 1, dpi=500, figsize=[7,7])
    # plt.tight_layout()
    fig2, ax2 = plt.subplots(1, 1, dpi=500, figsize=[7,7])
    # plt.tight_layout()

    # fig0.set_figsize_inches(7,7)
    # fig1.set_size_inches(7,7)
    # fig2.set_size_inches(7,7)

    im0 = ax0.imshow(amplpic[3:-3, 3:-3], cmap='jet', interpolation='None')
    im1 = ax1.imshow(fwhmpic[3:-3, 3:-3], cmap='BrBG', interpolation='None')
    im2 = ax2.imshow(midppic[3:-3, 3:-3] , cmap='BrBG', interpolation='None')

    ax0.set_title('AMPL')
    ax1.set_title('FWHM')
    ax2.set_title('MIDP')

    def fmt(x,  pos):
        a,  b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a,  b)

    fig0.subplots_adjust(right=0.8)
    fig1.subplots_adjust(right=0.8)
    fig2.subplots_adjust(right=0.8)

    cbar_ax0 = fig0.add_axes([0.85,  0.1,  0.05,  0.8])
    cbar_ax1 = fig1.add_axes([0.85,  0.1,  0.05,  0.8])
    cbar_ax2 = fig2.add_axes([0.85,  0.1,  0.05,  0.8])

    clb0 = fig0.colorbar(im0, cax=cbar_ax0)  # , format=ticker.FuncFormatter(fmt))
    clb1 = fig1.colorbar(im1, cax=cbar_ax1)  # , format=ticker.FuncFormatter(fmt))
    clb2 = fig2.colorbar(im2, cax=cbar_ax2)  # , format=ticker.FuncFormatter(fmt))

    # # clb0.set_clim(0., 200.)
    # clb2.set_clim(-0.1, 0.1)
    # clb1.set_clim(-0.03, 0.03)

    fig0.savefig('plots/ampl-map_%s.pdf' % (sampletitle))
    fig1.savefig('plots/fwhm-map_%s.pdf' % (sampletitle))
    fig2.savefig('plots/midp-map_%s.pdf' % (sampletitle))
    #return figstrain, axppos

def makeSmallArrays(gaussarray, midprange, fwhmrange, amplrange):
    amplpic = np.zeros((np.shape(gaussarray[:, :, 0])))
    midppic = np.zeros((np.shape(gaussarray[:, :, 0])))
    fwhmpic = np.zeros((np.shape(gaussarray[:, :, 0])))

    amplpic[:, :] = gaussarray[:, :, 0]
    midppic[:, :] = gaussarray[:, :, 1]
    fwhmpic[:, :] = 2*math.sqrt(2*math.log(2))*gaussarray[:, :, 2]

    if midprange != -1:
        midppic[midppic < midprange[0]] = midprange[0]
        midppic[midppic > midprange[1]] = midprange[1]
    if fwhmrange != -1:
        fwhmpic[fwhmpic < fwhmrange[0]] = fwhmrange[0]
        fwhmpic[fwhmpic > fwhmrange[1]] = fwhmrange[1]
    if amplrange != -1:
        amplpic[amplpic < amplrange[0]] = amplrange[0]
        amplpic[amplpic > amplrange[1]] = amplrange[1]

    return amplpic, fwhmpic, midppic


sampletitle = 'rollscan_center_roll'
gaussarray = np.load('reciprocal_data_gauss_fitted/rollscan_center_roll_data.npy')

midprange = [-0.05, 0.05]
fwhmrange = [-0.1, -0.02]
amplrange = [0, 100.]

amplpic, fwhmpic, midppic = makeSmallArrays(gaussarray, midprange, fwhmrange, amplrange)

makePlots(sampletitle, amplpic, midppic, fwhmpic)


sampletitle = 'rollscan_center_tt'
gaussarray = np.load('reciprocal_data_gauss_fitted/rollscan_center_tt_data.npy')

midprange = [-0.05, 0.05]
fwhmrange = -1  # [-0.1, -0.02]
amplrange = [0, 100.]

amplpic, fwhmpic, midppic = makeSmallArrays(gaussarray, midprange, fwhmrange, amplrange)

makePlots(sampletitle, amplpic, midppic, fwhmpic)


sampletitle = 'rollscan_top_tt'
gaussarray = np.load('reciprocal_data_gauss_fitted/rollscan_top_tt_data.npy')

midprange = -1  # [-0.05, 0.05]
fwhmrange = -1  # [-0.1, -0.02]
amplrange = [0, 100.]

amplpic, fwhmpic, midppic = makeSmallArrays(gaussarray, midprange, fwhmrange, amplrange)

makePlots(sampletitle, amplpic, midppic, fwhmpic)


sampletitle = 'rollscan_top_roll'
gaussarray = np.load('reciprocal_data_gauss_fitted/rollscan_top_roll_data.npy')

midprange = [0, 0.1]
fwhmrange = -1  # [-0.1, -0.02]
amplrange = [0, 200.]

amplpic, fwhmpic, midppic = makeSmallArrays(gaussarray, midprange, fwhmrange, amplrange)

makePlots(sampletitle, amplpic, midppic, fwhmpic)
