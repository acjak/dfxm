# !/bin/python
"""blah."""
import math
import numpy as np

import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def makeGvector(xylin, ang, phi, w):
    gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

    xcpos = len(gv[0, :, 0])/2
    ycpos = len(gv[0, 0, :])/2

    prefix = 2*math.pi/bv_111
    # prefix = 2*math.sin(ang)

    for i in range(len(xylin)):
        for j in range(len(xylin)):
            if i != xcpos or j != ycpos:
                eps_xz = (bv_110/(2*math.pi))*(xylin[j])/((xylin[i])**2 + (xylin[j])**2)
                eps_yz = (-bv_110/(2*math.pi))*(xylin[i])/((xylin[i])**2 + (xylin[j])**2)
                gv[0, i, j] = prefix*(eps_xz*math.cos(w)-eps_yz*math.sin(phi)*math.sin(w))
                gv[1, i, j] = prefix*(eps_yz*math.cos(phi))
                gv[2, i, i] = prefix*(eps_yz*math.cos(w)*math.sin(phi)+eps_xz*math.sin(w))
                # prefix*bv/(2*math.pi)*(xylin[j])/((xylin[i])**2 + (xylin[j])**2)
                # gv[1, i, j] = -prefix*bv/(2*math.pi)*(xylin[i])/((xylin[i])**2 + (xylin[j])**2)
                # gv[2, i, j] = prefix

    return gv


def makeGvectorDiamondPantleon(bv_111, bv_110, xylin, ang, phi, w):
    gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

    xcpos = len(gv[0, :, 0])/2
    ycpos = len(gv[0, 0, :])/2

    # prefix = 2*math.pi/bv_111
    prefix = 1
    # prefix = 2*math.sin(ang)

    sw = math.sin(w)
    cw = math.cos(w)

    for i in range(len(xylin)):
        for j in range(len(xylin)):
            if i != xcpos or j != ycpos:
                xc = -xylin[i]*sw
                eps_xz = (-bv_110/(2*math.pi))*(xc)/((xylin[j])**2 + (xc)**2)
                eps_yz = (bv_110/(2*math.pi))*(xylin[j])/((xylin[j])**2 + (xc)**2)
                gv[0, i, j] = prefix*(eps_yz*cw-sw)
                gv[1, i, j] = prefix*(eps_xz)
                gv[2, i, i] = prefix*(eps_yz*cw+sw)

    return gv


def makeGvectorDiamondSonja(bv_111, bv_110, xylin, ang, alpha, beta):
    gv = np.zeros((3, len(xylin), len(xylin)))  # G-vector

    xcpos = len(gv[0, :, 0])/2
    ycpos = len(gv[0, 0, :])/2

    sa = math.sin(math.radians(alpha))
    ca = math.cos(math.radians(alpha))
    sb = math.sin(math.radians(beta))
    cb = math.cos(math.radians(beta))

    z = 0

    prefix = (bv_110/bv_111)*ca

    for i in range(len(xylin)):
        for j in range(len(xylin)):
            if i != xcpos or j != ycpos:

                r2p = (ca*(xylin[i]*cb+xylin[j]*sb)-sa*z)**2 + (-xylin[i]*sb+xylin[j]*cb)**2
                gv[0, i, j] = prefix*(1/r2p)*(-sa*sb*z + ca*xylin[j])
                gv[1, i, j] = prefix*(1/r2p)*(sa*cb*z - ca*xylin[i])
                gv[2, i, j] = prefix*(1/r2p)*(sa*(sb*xylin[i] + cb*xylin[j]))+2*math.pi/bv_111

    return gv


def makeQvectorPrecise(omega, xi):
    sinc = np.sin(omega-xi/2.)
    cosc = np.cos(omega-xi/2.)
    sinx = np.sin(xi/2.)
    cosx = np.cos(xi/2.)

    q = np.array([np.swapaxes(2*sinx*sinc, 0, 1), np.zeros((np.shape(np.swapaxes(2*sinx*cosc, 0, 1)))), np.swapaxes(2*sinx*cosc, 0, 1)])

    return q


def makeQvector(omega, xi, dXi, dZeta):
    sinc = np.sin(omega-xi/2.)
    cosc = np.cos(omega-xi/2.)
    sinx = np.sin(xi/2.)
    cosx = np.cos(xi/2.)

    if rank == 0:
        print np.degrees(omega[4])
        print np.degrees(omega[4]-xi[3])
        print np.degrees(omega-xi)[4, 3]
        print np.swapaxes(np.degrees(omega-xi), 0, 1)[4, 3]

    # For x:
    x0 = np.swapaxes(2*(sinx*sinc + cosx*(dXi/2)*sinc + sinx*(dXi/2)*(-cosc)), 0, 1)
    x1 = np.swapaxes(2*(sinx*sinc + cosx*(dXi/2)*sinc - sinx*(dXi/2)*(-cosc)), 0, 1)
    x2 = np.swapaxes(2*(sinx*sinc - cosx*(dXi/2)*sinc + sinx*(dXi/2)*(-cosc)), 0, 1)
    x3 = np.swapaxes(2*(sinx*sinc - cosx*(dXi/2)*sinc - sinx*(dXi/2)*(-cosc)), 0, 1)

    # For y:
    y0 = 2*(sinx*(dZeta))
    y1 = 2*(sinx*(-dZeta))

    # For z:
    z0 = np.swapaxes(2*(sinx*cosc + cosx*(dXi/2)*cosc + sinx*(dXi/2)*(sinc)), 0, 1)
    z1 = np.swapaxes(2*(sinx*cosc + cosx*(dXi/2)*cosc - sinx*(dXi/2)*(sinc)), 0, 1)
    z2 = np.swapaxes(2*(sinx*cosc - cosx*(dXi/2)*cosc + sinx*(dXi/2)*(sinc)), 0, 1)
    z3 = np.swapaxes(2*(sinx*cosc - cosx*(dXi/2)*cosc - sinx*(dXi/2)*(sinc)), 0, 1)

    return x0, x1, x2, x3, y0, y1, z0, z1, z2, z3


def fullMatrixSPIPrecise(wl, gv, omega, xi, local_n):
    matrix_part = np.zeros((local_n, xi_len, len(xylin), len(xylin)))
    dist = np.zeros((3))

    # sigmax = 1E-3
    # # sigmay = 5.9666E-6
    # sigmay = 5.9666E-4
    # sigmaz = 8E-3

    sigma_a = 7.698e-5
    # sigma_a = 2.8e-4

    for m in range(local_n):
        if rank == 0:
            done = 100*float(m)/local_n
            print "Calculation is %g perc. complete..." % done
        for n in range(xi_len):
            q = [(2*math.pi/wl)*2*np.sin(xi[n]/2)*np.cos(omega[m]-xi[n]/2),
                 (2*math.pi/wl)*2*np.sin(xi[n]/2)*np.sin(omega[m]-xi[n]/2),
                 (2*math.pi/wl)*2*np.sin(xi[n]/2)*np.sin(omega[m]-xi[n]/2)]
                 # 0]
            print q

            ql = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)

            sigmax = 2*ql*sigma_a
            sigmay = ql*1e-5
            sigmaz = ql*(1/math.tan(np.radians(omega0+omega[m])))*sigma_a

            matrix_part[m, n, :, :] = 1
            matrix_part[m, n, :, :] *= np.exp(-(abs(gv[0, :, :] - q[0])**2)/(2*sigmax**2))
            matrix_part[m, n, :, :] *= np.exp(-(abs(gv[1, :, :] - q[1])**2)/(2*sigmay**2))
            matrix_part[m, n, :, :] *= np.exp(-(abs(gv[2, :, :] - q[2])**2)/(2*sigmaz**2))

            # for i in range(len(xylin)):
            #     for j in range(len(xylin)):
            #
            #         dist[0] = abs(abs(gv[0, i, j]) - abs(q[0]))
            #         dist[1] = abs(abs(gv[1, i, j]) - abs(q[1]))
            #         dist[2] = abs(abs(gv[2, i, j]) - abs(q[2]))
            #
            #         matrix_part[m, n, i, j] += math.exp(-(dist[1]**2)/(2*sigmay**2))
            #
            #         # for k in range(len(q)):
            #         #    dist[k] = abs(gv[k, i, j] - q[k])
            #         #    print gaus(1., 1., dist[k], -1.)
            #
            #         if rank == 0 and i == 0 and j == 0:
            #             print dist[0]
            #             print dist[1]
            #             if n == m:
            #                 print n, m
            #                 print xi[n], omega[n]
            #                 print math.exp(-((dist[2])**2)/(2*sigmaz**2))
            #                 print q[0], gv[0, i, j]
            #                 print q[1], gv[1, i, j]
            #                 print q[2], gv[2, i, j]
            #                 #     print math.exp(-(dist[2]**2)/(2*sigmaz**2)), dist[2], sigmaz
            #                 print ""

                    # if math.exp(-(dist[2]**2)/(2*sigmaz**2)) != 0:
                    # #     # print math.exp(-(dist[0]**2)/(2*sigmax**2)), dist[0], sigmax
                    # #     # print math.exp(-(dist[1]**2)/(2*sigmay**2)), dist[1], sigmay
                    # #
                    #     print q[2], gv[2, i, j]
                    #     print math.exp(-(dist[2]**2)/(2*sigmaz**2)), dist[2], sigmaz
                    #     print ""

                    # matrix_part[m, n, i, j] += math.exp(-(dist[0]**2)/(2*sigmax**2))

                    # matrix_part[m, n, i, j] *= math.exp(-(dist[2]**2)/(2*sigmaz**2))

                    # if i < 40 and i > 20 and j < 30 and j > 20:
                    #     matrix_part[m, n, i, j] = 0

                    # x, a, x0, sigma
                    # print dist
    matrix_part[0, 0, 0, 0] = rank
    return matrix_part


def fullMatrixSPI(q, gv, omega, local_n):
    matrix_part = np.zeros((local_n, xi_len, len(xylin), len(xylin)))

    for m in range(local_n):
        if rank == 0:
            done = 100*float(m)/local_n
            print "Calculation is %g perc. complete..." % done
        for n in range(xi_len):
            for i in range(len(xylin)):
                for j in range(len(xylin)):
                    if xi_len == 1:
                        qx_max = max(q[0][m], q[1][m], q[2][m], q[3][m])
                        qx_min = min(q[0][m], q[1][m], q[2][m], q[3][m])

                        qy_max = max(q[4], q[5])
                        qy_min = min(q[4], q[5])

                        qz_max = max(q[6][m], q[7][m], q[8][m], q[9][m])
                        qz_min = min(q[6][m], q[7][m], q[8][m], q[9][m])

                    else:
                        # if rank == 1:
                        #     print m, n, np.shape(q[0][:, :]), len(q[0][:, n]), om_len_part
                        qx_max = max(q[0][m, n], q[1][m, n], q[2][m, n], q[3][m, n])
                        qx_min = min(q[0][m, n], q[1][m, n], q[2][m, n], q[3][m, n])

                        qy_max = max(q[4][n], q[5][n])
                        qy_min = min(q[4][n], q[5][n])

                        qz_max = max(q[6][m, n], q[7][m, n], q[8][m, n], q[9][m, n])
                        qz_min = min(q[6][m, n], q[7][m, n], q[8][m, n], q[9][m, n])

                    if gv[0, i, j] <= qx_max and gv[0, i, j] >= qx_min \
                        and gv[1, i, j] <= qy_max and gv[1, i, j] >= qy_min \
                            and gv[2, i, j] <= qz_max and gv[2, i, j] >= qz_min:
                        matrix_part[m, n, i, j] += 5

                        # print m
                        # print n, m*rank, i, j, 'rank=', rank
                        # if rank == 0:
                        #    print q[6:9][m], gv[2, i, j]

                    # if gv[0, i, j] <= qx_max and gv[0, i, j] >= qx_min \
                    #         and gv[1, i, j] <= qy_max and gv[1, i, j] >= qy_min \
                    #             and gv[2, i, j] != qz_max and gv[2, i, j] != qz_min:
                    #         print n, m+rank, i, j, 'rank =', rank
                    #     # print q[6], q[7], q[8], q[9], gv[2, i, j]
                    #     # print 'max= ', max(q[4], q[5])
                    #     # print 'min= ', min(q[4], q[5]), '\n'
                    #     matrix_part[n, m, i, j] += 1
    matrix_part[0, 0, 0, 0] = rank
    return matrix_part


def makeMatrixPart(x0, x1, x2, x3, y0, y1, z0, z1, z2, z3, om_len):
    # Chose part of data set for a specific CPU (rank).
    local_n = int(om_len/float(size))
    istart = rank*local_n
    istop = (rank+1)*local_n

    x0 = x0[istart:istop, :]
    x1 = x1[istart:istop, :]
    x2 = x2[istart:istop, :]
    x3 = x3[istart:istop, :]

    z0 = z0[istart:istop, :]
    z1 = z1[istart:istop, :]
    z2 = z2[istart:istop, :]
    z3 = z3[istart:istop, :]

    local_q = [x0, x1, x2, x3, y0, y1, z0, z1, z2, z3]

    return local_q, istart, istop, local_n


def getSPILimits(size, rank, om_len):
    local_n = int(om_len/float(size))
    istart = rank*local_n
    istop = (rank+1)*local_n

    return istart, istop, local_n


def gaus(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def fitGaussian(x, y):
    from scipy.optimize import curve_fit
    # n = len(x)                          # the number of data
    # mean = sum(x*y)/n                   # note this correction
    # sigma = sum(y*(x-mean)**2)/n        # note this correction
    # sigma = -3.E-3

    try:
        popt, pcov = curve_fit(gaus, x, y, p0=[max(y), x[np.argmax(y)], -5.E-3], maxfev=100000000)
        return popt, pcov
    except RuntimeError:
        print "Error - curve_fit failed"


def getStrainCurve(line, omega, ang, xylin):
    offset_vector = []
    om = np.ogrid[10.9805:10.9995:0.0001]

    fig, ax1 = plt.subplots()
    for i in range(len(xylin)):
        popt, pcov = fitGaussian(np.degrees(omega), line[:, i])
        offset_vector.append(popt[1]/ang)

    ax1.plot(xylin, offset_vector)
    plt.show()


def makeT2TMap(matrix, om_len, xi_len, ang, omega, ptype):
    plt.figure(figsize=(16, 16))
    gs1 = matplotlib.gridspec.GridSpec(xi_len, om_len)
    if ptype == 'img':
        gs1.update(wspace=0.2,  hspace=0.2)

    if ptype == 'line':
        gs1.update(wspace=0.8,  hspace=0.8)

    sns.set_style("white")

    k = 0
    for n in range(xi_len):
        for m in range(om_len):
            axarr = plt.subplot(gs1[k])
            axarr.set_title('%.3f %.3f' % (math.degrees(omega[(om_len-1)-m]-ang), math.degrees(xi[n]-2*ang)))
            # axarr.set_title('%.3f' % np.max(matrix[m, n, :, :]))
            # if np.sum(matrix[m, n, :, :]) != 0:
            #    print m, n

            if ptype == 'img':
                axarr.imshow(matrix[(om_len-1)-m, n, :, :], cmap='Greens', vmin=0, vmax=1)
                axarr.xaxis.set_major_formatter(plt.NullFormatter())
                axarr.yaxis.set_major_formatter(plt.NullFormatter())
                axarr.text(10, 10, np.max(matrix[(om_len-1)-m, n, :, :]))
                # print np.max(matrix[(om_len-1)-m, n, :, :])

            if ptype == 'line':
                axarr.plot(xylin/1000, np.sum(matrix[(om_len-1)-m, n, :, :], axis=1))
                labels = axarr.get_xticks().tolist()
                axarr.set_xticklabels(labels, rotation=-90)
                # axarr.xaxis.set_major_formatter(plt.NullFormatter())
                # axarr.yaxis.set_major_formatter(plt.NullFormatter())

            k += 1

    plt.savefig('output/screwdisl_model_img_array.png')

# bv = 0.192  # 0.3596311555543606   # Burger's vector
bv_111 = 0.35668/math.sqrt(3)   # Burger's vector
bv_110 = 0.35668/math.sqrt(2)   # Burger's vector
# bv = 0.20087
wl = 0.07293  # Wavelength
# ang = 10.99  # Diffraction angle

beta = 0
alpha = 32.5

ptype = 'img'

step = 50.

omega0 = math.degrees(math.asin(wl/(2*bv_111)))
# omega0 = 10.94822
ang_range = 0.02
ang_steps = 10

xylin = np.ogrid[-7500.:7500.:step]

omega = np.radians(np.ogrid[omega0-ang_range/2:omega0+ang_range/2:ang_range/(ang_steps)])
xi = np.radians(np.ogrid[2*(omega0-ang_range/2):2*(omega0+ang_range/2):2*ang_range/(ang_steps)])

#
# if rank == 0:
#     print np.degrees(omega), np.degrees(xi)
#     print len(omega), len(xi)
# omega = np.radians(np.ogrid[10.9805:10.9995:0.002])
# xi = np.radians(np.ogrid[2*10.9805:2*10.9995:2*0.002])
xi = xi[:, np.newaxis]
# xi = np.radians(np.array([2*10.99]))

ang = math.radians(omega0)

dXi = np.radians(0.001)
dZeta = np.radians(0.01)

om_len = len(omega)
xi_len = len(xi)

bv_111 = 0.35668/math.sqrt(3)   # Burger's vector
bv_110 = 0.35668/math.sqrt(2)   # Burger's vector

# Make G vector.
gv = makeGvectorDiamondSonja(bv_111, bv_110, xylin, ang, alpha, beta)

# Make Q vector.
# x0, x1, x2, x3, y0, y1, z0, z1, z2, z3 = makeQvector(omega, xi, dXi, dZeta)

# Divide matrix into parts for each CPU.
# local_q, istart, istop, local_n = makeMatrixPart(x0, x1, x2, x3, y0, y1, z0, z1, z2, z3, om_len)

# Calculate strain on part of data set.
# matrix_part = fullMatrixSPI(local_q, gv, omega[istart:istop], local_n)

### TEST ###
istart, istop, local_n = getSPILimits(size, rank, om_len)
matrix_part = fullMatrixSPIPrecise(wl, gv, omega[istart:istop], xi, local_n)

print "Computation done."

# CPU 0 (rank 0) combines data parts from other CPUs.
if rank == 0:
    # Make empty arrays to fill in data from other cores.
    recv_buffer = np.zeros((local_n, xi_len, len(gv[0, :, 0]), len(gv[0, 0, :])))
    matrix = np.zeros((om_len, xi_len, len(gv[0, :, 0]), len(gv[0, 0, :])))

    datarank = matrix_part[0, 0, 0, 0]
    matrix_part[0, 0, 0, 0] = 0
    matrix[istart:istop, :, :, :] = matrix_part
    for i in range(1,  size):
        comm.Recv(recv_buffer,  MPI.ANY_SOURCE)
        datarank = recv_buffer[0, 0, 0, 0]
        # print datarank*local_n, (datarank+1)*local_n
        # if np.sum(recv_buffer) != 0:
        #     print datarank
        recv_buffer[0, 0, 0, 0] = 0
        matrix[datarank*local_n:(datarank+1)*local_n, :, :, :] = recv_buffer
else:
    # all other processes send their result
    comm.Send(matrix_part)

if rank == 0:
    img = np.zeros((np.shape(matrix[0, 0, :, :])))
    print "All data collected."

    print np.max(matrix[:, :, :, :])
    print np.argmax(matrix[1, 5, :, 0]), np.argmax(matrix[1, 5, 0, :])
    print np.shape(matrix[1, 5, :, :])

    makeT2TMap(matrix, om_len, xi_len, ang, omega, ptype)

    # line = np.sum(matrix[:, 3, :, :], axis=1)
    # getStrainCurve(line, omega, ang, xylin)

    plt.show()

    # q = makeQvectorPrecise(omega, xi)
