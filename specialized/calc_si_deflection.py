#!/usr/local/bin/python

import os, sys, math
import matplotlib.pylab as plt
import seaborn
import numpy as np

from matplotlib.patches import Rectangle


h = float(sys.argv[1])
dist_from_center = float(sys.argv[2])
l = 20000.
t = 500.

if dist_from_center > t/2:
    print "Position chosen is outside the sample."
    print dist_from_center, "um > ", t/2, "um"

else:
    radius = (l**2 + (4/3) * h**2)/(8*h)
    strain_inner = -(t/2)/radius
    strain_outer = (t/2)/radius
    strain_position = (dist_from_center)/radius
    print "Strain", dist_from_center, "um from center:", strain_position
    print "Top strain: ", strain_inner
    print "Bottom strain: ", strain_outer

    position = np.arange(-t/2,t/2,5.)
    strain = position/radius

    fig,ax = plt.subplots()
    ax.plot(position,strain)

    ax.scatter(dist_from_center,strain_position, color="green", s=200)
    ax.text(dist_from_center,strain_position-0.2*strain_position,('%.2e' % strain_position),ha="center")

    ax.set_ylabel("strain")
    ax.set_xlabel("Distance from center of sample [um]")
    ax.set_xlim(-400,400)

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    someX, someY = 0, 0
    ax.autoscale(False)
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((someX - t/2, someY - .5), t, 1, facecolor="red", alpha = 0.3))
    ax.plot([-t/2,-t/2],[-1,1],color='red',alpha=0.5)
    ax.plot([t/2,t/2],[-1,1],color='red',alpha=0.5)
    plt.show()
