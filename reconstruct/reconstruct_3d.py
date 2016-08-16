import numpy as n
import fabio
import pylab as pl
import forward_projection
from xfab import tools


###NB
#Detector tilt angles in radians, remaining angle in degrees!
#For use with sin/cos/tan the transformation to radians are done within the relevant functions,
#thus all function input follow the above conventions.

### Sample specs, must eventually be read in from input file
#grain structure simulation (order xyz and smf - slow/medium/fast = upper/lower/Omega)
sg_no = 225
unit_cell = [4.,4.,4.,90.,90.,90.]
hkl = [2,0,0]
wavelength = 0.73
grain_pos = [0,0,0]
grain_dim = [11.,11.,11.]
grain_steps = [11,11,11]
ang_mean = [0,0]
ang_width = [.01,.01]
#image simulation
ang_steps = [2,2,18]
ang_types = ["phi_upper","phi_lower","Omega"]
ang_min = [0,0,-85.]
ang_max = [.01,.01,85.]
detx_size = 20
detz_size = 30
M=1.
t_x = 0
t_y = 0
t_z = 0
stem = "../test/sim"
format = "edf"

#stem = "../test/sim2"
#ang_steps = [2,2,90]
#ang_min = [0,0,-89.]
#ang_max = [.01,.01,89.]

# box simulations
grain_dim = [16.,16.,10.]
grain_steps = [16,16,10]
ang_mean = [0,0]
ang_width = [.01,.01]
ang_steps = [3,3,18]
ang_min = [-.01,-.01,-85.]
ang_max = [.01,.01,85.]
stem = "../test/box"


###



def read_image_stack(stem,ang_steps,format,ang_types):
    """
    Read the stack of images stem_000s_000m_000f.format into a matrix
    data[is,im,if,detx_size,detz_size]
    NB  This will not work for real data, but good for testing
        Eventually add background subtraction here
    """

    data = []
    slow = []
    med = []
    fast = []
    for si in range(ang_steps[0]):
        data.append([])
        for mi in range(ang_steps[1]):
            data[si].append([])
            for fi in range(ang_steps[2]):
                im = fabio.open("%s_%0.4d_%0.4d_%0.4d.%s" %(stem,si,mi,fi,format))
                ss = im.header["%s" %ang_types[0]]
                mm = im.header["%s" %ang_types[1]]
                ff = im.header["%s" %ang_types[2]]
                if si==0 and mi==0:
                    fast.append(eval(ff))
                if mi==0 and fi==0:
                    slow.append(eval(ss))
                if si==0 and fi==0:
                    med.append(eval(mm))
                data[si][mi].append(im.data)

    return n.array(data),slow,med,fast

def read_full_array():
    fullarray = n.load('/u/data/andcj/tmp/largetest.npy')
    fast = n.load('/u/data/andcj/tmp/alpha.npy')
    slow = n.load('/u/data/andcj/tmp/beta.npy')
    med = [0]

    data = []

    for i, slowstep in enumerate(slow):
        data.append([])
        for j, medstep in enumerate(med):
            data[i].append([])
            for k, faststep in enumerate(fast):
                data[i][j].append(fullarray[k, i, :, :])

    return n.array(data), slow, med, fast



def reconstruct(grain_dim,grain_steps,data,slow,med,fast,theta,M,t_x,t_y,t_z,mode="horizontal"):
    """
    Loop through virtual sample voxel-by-voxel and assign orientations based on forward projections onto read image stack.
    Done by finding the max intensity in a probability map prop[slow,med] summed over the fast coordinate.
    """
    grain_xyz = n.zeros(grain_steps+[3])
    grain_ang = n.zeros(grain_steps+[2])
    grain_dimstep = n.array(grain_dim)/n.array(grain_steps)
    grain_prop = n.zeros(grain_steps)

    l = 0
    for iz in range(grain_steps[2]):
        for ix in range(grain_steps[0]):
            for iy in range(grain_steps[1]):
                grain_xyz[ix,iy,iz] = n.array(grain_pos)+grain_dimstep*(n.array([ix,iy,iz])-0.5*(n.array(grain_steps)-1))
                prop = n.zeros((len(slow),len(med)))
                # prop = n.zeros((len(slow)))
                print "voxel", l
                l += 1
                for si in range(len(slow)):
                    for mi in range(len(med)):
                        # for fi in range(len(fast)):
                        xyz_d = forward_projection.forw_proj(grain_xyz[ix,iy,iz],slow[si],med[mi],fast,theta,M,t_x,t_y,t_z,mode)
                        # xyz_d = forward_projection.forw_proj(grain_xyz[ix,iy,iz],slow[si],[0],fast[fi],theta,M,t_x,t_y,t_z,mode)
                        detx = int(round(xyz_d[0]))+n.shape(data)[3]/2
                        detz = int(round(xyz_d[2]))+n.shape(data)[4]/2
                        # print detx, detz
                        #prop[si] = prop[si] + data[si,fi,detx,detz]
                        try:
                            prop[si,mi] = prop[si,mi] + data[si,mi,fi,detx,detz]
                        except IndexError:
                            pass
                            # print "Outside image."
                grain_ang[ix,iy,iz,0] = slow[n.where(prop==n.max(prop))[0][0]]
                grain_ang[ix,iy,iz,1] = med[n.where(prop==n.max(prop))[1][0]]
                grain_prop[ix,iy,iz] = n.max(prop)#/n.sum(prop)


    return grain_xyz,grain_ang,grain_prop





# actual program begins here
# data,slow,med,fast = read_image_stack(stem,ang_steps,format,ang_types)
data,slow,med,fast = read_full_array()
print n.shape(data)
print "slow",slow
print "med",med
print "fast",fast
#pl.imshow(data[0,0,0],interpolation="none")
#pl.show()
tth = tools.degrees(tools.tth(unit_cell,hkl,wavelength))
print "Unit cell:", unit_cell
print "Wavelength: %0.2f A" %wavelength
print "Diffraction angle: %0.2f deg" %tth
theta = 0.5*tth
grain_xyz,grain_ang,grain_prop=reconstruct(grain_dim,grain_steps,data,slow,med,fast,theta,M,t_x,t_y,t_z,mode="horizontal")
pl.figure(0)
pl.imshow(grain_prop[:,:,0],interpolation="none")#,clim=(1/3.,.4))
print n.min(grain_prop),n.max(grain_prop),n.median(grain_prop),n.mean(grain_prop),n.sqrt(n.var(grain_prop))
pl.figure(1)
forward_projection.display_grain_map(grain_ang,ang_max,ang_min,weight=grain_prop,title="recon")
forward_projection.display_grain_map_3d(grain_ang,ang_max,ang_min,weight=grain_prop,title="recon")
