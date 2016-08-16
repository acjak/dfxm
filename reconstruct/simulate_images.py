import numpy as n
from scipy import sparse
import pylab as pl
from xfab import tools
from fabio import edfimage
import forward_projection
import time


###NB
#Detector tilt angles in radians, remaining angle in degrees!
#For use with sin/cos/tan the transformation to radians are done within the relevant functions,
#thus all function input follow the above conventions. 

### Sample specs, must eventually be read in from input file
#grain structure simulation (order xyz and smf - slow/medium/fast = upper/lower/omega)
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
grain_dim = [10.,10.,10.]
grain_steps = [10,10,10]
ang_mean = [0,0]
ang_width = [.01,.01]
ang_steps = [3,3,18]
ang_min = [-.01,-.01,-85.]
ang_max = [.01,.01,85.]
stem = "../test/box"


###


def generate_grain(grain_pos,grain_dim,grain_steps,ang_mean,ang_width):
    """
    Generate a grain matrices: 
    grain_xyz[ix,iy,iz] = [x,y,z] and grain_ang[ix,iy,iz] = [phi_up,phi_lo]
    
    Usage: grain_xyz,grain_ang = generate_grain(grain_pos,grain_dim,grain_steps,ang_mean,ang_width)
    jeto@fysik.dtu.dk 23/06/16
    """
    
    grain_xyz = n.zeros(grain_steps+[3])    
    grain_ang = n.ones(grain_steps+[2])*n.array(ang_mean)
    grain_dimstep = n.array(grain_dim)/n.array(grain_steps)
    
    for ix in range(grain_steps[0]):
        for iy in range(grain_steps[1]):
            for iz in range(grain_steps[2]):
                grain_xyz[ix,iy,iz] = n.array(grain_pos)+grain_dimstep*(n.array([ix,iy,iz])-0.5*(n.array(grain_steps)-1))
                #diagonally striped grain with 3 different orientations
#                choice = (ix+iy+iz)%3
#                ang_add = n.zeros(2)
#                if choice!=2:
#                    ang_add[choice] = ang_width[choice]
#                grain_ang[ix,iy,iz] = ang_mean + ang_add
                #end of diagonal stripes, can eg do random(ly clustered)/continuously varying instead 
                #cluster of orientations
                if ix < grain_steps[0]/2:
                	xsign = -1
                else:
                	xsign = 1
                if iy < grain_steps[1]/2:
                	ysign = -1
                else:
                	ysign = 1
                if iz < grain_steps[2]/2:
                	zsign = -1
                else:
                	zsign = 1
                grain_ang[ix,iy,iz] = ang_mean + zsign*n.array([xsign,ysign])*ang_width
                #end of cluster of orientations
    
    return grain_xyz,grain_ang      


def generate_images(grain_xyz,grain_ang,theta,ang_steps,ang_min,ang_max,M,t_x,t_y,t_z,detx_size,detz_size,mode="horizontal"):
    """
    Generate series of diffractions images sim_000X_000Y_000Z.edf where
    Z is the fastest moving motor (here omega)
    Y is the medium (typically phi_low)
    X is the slow motor (phi_high)
    
    Usage:
    jeto@fysik.dtu.dk 23/06/16
    """
    
    #create image stack as 5D matrix: images[is,im,if,detx_size,detz_size]
    images = n.zeros(ang_steps+[detx_size,detz_size])
                
    #set up angle tolerances to determine which image should be appended the intensity
    grain_steps = n.shape(grain_xyz)
    ang_val_s = n.array(range(ang_steps[0]))*(ang_max[0]-ang_min[0])/(ang_steps[0]-1)+float(ang_min[0])
    ang_val_m = n.array(range(ang_steps[1]))*(ang_max[1]-ang_min[1])/(ang_steps[1]-1)+float(ang_min[1])
    omega_val = n.array(range(ang_steps[2]))*(ang_max[2]-ang_min[2])/(ang_steps[2]-1)+float(ang_min[2])
    ang_tol_s = 0.5*(ang_val_s[1]-ang_val_s[0])
    ang_tol_m = 0.5*(ang_val_m[1]-ang_val_m[0])
    print ang_tol_s,ang_val_s
    print ang_tol_m,ang_val_m
    print omega_val
                
    #loop through sample to add intensities to imagestack
    for ix in range(grain_steps[0]):
        for iy in range(grain_steps[1]):
            for iz in range(grain_steps[2]):
                for io in range(ang_steps[2]):
                    xyz_d = forward_projection.forw_proj(grain_xyz[ix,iy,iz],grain_ang[ix,iy,iz,0],grain_ang[ix,iy,iz,1],omega_val[io],theta,M,t_x,t_y,t_z,mode)
                    for si in range(ang_steps[0]):
                        if abs(grain_ang[ix,iy,iz,0]-ang_val_s[si])<ang_tol_s:
                            slow = si
                            break
                    for mi in range(ang_steps[1]):
                        if abs(grain_ang[ix,iy,iz,1]-ang_val_m[mi])<ang_tol_m:
                            medium = mi
                            break
                    detx = int(round(xyz_d[0]))+detx_size/2
                    detz = int(round(xyz_d[2]))+detz_size/2
                    try:
                        images[slow,medium,io,detx,detz] = images[slow,medium,io,detx,detz] + 1
                    except:
                        pass
                
    return images,ang_val_s,ang_val_m,omega_val 
    
    
def write_edf(slow,medium,omega,ang_val_s,ang_val_m,omega_val,frame,stem,format):
    e=edfimage.edfimage()
    e.data=frame
    e.dim2,e.dim1=frame.shape
    e.header = {}
    e.header['origin']='dfxrm'
    e.header['Dim_1']=e.dim1
    e.header['Dim_2']=e.dim2
    e.header['col_end']=e.dim1-1
    e.header['row_end']=e.dim2-1
    e.header['DataType']='UnsignedShort'
    e.header['Image']=1
    e.header['ByteOrder']='Low'
    e.header['time']=time.asctime()
    e.header['Omega']= omega_val[omega]
    e.header['OmegaStep']=omega_val[1]-omega_val[0]
    e.header['phi_lower']= ang_val_m[medium]
    e.header['phi_lower_step']=ang_val_m[1]-ang_val_m[0]
    e.header['phi_upper']= ang_val_s[slow]
    e.header['phi_upper_step']=ang_val_s[1]-ang_val_s[0]
    e.write('%s_%0.4d_%0.4d_%0.4d.%s' %(stem,slow,medium,omega,format))

    
    
    


#Begin actual program calling the functions defined above
tth = tools.degrees(tools.tth(unit_cell,hkl,wavelength))
print "Unit cell:", unit_cell
print "Wavelength: %0.2f A" %wavelength
print "Diffraction angle: %0.2f deg" %tth
theta = 0.5*tth
grain_xyz,grain_ang = generate_grain(grain_pos,grain_dim,grain_steps,ang_mean,ang_width)
forward_projection.display_grain_map(grain_ang,ang_max,ang_min)
forward_projection.display_grain_map_3d(grain_ang,ang_max,ang_min)
image_stack,ang_val_s,ang_val_m,omega_val = generate_images(grain_xyz,grain_ang,theta,ang_steps,ang_min,ang_max,M,t_x,t_y,t_z,detx_size,detz_size,mode="horizontal")



for si in range(ang_steps[0]): 
    for mi in range(ang_steps[1]): 
        for fi in range(ang_steps[2]): 
            write_edf(si,mi,fi,ang_val_s,ang_val_m,omega_val,image_stack[si,mi,fi],stem,format)
            if image_stack[si,mi,fi].max()>0 and si==0 and mi==0:
                pl.figure(1e8*si+1e4*mi+fi)
                pl.imshow(n.transpose(image_stack[si,mi,fi]),interpolation="none",extent=[1,detx_size,1,detz_size],origin="lower")
                pl.xlabel("x")
                pl.ylabel("z")
                
pl.show()