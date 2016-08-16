import numpy as n
from xfab import tools,detector
import pylab as pl
import mpl_toolkits.mplot3d as p3

### 

###

def trans_sample_to_image(phi_up,phi_lo,omega,theta,mode="horizontal"):
    """
    Set up the matrix to transform from the sample to the image system
    Depends on whether we are in horizontal or vertical scattering geometry
    Usage T_s2i = trans_sample_to_image(theta,omega,phi_lo,phi_up,mode="horizontal")
    """
    up = n.pi*phi_up/180.
    lo = n.pi*phi_lo/180.        
    om = n.pi*omega/180.
    th = n.pi*theta/180.
    
    if mode=="horizontal":
        T_up = n.array([[  n.cos(up), -n.sin(up),              0],
                        [  n.sin(up),  n.cos(up),              0],
                        [          0,              0,              1]])
        T_lo = n.array([[  n.cos(lo),              0,  n.sin(lo)],
                        [          0,              1,          0],
                        [ -n.sin(lo),              0,  n.cos(lo)]])
        Omega = n.array([[              1,             0,             0],
                         [              0,  n.cos(om), -n.sin(om)],
                         [              0,  n.sin(om),  n.cos(om)]])
        Theta = n.array([[  n.cos(th), -n.sin(th),              0],
                         [  n.sin(th),  n.cos(th),              0],
                         [          0,          0,              1]])
    else:
        print "error"
        
    T_s2i = n.dot(Theta,n.dot(Omega,n.dot(T_lo,T_up)))
    
    return T_s2i

        
    

def forw_proj(xyz_s,phi_up,phi_lo,omega,theta,M,t_x,t_y,t_z,mode="horizontal",detect="ideal"):
    """
    Calculate the detector coordinate [x_d,0,z_d] hit by the diffracted beam from a given volume element 
    [x_s,y_s,z_s] with a specific orientation [omega,phi_lo,phi_up] at diffraction angle 2theta using a 
    setup of type mode ("horizontal"/"vertical") with magnification M and detector tilt [t_x,t_y,tz] on 
    a detector that is either in the ideal or real system
    Usage: [x_d,0,z_d] = forw_proj([x_s,y_s,z_s],omega,alpha,beta,M,t_x,t_y,t_z,mode="horizontal")
    """
    det_tilt = tools.detect_tilt(t_x,t_y,t_z)
    th = n.pi*theta/180.

    if detect=="ideal":
        det_tilt_frame = n.eye(3)
    elif mode=="horizontal":
        det_tilt_frame = n.array([[  n.cos(th), -n.sin(th),              0],
                                  [  n.sin(th),  n.cos(th),              0],
                                  [          0,          0,              1]])
    elif mode=="vertical":
        det_tilt_frame = n.array([[              1,          0,          0],
                                  [              0,  n.cos(th), -n.sin(th)],
                                  [              0,  n.sin(th),  n.cos(th)]])
                          
    
    xyz_d = -M*n.dot(det_tilt,n.dot(det_tilt_frame,n.dot(trans_sample_to_image(phi_up,phi_lo,omega,theta,mode),n.array(xyz_s))))
    return [xyz_d[0],0,xyz_d[2]]
    
    
def display_grain_map(grain_ang,ang_max,ang_min,weight=None,title=None):
    """
    Function to display a grain map colour coded by two tilt angles
    """
    pl.figure(2)
    slice = n.shape(grain_ang)[2]/2
    map = n.zeros((n.shape(grain_ang)[0],n.shape(grain_ang)[1],3))
    map[:,:,0:2] = (grain_ang[:,:,slice]-ang_min[0:2])/(n.array(ang_max[0:2])-ang_min[0:2]) 
    map[:,:,2] = n.ones((n.shape(grain_ang)[0],n.shape(grain_ang)[1]))-n.maximum.reduce([map[:,:,0],map[:,:,1]])
    if weight != None:
	    map[:,:,0] = map[:,:,0]*weight[:,:,slice]/n.max(weight)
	    map[:,:,1] = map[:,:,1]*weight[:,:,slice]/n.max(weight)
	    map[:,:,2] = map[:,:,2]*weight[:,:,slice]/n.max(weight)
    pl.imshow(map,interpolation="None")
    pl.title(title)
    pl.show()        

def display_grain_map_3d(grain_ang,ang_max,ang_min,weight=None,title=None):
    """
    Function to display a grain map colour coded by two tilt angles
    """
    fig = pl.figure(3,figsize=pl.figaspect(1.0))
    ax = p3.Axes3D(fig)
    map = n.zeros((n.shape(grain_ang)[0],n.shape(grain_ang)[1],n.shape(grain_ang)[2],3))
    map[:,:,:,0:2] = (grain_ang[:,:,:]-ang_min[0:2])/(n.array(ang_max[0:2])-ang_min[0:2]) 
    map[:,:,:,2] = n.ones((n.shape(grain_ang)[0],n.shape(grain_ang)[1],n.shape(grain_ang)[2]))-n.maximum.reduce([map[:,:,:,0],map[:,:,:,1]])
    if weight != None:
	    map[:,:,:,0] = map[:,:,:,0]*weight[:,:,:]/n.max(weight)
	    map[:,:,:,1] = map[:,:,:,1]*weight[:,:,:]/n.max(weight)
	    map[:,:,:,2] = map[:,:,:,2]*weight[:,:,:]/n.max(weight)
    for ix in range(n.shape(grain_ang)[0]):
        for iy in range(n.shape(grain_ang)[1]):
            for iz in range(n.shape(grain_ang)[2]):
                cax=ax.scatter3D(ix,iy,iz,s=100,c=map[ix,iy,iz])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pl.title(title)
    pl.show()        


    