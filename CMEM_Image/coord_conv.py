#THIS WILL CONTAINS FUNCTIONS TO CONVERT BETWEEN COORDINATE SYSTEMS. 
import numpy as np
#Cartesian-SHUE coordinates. 

def convert_xyz_to_shue_coords(x, y, z):
    '''This will convert the x,y,z coordinates to those used in the Shue model 
     of the magnetopause and bowshock. 

    Parameters
    ---------
    x, y, z - now 3D. Must be entered as numpy arrays.  

    Returns
    -------
    r, theta (rad) and phi (rad)
    '''

   
    # r 
    r = (x**2 + y**2 + z**2)**0.5
       
    # theta - only calc. where coordinate singularities won't occur. 
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] =  np.arccos(x[i]/r[i])

    # phi - only calc. where coordinate singularities won't occur. 
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(z[i], y[i]) 
     
    return r, theta, phi
    
def convert_shue_to_xyz_coords(r, theta, phi):
    '''This will convert the Shue coordinates back to xyz coordinates. 
        
    Parameters
    ----------
    r, theta (rad), phi (rad)
        
    Returns
    -------
    x,y,z
    '''

    x = r*np.cos(theta)
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(theta)*np.sin(phi)

    return x,y,z 



