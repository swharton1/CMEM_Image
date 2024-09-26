#This will work out the CMEM image along the track of an orbit using a range of true anomalies.

import CMEM_Image
import numpy as np
import os

#Get the PPMLR simulation object. 
filename = "S05D05V400B0000-05rad.fits"
ppmlr = CMEM_Image.ppmlr_fits.read_ppmlr_fits(filename=filename)

#Get the elliptical orbit object. 
nu = np.arange(360)
rp=2
ra=19
inc=70
raan=80
omega=300
ellipse = CMEM_Image.ellipse.ellipse(nu, rp=rp, ra=ra, inc=inc, raan=raan, omega=omega)

#Extract key information. 
r = ellipse.r/ellipse.RE
x = ellipse.x3/ellipse.RE
y = ellipse.y3/ellipse.RE
z = ellipse.z3/ellipse.RE

#Only select values with an r value greater than 8. 
view = np.where(r > 8) 
rv = r[view]
xv = x[view]
yv = y[view]
zv = z[view]

#Make a directory to put the images in. 
orbit_folder_name = "ra{}_rp{}_inc{}_raan{}_omega{}/".format(rp, ra, inc, raan, omega) 
plot_path = os.path.join(os.environ.get("PLOT_PATH"), 'ppmlr_orbit_sim/')
if not os.path.isdir(os.path.join(plot_path, orbit_folder_name)):
    print ('Making folder {}:'.format(os.path.join(plot_path, orbit_folder_name))) 
    os.mkdir(os.path.join(plot_path, orbit_folder_name)) 
if not os.path.isdir(os.path.join(plot_path, orbit_folder_name, filename)):
    print ('Making folder {}:'.format(os.path.join(plot_path, orbit_folder_name, filename))) 
    os.mkdir(os.path.join(plot_path, orbit_folder_name, filename))
full_path = os.path.join(plot_path, orbit_folder_name, filename)

#Set the parameters for the image. 
n_pixels=100
m_pixels=50
target_loc=(10,0,0)
#Loop through every 10 points along the viewable part of the orbit and work out the SMILE object and then the image through the simulation. 


for i in range(len(rv)):
    if i%5 == 0:
        print (i)
        smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=(xv[i], yv[i], zv[i]), target_loc=target_loc) 
        
        #To get an image through the ppmlr datacube. 
        ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile)
        
        #Make and save the image. 
        image_name = os.path.join(full_path,'nxm_{}_{}_target_({},{},{})_{:0>2d}'.format(n_pixels, m_pixels, target_loc[0], target_loc[1], target_loc[2], i))
        ppmlr_image.plot_image(los_max=6, ellipse=ellipse, image_name=image_name)
 

