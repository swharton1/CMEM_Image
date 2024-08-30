#This will make an elliptical orbit in GSM and calculate a CMEM image for each one. 
#The images will be stored and made into a gif. 

import CMEM_Image
import numpy as np
import os
import matplotlib.pyplot as plt 

#Make time array to go into orbit simulator. 
t = np.linspace(0,48,97)
rp=2
ra=19
inc=70
raan=180
omega=300
ellipse = CMEM_Image.ellipse2.ellipse(t, rp=rp, ra=ra, inc=inc, raan=raan, omega=omega)

#Extract key information. GSM coordinates. 
x = ellipse.coords_gsm.x/ellipse.RE
y = ellipse.coords_gsm.y/ellipse.RE
z = ellipse.coords_gsm.z/ellipse.RE
r = (x**2 + y**2 + z**2)**0.5 

#Only select values with an r value greater than 8. 
view = np.where(r > 8) 
rv = r[view]
xv = x[view]
yv = y[view]
zv = z[view]
print ('xv.size = ', xv.size)

sim_folder = 'sim3/'
#Make a directory to put the images in. 
#orbit_folder_name = "ra{}_rp{}_inc{}_raan{}_omega{}/".format(rp, ra, inc, raan, omega) 
plot_path = os.path.join(os.environ.get("PLOT_PATH"), 'outreach_sim/')
#if not os.path.isdir(os.path.join(plot_path, orbit_folder_name)):
#	print ('Making folder {}:'.format(os.path.join(plot_path, orbit_folder_name))) 
#	os.mkdir(os.path.join(plot_path, orbit_folder_name)) 
if not os.path.isdir(os.path.join(plot_path, sim_folder)):
	print ('Making folder {}:'.format(os.path.join(plot_path, sim_folder))) 
	os.mkdir(os.path.join(plot_path, sim_folder))
full_path = os.path.join(plot_path, sim_folder)

#Set the parameters for the image. 
n_pixels=100
m_pixels=50
target_loc=(9,0,0)

#Loop through all points on the orbit. 
for i in range(len(rv)):
	print ('Image = ',i)
	smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=(xv[i], yv[i], zv[i]), target_loc=target_loc)
	
	#Create CMEM image. 
	image = CMEM_Image.sim_image.image_sim(smile, model='cmem', density=20)
	
	#If you want to plot the whole space. 
	image.get_emissivity_all_coords()
	
	#Use outreach_plot1 for just FOV. outreach_plot2 for whole space. 
	image.outreach_plot2(colour_cap=2, elev=45, azim=45, cmap='bone', fname=full_path+'sim2_{:0>2d}.png'.format(i), max_alpha=0.3)
	
	plt.close("all")


	

