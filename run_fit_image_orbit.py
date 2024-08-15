#This will work out the ppmlr image along an orbit and try to fit a CMEM image to it all the way along. It will save the images and create an output file with rmp_sub and the model parameters. 

import CMEM_Image
import numpy as np
import os
import pickle

#Get the PPMLR simulation object. 
filename = "S05D05V400B0000-05rad.dat"
ppmlr = CMEM_Image.read_ppmlr.read_ppmlr_cube(filename=filename)

#Get the elliptical orbit object. 
nu = np.arange(360)
rp=2
ra=19
inc=70
raan=260
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
nuv = nu[view]

#Make a directory to put the images in. 
orbit_folder_name = "ra{}_rp{}_inc{}_raan{}_omega{}/".format(rp, ra, inc, raan, omega) 
plot_path = os.path.join(os.environ.get("PLOT_PATH"), 'fitted_orbit_sim/')
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

#Make lists to return values to. 
rmp_sub_list = []
params_list = []
min_cost_list = []
nu_list = []
#Loop through every 5/10 points along the viewable part of the orbit and work out the SMILE object and then the image through the simulation. 


for i in range(len(rv)):
	if i%5 == 0:
		print (i)
		smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=(xv[i], yv[i], zv[i]), target_loc=target_loc) 
		
		#To get an image through the ppmlr datacube. 
		ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile)
		
		#Fit to the image. 
		fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
		fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func='absolute') 
		fit.write_pickle()
		
		#Get the filename of the pickle file. 
		fname = "fit_image_n_{}_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{}_{}_{}_nxm_{}_{}_{}_{}_im{}_{}.pkl".format(5.0, smile.smile_loc[0], smile.smile_loc[1], smile.smile_loc[2], smile.target_loc[0], smile.target_loc[1], smile.target_loc[2], smile.n_pixels, smile.m_pixels, 'cmem', 'absolute', 2, '')
		
		figname = os.path.join("fitted_orbit_sim", orbit_folder_name, filename, "fit_image_n_{}_{:03d}_Target_{}_{}_{}_nxm_{}_{}_{}_{}_im{}{}_fov.png".format(5.0, i, smile.target_loc[0], smile.target_loc[1], smile.target_loc[2], smile.n_pixels, smile.m_pixels, 'cmem', 'absolute', 2, ''))
		
		#Now make the plots to show the output of the fitting process. 
		analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=fname, model='cmem')
		analysis.plot_images(save=True, add_mp_projection=True, fname=figname, add_fov_projection=True, ellipse=ellipse)
		print ("Image made...") 
		
		#Get subsolar magnetopause out. 
		rmp_sub = analysis.rmp_sub 
		rmp_sub_list.append(rmp_sub)
		
		#Get all parameters out. 
		params = analysis.model['params best nm'] 
		params_list.append(params) 
		
		#Get minimum cost out. 
		min_cost = analysis.model['min cost']
		min_cost_list.append(min_cost)
		
		#Get true anomaly for this point on the orbit. 
		nu_i = nuv[i] 
		nu_list.append(nu_i) 
		
#Save the output from the full orbital run in a more convenient file. 
rmp_output_file = "parameter_output.pkl" 
output_fullpath = os.path.join(plot_path, orbit_folder_name, filename, rmp_output_file)
output_dict = {'rmp': np.array(rmp_sub_list), 'params':np.array(params_list), 'min cost':np.array(min_cost_list), 'nu':np.array(nu_list)}

with open(output_fullpath, 'wb') as f: 
	pickle.dump(output_dict, f)
print ('Pickled: ', output_fullpath)  


