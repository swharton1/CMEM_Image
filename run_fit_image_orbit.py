#This will work out the ppmlr image along an orbit and try to fit a CMEM image to it all the way along. It will save the images and create an output file with rmp_sub and the model parameters. 

import CMEM_Image
import numpy as np
import os
import pickle

#Get the PPMLR simulation object. 
ppmlr_filename = "S05D20V400B0000-05rad.fits"
ppmlr = CMEM_Image.ppmlr_fits.read_ppmlr_fits(filename=ppmlr_filename)

#Get the elliptical orbit object. 
#nu = np.arange(360)
#rp=2
#ra=19
#inc=70
#raan=260
#omega=300
#ellipse = CMEM_Image.ellipse.ellipse(nu, rp=rp, ra=ra, inc=inc, raan=raan, omega=omega)

#Make time array to go into orbit simulator using new ellipse code.  
t = np.linspace(0,48,49)
rp=2
ra=19
inc=70
raan=180
omega=300
ellipse = CMEM_Image.ellipse2.ellipse(t, rp=rp, ra=ra, inc=inc, raan=raan, omega=omega)

#Extract key information. 
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
tv = t[view]

#Make a directory to put the images in. 
orbit_folder_name = "ra{}_rp{}_inc{}_raan{}_omega{}/".format(rp, ra, inc, raan, omega) 
plot_path = os.path.join(os.environ.get("PLOT_PATH"), 'fitted_orbit_sim/')
if not os.path.isdir(os.path.join(plot_path, orbit_folder_name)):
    print ('Making folder {}:'.format(os.path.join(plot_path, orbit_folder_name))) 
    os.mkdir(os.path.join(plot_path, orbit_folder_name)) 
if not os.path.isdir(os.path.join(plot_path, orbit_folder_name, ppmlr_filename)):
    print ('Making folder {}:'.format(os.path.join(plot_path, orbit_folder_name, ppmlr_filename))) 
    os.mkdir(os.path.join(plot_path, orbit_folder_name, ppmlr_filename))
full_path = os.path.join(plot_path, orbit_folder_name, ppmlr_filename)

#Set the parameters for the image. 
n_pixels=100
m_pixels=50
target_loc=(9.5,0,0)

#Make lists to return values to. 
output_dict = {}
output_dict['maxIx'] = ppmlr.maxIx
output_dict['maxdIx'] = ppmlr.maxdIx
output_dict['f.25'] = ppmlr.f[0]
output_dict['target_loc'] = target_loc
output_dict['cmem_mp'] = []
output_dict['params'] = []
output_dict['min_cost'] = []
output_dict['t_list'] = tv
output_dict['inout'] = []
output_dict['smile_loc'] = []
output_dict['n_pixels'] = n_pixels
output_dict['m_pixels'] = m_pixels
output_dict['density'] = ppmlr.density 
output_dict['x_gsm'] = xv
output_dict['y_gsm'] = yv
output_dict['z_gsm'] = zv 
output_dict['rp'] = rp
output_dict['ra'] = ra
output_dict['inc'] = inc
output_dict['omega'] = omega
output_dict['raan'] = raan
output_dict['dtimes'] = ellipse.dt_list 

#Loop through every 5/10 points along the viewable part of the orbit and work out the SMILE object and then the image through the simulation. 


for i in range(len(rv)):
    
    print (i)
    smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=(xv[i], yv[i], zv[i]), target_loc=target_loc) 
        
    #To get an image through the ppmlr datacube. 
    #ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile)
        
    #Fit to the image. 
    #fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
    #fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func='absolute') 
    #fit.write_pickle()
        
    #Get the filename of the pickle file. 
    pkl_filename = 'fit_image_n_{:.1f}_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_cmem_absolute_im2_.pkl'.format(ppmlr.density, smile.smile_loc[0], smile.smile_loc[1], smile.smile_loc[2], target_loc[0], target_loc[1], target_loc[2], n_pixels, m_pixels)
        
    figname = os.path.join("fitted_orbit_sim", orbit_folder_name, ppmlr_filename, "fit_image_n_{:.1f}_{:0>2d}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_{}_{}_im{}{}_fov.png".format(ppmlr.density, i, smile.target_loc[0], smile.target_loc[1], smile.target_loc[2], smile.n_pixels, smile.m_pixels, 'cmem', 'absolute', 2, ''))
        
    #Now make the plots to show the output of the fitting process. 
    analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=pkl_filename, model='cmem')
    analysis.plot_images(save=True, add_mp_projection=True, fname=figname, add_fov_projection=True, ellipse=ellipse, los_max=60)
    print ("Image made...") 
        
    #Get subsolar magnetopause out. 
    output_dict['cmem_mp'].append(analysis.rmp_sub)
        
    #Get all parameters out. 
    output_dict['params'].append(analysis.model['params best nm'] ) 
    
    #Get minimum cost out. 
    output_dict['min_cost'].append(analysis.model['min cost'])
    
    #Get info on whether s/c is inside or outside MP.
    output_dict['inout'].append(analysis.inout)
    
    #Get Smile position. 
    output_dict['smile_loc'].append(smile.smile_loc) 
        
#Save the output from the full orbital run in a more convenient file. 
rmp_output_file = "parameter_output.pkl" 
output_fullpath = os.path.join(plot_path, orbit_folder_name, ppmlr_filename, rmp_output_file)
#output_dict = {'rmp': np.array(rmp_sub_list), 'params':np.array(params_list), 'min cost':np.array(min_cost_list), 'nu':np.array(nu_list)}

with open(output_fullpath, 'wb') as f: 
    pickle.dump(output_dict, f)
print ('Pickled: ', output_fullpath)  


