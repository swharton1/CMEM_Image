#This will work out the ppmlr image along an orbit and try to fit a CMEM image to it all the way along. It will save the images and create an output file with rmp_sub and the model parameters. 

#This file will use the updated smile object to get a realistic pointing direction and read in the positions from Yasir's GSE orbital file. 

import CMEM_Image
import numpy as np
import os
import pickle
import datetime as dt


#Get the orbit data. 
print ('Get orbit data...') 
#stime=(2025,10,1)
#etime=(2025,10,3,3)
stime=(2026,4,2)
etime=(2026,4,4,3)
orbit = CMEM_Image.load_ephemeris_vlocal.orbit(stime=stime, etime=etime, calc_gsm=True)
orbit_data = orbit.new_data 

#Extract key information. 
t = orbit_data['dtime'] 
x_gse = orbit_data['x_gse']
y_gse = orbit_data['y_gse']
z_gse = orbit_data['z_gse']
rv = (x_gse**2 + y_gse**2 + z_gse**2)**0.5 

#Only select values with an r value greater than 8. 
#print ('Only use values above 8 RE')
#view = np.where(r > 8) 
#rv = r[view]
#xv = orbit_data['x_gsm'][view]
#yv = orbit_data['y_gsm'][view]
#zv = orbit_data['z_gsm'][view]
#tv = t[view]

xv = orbit_data['x_gsm']
yv = orbit_data['y_gsm']
zv = orbit_data['z_gsm']
tv = orbit_data['dtime'] 


print ('Number of images:', len(tv)/6)
        

#Get the PPMLR simulation object. 
ppmlr_filename = "S05D20V400B0000-05rad.fits"
ppmlr = CMEM_Image.ppmlr_fits.read_ppmlr_fits(filename=ppmlr_filename)




#Make a directory to put the images in. 
orbit_folder_name = "{}-{}".format(tv[0].strftime("%Y%m%d"), tv[-1].strftime("%Y%m%d")) 
plot_path = os.path.join(os.environ.get("PLOT_PATH"), 'fitted_orbit_sim_limb/')
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


#Make lists to return values to. 
output_dict = {}
output_dict['maxIx'] = ppmlr.maxIx
output_dict['maxdIx'] = ppmlr.maxdIx
output_dict['f.25'] = ppmlr.f[0]
output_dict['target_loc'] = []
output_dict['cmem_mp'] = []
output_dict['params'] = []
output_dict['min_cost'] = []
output_dict['t_list'] = []
output_dict['inout'] = []
output_dict['smile_loc'] = []
output_dict['n_pixels'] = n_pixels
output_dict['m_pixels'] = m_pixels
output_dict['density'] = ppmlr.density 
output_dict['x_gsm'] = []
output_dict['y_gsm'] = []
output_dict['z_gsm'] = [] 
#output_dict['tilt'] = [] 


#Loop through every hour along the viewable part of the orbit and work out the SMILE object and then the image through the simulation. 


for i in range(len(rv)):
    
    print (i%6)
    if i%6 == 0:
        smile = CMEM_Image.smile_fov_limb.smile_limb(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=(xv[i], yv[i], zv[i])) 
        
        #Put constraint conditions here to only calculated when the radial and solar conditions are met. 
        if (smile.radial_constraint is True) & (smile.solar_constraint is True): 
        
        
            #To get an image through the ppmlr datacube. 
            #ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile)
        
            #Get the filename of the pickle file. 
            pkl_filename = 'orbit_lim/fit_image_n_{:.1f}_SMILE_limb_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_cmem_absolute_im2.pkl'.format(ppmlr.density, smile.smile_loc[0], smile.smile_loc[1], smile.smile_loc[2], n_pixels, m_pixels)
        
            #Fit to the image. 
            #fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
            #fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func='absolute') 
            #fit.write_pickle(fname=pkl_filename)
        
        
            figname = os.path.join("fitted_orbit_sim_limb", orbit_folder_name, ppmlr_filename, "fit_image_n_{:.1f}_{:0>3d}_nxm_{}_{}_{}_{}_im{}{}_fov_newsmile.png".format(ppmlr.density, i, smile.n_pixels, smile.m_pixels, 'cmem', 'absolute', 2, ''))
        
            #Now make the plots to show the output of the fitting process. 
            analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=pkl_filename, model='cmem')
            analysis.plot_images(save=True, add_mp_projection=True, fname=figname, add_fov_projection=True, los_max=60)
            print ("Image made...") 
        
            #Add times and positions. 
            output_dict['t_list'].append(tv[i]) 
            output_dict['x_gsm'].append(xv[i])
            output_dict['y_gsm'].append(yv[i])
            output_dict['z_gsm'].append(zv[i]) 
        
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

            #Get Target position. 
            output_dict['target_loc'].append(smile.target_loc) 
        
            #Get Tilt angle. 
            #output_dict['tilt'].append(smile.sxi_tilt)    
            
#Save the output from the full orbital run in a more convenient file. 
rmp_output_file = "orbital_limb_output_newsmile.pkl" 
output_fullpath = os.path.join(plot_path, orbit_folder_name, ppmlr_filename, rmp_output_file)


with open(output_fullpath, 'wb') as f: 
    pickle.dump(output_dict, f)
print ('Pickled: ', output_fullpath)  


