#This will run from a constant spacecraft position but will vary the number of pixels to 
#see how this affects the determination of the magnetopause radius. 
#One experiment for my proposed CMEM Image paper. 

import numpy as np 
import CMEM_Image
import os
import pickle 

#Get main plot path. 
plot_path = os.environ.get('PLOT_PATH') 

#Get the SMILE FOV Object. 

smile_loc=(6.57,-5.94,17.33)
target_loc=(9.5,0,0)

#Get range of pixel numbers. 
n_pixels = [10,20,30,40,50,60,70,80,90,100]
m_pixels = [5,10,15,20,25,30,35,40,45,50]
print (n_pixels)
print (m_pixels)

cost_func='normalised'

ppmlr_files = ["S05D05V400B0000-05rad.fits",
"S05D7.5V400B0000-05rad.fits",
"S05D12.25V400B0000-05rad.fits",
"S05D20V400B0000-05rad.fits",
"S05D25V400B0000-05rad.fits",
"S05D35V400B0000-05rad.fits"]

#Set up loop to go through all six PPMLR cubes. 
for ecube in ppmlr_files:

    #Get the PPMLR Object. 
    ppmlr = CMEM_Image.ppmlr_fits.read_ppmlr_fits(filename=ecube)

    #Make a dictionary to save all the key information so you can 
    #complete the analysis. 
    target_dict = {}
    target_dict['maxIx'] = ppmlr.maxIx
    target_dict['maxdIx'] = ppmlr.maxdIx
    target_dict['f.25'] = ppmlr.f
    target_dict['target_loc'] = target_loc
    target_dict['cmem_mp'] = []
    target_dict['inout'] = []
    target_dict['params'] = []
    target_dict['smile_loc'] = smile_loc
    target_dict['n_pixels'] = n_pixels
    target_dict['m_pixels'] = m_pixels
    target_dict['density'] = ppmlr.density 
    target_dict['cost_func'] = cost_func 


    #Loop through each set of pixels. 
    for p in range(len(n_pixels)): 

        print (p)
    
        #Get the SMILE Object. 
        smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels[p], m_pixels=m_pixels[p], theta_fov=27, phi_fov=16, smile_loc=smile_loc, target_loc=target_loc) 
    
        #To get an image through the ppmlr datacube. 
        ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile) 
    
        #First, you need the pickled filename.
        pickle_filename = 'fit_image_n_{:.1f}_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_cmem_{}_im2_.pkl'.format(ppmlr.density, *smile_loc, *target_loc, n_pixels[p], m_pixels[p], cost_func)
    
        #Now run the code to fit a model image to the PPMLR image. 
        fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
        fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func=cost_func) 
        fit.write_pickle(fname=pickle_filename)
    
        figname ='pixel_variation/Pixel_var_n_{:.1f}_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}.png'.format(ppmlr.density, *smile_loc, *target_loc, n_pixels[p], m_pixels[p])
    
        analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=pickle_filename, model='cmem', limb=False)

        #To plot the ppmlr image next to the fitted model image. 
        analysis.plot_images(save=True, los_max=60, add_fov_projection=True, fname=figname)
    
        #Start adding key information to the target_dict from each run. 
        #target_dict['target_x'].append(t)
        target_dict['cmem_mp'].append(analysis.rmp_sub)
        target_dict['inout'].append(analysis.inout)
        target_dict['params'].append(analysis.model['params best nm'])

    #Now save the dictionary as a pickle file. 
    target_output_file = "pixel_variation_output_n_{:.1f}.pkl".format(ppmlr.density)  
    output_fullpath = os.path.join(plot_path, 'pixel_variation/', target_output_file)

    with open(output_fullpath, 'wb') as f: 
        pickle.dump(target_dict, f)
    print ('Pickled: ', output_fullpath)  
