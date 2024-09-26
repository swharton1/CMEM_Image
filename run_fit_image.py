#This file will run the code to fit a model image to a PPMLR image. 
import numpy as np
import CMEM_Image 

#To use same coordinates at Samsonov et al. (2022b)
#smile_loc = (6.57, -5.94, 17.33)
#target_loc = (9.5,0,0)

#Get the SMILE FOV Object. 
n_pixels=100
m_pixels=50
smile_loc=(6.57,-5.94,17.33)
target_loc=(5,0,5)
smile = CMEM_Image.smile_fov.smile_fov(n_pixels=n_pixels, m_pixels=m_pixels, theta_fov=27, phi_fov=16, smile_loc=smile_loc, target_loc=target_loc) 

#Get the PPMLR Object. 
ppmlr = CMEM_Image.ppmlr_fits.read_ppmlr_fits(filename="S05D20V400B0000-05rad.fits")

#To get an image through the ppmlr datacube. 
ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile) 

#Now run the code to fit a model image to the PPMLR image. 
fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func='absolute') 
fit.write_pickle()

#Use this code if you want to plot the output. 

#First, you need the pickled filename.
filename = 'fit_image_n_20.0_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_cmem_absolute_im2_.pkl'.format(smile_loc[0], smile_loc[1], smile_loc[2], target_loc[0], target_loc[1], target_loc[2], n_pixels, m_pixels)
analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=filename, model='cmem')

#To plot the parameter and cost variation with iteration.
analysis.plot_change_in_parameters(save=True)

#To plot the ppmlr image next to the fitted model image. 
analysis.plot_images(save=True, los_max=60, add_fov_projection=True)
analysis.plot_images(save=True, los_max=60, add_fov_projection=False) 
