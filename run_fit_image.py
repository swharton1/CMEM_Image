#This file will run the code to fit a model image to a PPMLR image. 
import numpy as np
import CMEM_Image 

#Get the SMILE FOV Object. 
smile = CMEM_Image.smile_fov.smile_fov(n_pixels=100, m_pixels=50, theta_fov=27, phi_fov=16, smile_loc=(10, -30, 0), target_loc=(10,0,0)) 

#Get the PPMLR Object. 
ppmlr = CMEM_Image.read_ppmlr.read_ppmlr_cube(filename="S05D05V400B0000-05rad.dat")

#To get an image through the ppmlr datacube. 
ppmlr_image = CMEM_Image.ppmlr_image.ppmlr_image(ppmlr, smile) 

#Now run the code to fit a model image to the PPMLR image. 
fit = CMEM_Image.fit_model_image_to_ppmlr_image.fit_image(ppmlr_image, smile) 
fit.fit_function_with_nelder_mead(model='cmem', init_method=2, params0=None, cost_func='absolute') 
fit.write_pickle()

#Use this code if you want to plot the output. 

#First, you need the pickled filename.
filename = 'fit_image_n_5.0_SMILE_10_-30_0_Target_10_0_0_nxm_100_50_cmem_absolute_im2_.pkl'
analysis = CMEM_Image.visualise_image_fit.analyse_fit(filename=filename, model='cmem')

#To plot the parameter and cost variation with iteration.
analysis.plot_change_in_parameters(save=True)

#To plot the ppmlr image next to the fitted model image. 
analysis.plot_images(save=True, los_max=12) 
