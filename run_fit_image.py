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
