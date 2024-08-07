#This file will run the CMEM_Image module. 
#It's what you would type in ipython3. 
import CMEM_Image
ppmlr = CMEM_Image.read_ppmlr.read_ppmlr_cube()

smile = CMEM_Image.smile_fov.smile_fov(n_pixels=100, m_pixels=50, theta_fov=27, phi_fov=16, smile_loc=(10, -30, 0), target_loc=(6,0,0)) 

image = CMEM_Image.sim_image.image_sim(smile, ppmlr, model="jorg")
image.plot_image() 
