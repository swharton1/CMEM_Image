#This file will run the CMEM_Image module. 
#It's what you would type in ipython3. 
import CMEM_Image

smile = CMEM_Image.smile_fov.smile_fov(n_pixels=100, m_pixels=50, theta_fov=27, phi_fov=16, smile_loc=(10, -30, 0), target_loc=(6,0,0)) 

image = CMEM_Image.sim_image.image_sim(smile, model="cmem", init_method=2, params0=None, temp=200000, density=5, vx=400, vy=0, vz=0, bx=0, by=0, bz=5, dipole=0)
image.plot_image() 
