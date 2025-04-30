# CMEM_Image
This can simulate line of sight intensity images from SMILE through the Jorgensen or CMEM emissivity models. 

All Figures in the CMEM_Image paper can be made with this code. All the code is in the CMEM_Image folder. 

The code has two parts currently. 
First, you need to create a smile object which contains all the geometry on the field of view. 

Then, you create an image simulation object with the desired model parameters and input the smile object. The output is a png image file. 

The code can be ran from ipython3 or the terminal. Use run_cmem_image.py to see what commands you can use. 

THINGS YOU NEED TO DO: 

1. Install SXI_Core first from my github site. Follow those instructions. This contains some core, common functions across all of my projects. 

2. You may want to add CMEM_Image to your PYTHONPATH variable, say in your .bashrc file. Then you can call it from anywhere. Mine looks like:
PYTHONPATH=$PYTHONPATH:~/Code/CMEM_Image/

3. There are a series of paths you will need to set in CMEM/__init__.py. read_fits_cube uses paths stored as environment variables to find the emissivity cubes. You need to set these in the init file. Change whichever of these variables you need too: PPMLR_PATH, OPENGGCM_PATH or BATSRUS_PATH. You will also need to set the default location for your PLOT_PATH and PICKLE_PATH variables, where the output goes. 

