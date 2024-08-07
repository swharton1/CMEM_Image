# CMEM_Image
This can simulate line of sight intensity images from SMILE through the Jorgensen or CMEM emissivity models. 

The code has two parts currently. 
First, you need to create a smile object which contains all the geometry on the field of view. 

Then, you create an image simulation object with the desired model parameters and input the smile object. The output is a png image file. 

The code can be ran from ipython3 or the terminal. Use run_cmem_image.py to see what commands you can use. 


