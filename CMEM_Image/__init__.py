from . import read_ppmlr
from . import smile_fov 
from . import boundary_emissivity_functions
from . import sim_image 
from . import get_names_and_units
from . import set_initial_params
from . import ppmlr_image_old
from . import ppmlr_image
from . import ellipse
from . import ellipse2
from . import fit_model_image_to_ppmlr_image
from . import visualise_image_fit 
from . import analyse_target_variation 
from . import analyse_pixel_variation 
from . import analyse_orbit_variation
from . import analyse_density_variation
from . import analyse_tilt_variation
from . import get_meridians
from . import smile_fov_limb 
from . import smile_fov_limb_old
from . import load_ephemeris_vlocal 
from . import ppmlr_fits
from . import orbit_variation 

import os 
# This is a backup to set the environment variables if they haven't been set externally. 
#This is mostly so I can do testing in ipython3. 
if "PLOT_PATH" not in os.environ: 
    os.environ["PLOT_PATH"] = "/home/s/sw682/Code/plots/CMEM_Image_plots/"
if "PICKLE_PATH" not in os.environ:
    os.environ["PICKLE_PATH"] = "/home/s/sw682/Code/pickled_files/CMEM_Image_pickled_models/" 
if "PPMLR_PATH" not in os.environ:
    os.environ["PPMLR_PATH"] = "/data/sol-ionosphere/SMILE/PPMLR/"
