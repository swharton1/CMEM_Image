
from . import smile_fov 
from . import boundary_emissivity_functions
from . import sim_image 
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
from . import smile_fov_limb 
from . import smile_fov_limb_old
from . import load_ephemeris_vlocal 
from . import overlays 
from . import orbit_variation 
from . import load_orbit_data

import os 
# This is a backup to set the environment variables if they haven't been set externally. 
#This is mostly so I can do testing in ipython3. 
if "PLOT_PATH" not in os.environ: 
    os.environ["PLOT_PATH"] = "/home/s/sw682/Code/plots/CMEM_Image_plots/"
if "PICKLE_PATH" not in os.environ:
    os.environ["PICKLE_PATH"] = "/data/sol-ionosphere/sw682/pickled_files/CMEM_Image_pickled_models/" 
if "PPMLR_PATH" not in os.environ:
    os.environ["PPMLR_PATH"] = "/data/sol-ionosphere/SMILE/PPMLR/"
if "BATSRUS_PATH" not in os.environ:
    os.environ["BATSRUS_PATH"] = "/data/smile/BATSRUS/"
if "OPENGGCM_PATH" not in os.environ:
    os.environ["OPENGGCM_PATH"] = "/data/smile/OpenGGCM/"
