#This will attempt to fit a model image from CMEM to a PPMLR image for a given SMILE FOV. 

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from time import process_time
import pickle
import os

from . import set_initial_params as sip 
from . import boundary_emissivity_functions as bef 

class fit_image():
	'''This class will try to fit a model image to a PPMLR image. '''
	
	def __init__(self, ppmlr_image, smile):
		'''This takes in the ppmlr_image and smile objects and attaches them. It also works out the shue coordinates.
		
		Parameters
		----------
		ppmlr_image - ppmlr_image object. must be made using the smile_fov object
		smile - smile_fov object. 
		
		'''
		
		self.ppmlr_image = ppmlr_image 
		self.smile = smile
		
		# Extract any useful solar wind parameters
		self.temp = ppmlr_image.temp
		self.density = ppmlr_image.density
		self.vx = ppmlr_image.vx
		self.vy = ppmlr_image.vy
		self.vz = ppmlr_image.vz
		self.bx = ppmlr_image.bx
		self.by = ppmlr_image.by
		self.bz = ppmlr_image.bz
		self.pdyn = self.calc_dynamic_pressure()
		self.pmag = self.calc_magnetic_pressure()
		self.dipole = 0 
		
		# Get the r, theta and phi coordinates. 
		# Convert to Shue cooords. to calculate the function. 
		print ("Calculating shue coordinates:")
		ts = process_time()
		self.r, self.theta, self.phi = self.convert_xyz_to_shue_coords(self.smile.xpos, self.smile.ypos, self.smile.zpos)
		te = process_time()
		print ("Time = {:.1f}s".format(te-ts))   
		
		
	def __repr__(self):
		# This is what will be printed about the model 
		#if you print the object to the terminal. 
		# It's meant to be better than using __string__(self):
		return f"fit image object. Will fit a model image ('jorg' or 'cmem') to a ppmlr image.  "
    
	def calc_dynamic_pressure(self):
		'''Calculate this as it's a parameter in some models.'''

		# Assume average ion mass = mass of proton. 
		mp = 0.00000000000000000000000000167
        
		# Calculate v in m/s 
		v = (self.vx**2 + self.vy**2 + self.vz**2)**0.5
		v = v*1000 

		# Convert number of particles /cm^3 to /m^3. 
		n = self.density*1000000

		# Calculate dynamic pressure first in Pascals, then nPa. 
		dyn_pressure = 0.5*mp*n*(v**2)
		dyn_pressure = dyn_pressure*1000000000

		return (dyn_pressure)

	def calc_magnetic_pressure(self):
		'''Calculate the magnetic pressure'''

		# Calculate magnitude of B in T. 
		B = (self.bx**2 + self.by**2 + self.bz**2)**0.5
		B = B*0.000000001

		# mu0
		mu0 = 4*np.pi*0.0000001

		# Calculate magnetic pressure in Pa, then nPa. 
		mag_pressure = (B**2)/(2*mu0)
		mag_pressure = mag_pressure*1000000000

		return mag_pressure
	
	def convert_xyz_to_shue_coords(self, x, y, z):
		'''This will convert the x,y,z coordinates to those used in the Shue model 
		of the magnetopause and bowshock. 

		Parameters
		----------
		x, y, z - now 3D.  

		Returns
		-------
		r, theta (rad) and phi (rad)
		'''

		# r 
		r = (x**2 + y**2 + z**2)**0.5
        
		# theta - only calc. where coordinate singularities won't occur. 
		theta = np.zeros(r.shape)
		i = np.where(r != 0)
		theta[i] =  np.arccos(x[i]/r[i])

		# phi - only calc. where coordinate singularities won't occur. 
		phi = np.zeros(r.shape)
		j = np.where((y**2 + z**2) != 0)
		phi[j] = np.arccos(y[j]/((y[j]**2 + z[j]**2)**0.5))
        
		return r, theta, phi
		
	#COST FUNCTIONS
	###############

	def get_cost_function(self):
		'''This returns the cost function that calculates the misfit/n.
        
		Parameters
		----------
		self - variable that contains the data.  
        
		Returns
		-------
		Cost Function. 
			- if self.cost_func == "sum squares", it will return the cost function using squared deviations.  
			- elif self.cost_func == "absolute", it will return the cost function using absolute deviations. 
         
		'''

		# To minimise with Nelder-Mead, you need a cost 
		#function, or as Matt J. called it, the misfit 
		#function. See his FitPlasma.py file in WaveHarmonics. 

		if self.cost_func.lower() == "sum squares":
			def cost_func_sum_squared_differences_by_n(params):
				# This cost function is the sum of 
				#the squared differences divided by n. 
                
				# Calculate Model values of eta with 
				#a given set of parameters.
				eta_model = self.get_eta_model(params) 
                
                #Then calculate the LOS intensity image. 
                #Calculate the LOS intensity. 
				los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
		
				# For each pixel: 
				for i in range(self.smile.n_pixels):
					for j in range(self.smile.m_pixels):
			
						#Added unit conversion factor from ev.RE to kev.cm
						los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, eta_model[i][j]))*637100
							
				# Now get the misfit, (model - observed)**2
				sq_diff = (los_intensity - self.ppmlr_image.los_intensity)**2
				cost = sq_diff.sum()/self.ppmlr_image.los_intensity.size
				self.cost_per_iteration.append(cost)
				self.param_list.append(params)
				print (cost)
               
				return cost
			return cost_func_sum_squared_differences_by_n
        
		elif self.cost_func.lower() == "absolute":
			def cost_func_sum_absolute_deviations_by_n(params):
				# This cost function is the sum of 
				#the absolute deviations divided by n. 
				# This is the cost function used in 
				#Jorgensen et al. (2019b).

				# Calculate Model values of eta with a 
				#given set of parameters.
				eta_model = self.get_eta_model(params) 
				
				#Then calculate the LOS intensity image. 
                #Calculate the LOS intensity. 
				los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
		
				# For each pixel: 
				for i in range(self.smile.n_pixels):
					for j in range(self.smile.m_pixels):
			
						#Added unit conversion factor from ev.RE to kev.cm
						los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, eta_model[i][j]))*637100

				#  Now get the misfit, abs(model - observed)
				abs_diff = abs(los_intensity - self.ppmlr_image.los_intensity)
				cost = abs_diff.sum()/self.ppmlr_image.los_intensity.size
				self.cost_per_iteration.append(cost)
				self.param_list.append(params)
				print (cost)

				return cost
			return cost_func_sum_absolute_deviations_by_n
            
		elif self.cost_func.lower() == "normalised":
			def cost_func_sum_squares_by_sum_observed(params):
				# This cost function is the sum of the 
				#squared deviations normalised by the data value
				# and the sum of the observed emissivities. 

				# Calculate Model values of eta with 
				#a given set of parameters.
				eta_model = self.get_eta_model(params) 
                
                #Then calculate the LOS intensity image. 
                #Calculate the LOS intensity. 
				los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
		
				# For each pixel: 
				for i in range(self.smile.n_pixels):
					for j in range(self.smile.m_pixels):
			
						#Added unit conversion factor from ev.RE to kev.cm
						los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, eta_model[i][j]))*637100
						
				#  Now get the misfit, (model - observed)**2/observed**2 
				sq_diff_norm = (los_intensity - self.ppmlr_image.los_intensity)**2
				cost = sq_diff_norm.sum()/(self.ppmlr_image.los_intensity**2).sum()
				self.cost_per_iteration.append(cost)
				self.param_list.append(params)
				print (cost)

				return cost
			return cost_func_sum_squares_by_sum_observed
		else:
			raise ValueError("Invalid cost function chosen. Select either 'sum squares', 'absolute' or 'normalised'.") 


	def get_eta_model(self, params):
		'''This function calculates the eta model values for each iteration. This function is intended to be run from the cost function. 
		Parameters
		----------
		params - tuple of the model parameters for the chosen model. '''
		
		if self.current_model == "jorg": 
			eta_model = self.current_func(self.r, self.theta, self.phi, *params)
		elif self.current_model == "cmem":
			eta_model = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *params)
		else: 
			raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
		
		return eta_model
    	
	def trapezium_rule(self, p_spacing, eta_LOS):
		'''This will integrate a function using the trapezium rule. '''

		return (p_spacing/2)*(eta_LOS[0] + eta_LOS[-1] + 2*sum(eta_LOS[1:-1]))
   
	#FITTING THE MODEL TO THE DATA
	##############################
		   
	def fit_function_with_nelder_mead(self, model = "cmem", params0 = None, cost_func="absolute", init_method=2):
		'''This uses a Nelder-Mead minimisation technique to find the best 
		parameters for the chosen model to the data. 
        
		Parameters
		----------
		model - which model to fit. "jorg" or "cmem" 
		params - (a0, b0, ...) - Tuple containing initial guesses for the model parameter.
			def = None. Will default to values inside the program for each model unless specified.  
		cost - Type of cost function to use. "sum squares" (def), "absolute" or "normalised"
			- Sum Squares will calculate the sum of the squared deviations/n
			- Absolute will calculate the sum of the absolute deviations/n 
			- Normalised will calculate the sum of the (squared deviations) /(n*sum observed)
		init_method - Boolean to use either method 1 or method 2 from the CMEM paper to set the initial model parameters. 

		'''

        
		# Use this tutorial: https://machinelearningmastery.com/how-to-use-nelder-mead-optimization-in-python/#:~:text=The%20Nelder%2DMead%20optimization%20algorithm%20can%20be
		#%20used%20in%20Python,initial%20point%20for%20the%20search.
		# Nelder-Mead does not use gradient methods to find the best fit. 
		# It requires a starting point. 
		# It can be used for multidimensional functions with multiple parameters. 
		# It can be applied to multi-modal functions. 
		# The correct reference is at: https://academic.oup.com/comjnl/article-abstract/7/4/308/354237?login=true

		# First, select the correct model to fit. 
		self.current_model = model.lower() 
		self.current_func = bef.get_model_func(self.current_model)
		self.init_method = init_method

		#GET THE INITIAL PARAMETERS (IF NOT GIVEN)
		##########################################
		
		# Get Lin model coefficients. 
		if self.current_model == "cmem":
			self.lin_coeffs = bef.get_lin_coeffs(self.dipole, self.pdyn, self.pmag, self.bz)
			self.r0_lin = self.lin_coeffs[-1] 
        	
		#Get initial parameters. 
		if params0 is None: 
			if self.current_model == "jorg":
        		
				self.params0 = sip.get_init_params(self.current_model, self.init_method, self.bz, self.pdyn, self.density) 
        
			elif self.current_model == "cmem":
				self.params0 = sip.get_init_params(self.current_model, self.init_method, self.bz, self.pdyn, self.density, self.r0_lin) 
		
			else:
				raise ValueError("{} not a valid model. Choose 'cmem' or 'jorg'".format(self.current_model))
		else:
			self.params0 = params0 
		
		#GET COST FUNCTION AND MINIMISE IT.
		###################################

		# Get cost calculation function. 
		self.cost_func = cost_func.lower()
		Calculate_cost = self.get_cost_function()
        
		# Set arrays to record info as optimisation progresses. 
		self.cost_per_iteration = []
		self.param_list = []
		
		# The minimize function takes the function and 
		#initial parameters to search.
		# There is an option to add in boundaries for each 
		#parameter as (min, max) pairs.
		# It returns an OptimizeResult object (see scipy docs). 
		print ("Minimising function:")
		ts = process_time() 
		self.result = minimize(Calculate_cost, self.params0, method='nelder-mead', bounds=None)
		te = process_time()
		self.opt_time = te-ts
		print ("Time = {:.1f}s".format(self.opt_time)) 

		# Extract the best fit parameters of a and b, 
		#as well as the final cost and number of iterations. 
		self.params_best_nm = self.result.x
		self.minimum_cost = self.result.fun 
        
		# This is not the number of times it went 
		#through the cost function... 
		self.iterations = self.result.nfev 
		self.cost_per_iteration = np.array(self.cost_per_iteration)
	
		#Calculate the final los intensity. 
		self.eta_model = self.get_eta_model(self.params_best_nm) 
		
		#Then calculate the LOS intensity image. 
        #Calculate the LOS intensity. 
		self.model_los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
		
		# For each pixel: 
		for i in range(self.smile.n_pixels):
			for j in range(self.smile.m_pixels):
			
				#Added unit conversion factor from ev.RE to kev.cm
				self.model_los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.eta_model[i][j]))*637100
				

	def write_pickle(self, savetag=""):
		'''This will create a pickle file of all the information that would be needed for plotting.
		This is to save an object already created. 
        
        
		'''

		# Name and locate the pickle file. 
		pickle_path = os.environ.get("PICKLE_PATH")
        
		fname = "fit_image_n_{}_SMILE_{:.2f}_{:.2f}_{:.2f}_Target_{:.2f}_{:.2f}_{:.2f}_nxm_{}_{}_{}_{}_im{}_{}.pkl".format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2], self.smile.n_pixels, self.smile.m_pixels, self.current_model, self.cost_func, self.init_method, savetag)

		# Add all the desired information to a dictionary. 
		pickle_dict = {
						"cost func":self.cost_func,
						"min cost":self.minimum_cost,
						"param list":self.param_list,
						"cost per it":self.cost_per_iteration,
						"opt time":self.opt_time,
						"model los intensity":self.model_los_intensity,
						"ppmlr los intensity":self.ppmlr_image.los_intensity,
						"params0":self.params0,
						"params best nm":self.params_best_nm,
						"filename":self.ppmlr_image.ppmlr.filename,
						"model":self.current_model,
						"r0lin":self.r0_lin if self.current_model == "cmem" else 0,
						"temp":self.temp,
						"density":self.density,
						"vx":self.vx,
						"vy":self.vy,
						"vz":self.vz,
						"bx":self.bx,
						"by":self.by,
						"bz":self.bz,
						"pdyn":self.pdyn,
						"pmag":self.pmag,
						"dipole":self.dipole,
						"init method":self.init_method,
						"smile_loc":self.smile.smile_loc,
						"target_loc":self.smile.target_loc,
						"n_pixels":self.smile.n_pixels,
						"m_pixels":self.smile.m_pixels,
						"theta_fov":self.smile.theta_fov,
						"phi_fov":self.smile.phi_fov,
						"sxi_theta":self.smile.sxi_theta,
						"sxi_phi":self.smile.sxi_phi,
						"eta_model":self.eta_model,
						"xpos":self.smile.xpos,
						"ypos":self.smile.ypos,
						"zpos":self.smile.zpos,
						"maxIx":self.ppmlr_image.ppmlr.maxIx,
						"maxdIx":self.ppmlr_image.ppmlr.maxdIx,
						"f.25":self.ppmlr_image.ppmlr.f
						}
        
		with open(os.path.join(pickle_path, self.current_model+"_optimised", fname), "wb") as f:
			pickle.dump(pickle_dict, f)

		print ("Pickled: {}".format(fname))
