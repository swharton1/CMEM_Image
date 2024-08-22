#This is a script containing functions to set the initial parameters of either the jorg or cmem models. 
#It is not intended this will be run directly. 

#METHOD 1 INITIALISING FUNCTIONS 
################################
	
def get_initial_magnetopause(bz, pdyn):
    '''This uses equation 12 in Shue et al. (1997) to estimate the initial 
    subsolar magnetopause position from Bz and Dp, which are both in the ppmlr object. '''

    if bz >= 0: 
        return (11.4 + 0.013*bz)*(pdyn**(-1.0/6.6))
    else:
        return (11.4 + 0.14*bz)*(pdyn**(-1.0/6.6))

def get_initial_alpha(bz, pdyn):
    '''This uses equation 13 in Shue et al. (1997) to estimate the initial value of alpha
    in the Shue magnetopause model. Assumes all initial alpha values will be equal to start with.'''

    return (0.58 - 0.010*bz)*(1 + 0.010*pdyn)
   

#METHOD 2 INITIALISING FUNCTIONS 
################################
	
def get_initial_mp_method2(density):
 	'''Gets mp for method 2'''
 	
 	return -0.10*density + 10.28
	
def get_initial_bs_method2(density):
    '''Gets bs for method 2'''
    
    return -0.12*density + 13.24 

def get_initial_A1_method2(density):
    '''This function estimates the initial value of the parameter A1 for the Jorgensen model. '''
    
    return 0.0000027*density - 0.0000063
    
def get_initial_A2_method2(density):
    '''This function estimates the initial value of the parameter A2. '''
    
    return 0.0000009*density - 0.0000010

def get_initial_p0_method2(density):
 	'''Gets p0 for CMEM for method 2'''
 	
 	return  0.0022*density + 0.7753

#FUNCTION TO SELECT THE INITIAL PARAMETERS
##########################################

def get_init_params(current_model, init_method, bz, pdyn, density, r0_lin=None):
	'''This will get the initial parameters of either the jorg or cmem model. Only run if params0 is not already defined. 
	
	Parameters
	----------
	current_model - 'jorg' or 'cmem'
	init_method - 1 or 2. both are described in the CMEM paper. 
	bz - IMF Bz component in nT
	pdyn - SW dynamic pressure in nPa.
	density - SW proton density in cm-3.
	r0_lin - r0_lin for CMEM model. Only needs filling in for this model.  
	
	Returns
	-------
	params0 - tuple of parameters for the current_model.
	
	'''
	
	# Sort out the initial parameters if not specified. 
	if current_model == "jorg":
        # mp, bs, A1, A2, B, alpha, beta, 
        #ay_mp, az_mp, ay_bs, az_bs.    
	    
	    
	    if init_method == 1: 
	    	# These values are rounded values taken from a 
	    	#model run in Jorgensen et al, except those 
	    	#calculate using the functions below. 
	    		
	    	# Get initial Mp using Shue et al. (1997) 
	    	#formula. Initial Bs is Mp + 3. 
	    	mp_i = get_initial_magnetopause(bz, pdyn) 
	    		
	    	# Get initial alpha values for Mp. 
	    	#Bs values are Mp + 0.2. 
	    	alpha_i = get_initial_alpha(bz, pdyn)
	    		
	    	params0 = (mp_i,mp_i+3, 0.000032, 0.000013, -0.000018, 2.5, -1.6, alpha_i, alpha_i, alpha_i+0.2, alpha_i+0.2)
	    		
	    elif init_method == 2: 
	    	mp = get_initial_mp_method2(density)
	    	bs = get_initial_bs_method2(density)
	    	A1 = get_initial_A1_method2(density)
	    	A2 = get_initial_A2_method2(density)
	    		
	    	# Get initial alpha values for Mp. Bs values are Mp + 0.2. 
	    	alpha_i = get_initial_alpha(bz, pdyn)
	    		
	    	params0 = (mp, bs, A1, A2, -0.000018, 2.5, -1.6, alpha_i, alpha_i, alpha_i+0.2, alpha_i+0.2)
                
                
	elif current_model == "cmem":
		# p0, bs, A1, A2, B, alpha, beta, 
		#p1, p2, p3, ay_bs, az_bs. 
			
		if init_method == 1:
			# A1, A2, B, alpha and beta are rounded values 
			#from Jorgensen et al. (2019).
			# p values to scale magnetopause are from 
			#inspection of an example. (1,1,3,4) 
				
			# Initial Bs is Mp + 3. 
			bs_i = r0_lin + 3
				
			# Get initial alpha values for Mp. 
			#Bs values are Mp + 0.2. 
			bs_alpha_i = get_initial_alpha(bz, pdyn) + 0.2
				
			params0 = (1, bs_i, 0.000015, 0.000013, 2, 2.5, -1.6, 1, 3, 4, bs_alpha_i, bs_alpha_i)
				
		elif init_method == 2: 
			p0 = get_initial_p0_method2(density)
			bs = get_initial_bs_method2(density)
			A1 = get_initial_A1_method2(density)
			A2 = get_initial_A2_method2(density)
				
			# Get initial alpha values for Mp. 
			#Bs values are Mp + 0.2. 
			bs_alpha_i = get_initial_alpha(bz, pdyn) + 0.2
			
			#Commented out line are the parameters used in the CMEM paper. 	
			#params0 = (p0, bs, A1, A2, 2, 2.5, -1.6, 1, 3, 4, bs_alpha_i, bs_alpha_i)
			params0 = (p0, bs, A1, A2, 2, 2.5, -1.6, 1, 3, 4, bs_alpha_i, bs_alpha_i)
				
	else:
		raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(current_model))
	
	print ('Initial parameters are: ', params0)
	
	return params0

