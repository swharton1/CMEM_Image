#This class will visualise the output of the fit image class. 

import numpy as np
import matplotlib.pyplot as plt 
import pickle 
import os

from . import get_names_and_units as gnau
from . import boundary_emissivity_functions as bef

class analyse_fit():
    '''This class will contain the plotting functions that 
    the CMEM/visualise_models.py analysis class used.'''
    
    def __init__(self, filename='target/fit_image_n_20.0_SMILE_6.57_-5.94_17.33_Target_10.00_0.00_0.00_nxm_100_50_cmem_absolute_im2_.pkl', model='cmem', limb=True): 
        '''This takes in the filename for the pickle file.
        
        filename - pickle filename, inc. first folder. 
        model - 'cmem' or something else. 
        limb - Boolean to state whether the limb object was used or not. Needed to plot correctly. Basically a workaround for a minus sign error! ''' 
        
        
        self.current_model = model 
        self.filename = filename
        self.limb = limb 
        
        if self.current_model == "cmem":
            self.image_tag = "CMEM"
        else:
            self.image_tag = self.current_model.capitalize()
            
        # Name and locate the pickle file. 
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")
        
        self.full_filename = os.path.join(self.pickle_path, self.current_model+"_optimised", filename) 
        
        #Read in the pickle file. 
        self.model = self.read_pickle(self.full_filename) 
        
        #Get the names of the variables and units for plotting. 
        info = gnau.get_parameter_info(model=self.model['model'])
        self.parameter_names = [info[i][0] for i in info.keys()]
        self.parameter_units = [info[i][1] for i in info.keys()]
        
        #Get magnetopause projection information. 
        self.get_magnetopause_projection() 
        
        #Find whether SMILE is inside or outside the magnetopause. 
        self.smile_in_or_out()
        
    def get_magnetopause_projection(self):
        '''This will do all the calculations for the optimised magnetopause in order to project it into an image.'''
        
        #Get model magnetopause positions. 
        #Define ranges for theta and phi in degrees. 
        theta = np.linspace(0,120,121)
        phi = np.linspace(0,360,181)
        
        self.get_model_magnetopause(theta, phi) 
        
        #Get xyz coordinates of the magnetopause points. These are the 'magnetopause' vectors. 
        self.mpx, self.mpy, self.mpz = self.convert_shue_to_xyz_coords(self.rmp, self.theta2d, self.phi2d)
        
        
        #Get magnetopause unit vectors. 
        self.mpx_unit = self.mpx/self.rmp
        self.mpy_unit = self.mpy/self.rmp
        self.mpz_unit = self.mpz/self.rmp 
        
        #Now get the magnetopause p vectors, from the point of view of smile. 
        self.px = self.mpx - self.model['smile_loc'][0]
        self.py = self.mpy - self.model['smile_loc'][1]
        self.pz = self.mpz - self.model['smile_loc'][2] 
        
        #Get magnitude of these vectors. 
        self.pmag = np.sqrt(self.px**2 + self.py**2 + self.pz**2) 
        
        #Get the unit p vectors. 
        self.px_unit = self.px/self.pmag
        self.py_unit = self.py/self.pmag
        self.pz_unit = self.pz/self.pmag 
        
        #Now get the theta and phi of the unit p vectors. 
        self.thetap = np.arccos(self.pz_unit)
        self.phip = np.arccos(self.px_unit/(self.px_unit**2 + self.py_unit**2)**0.5) 
        
        #Now get the difference in theta and the difference 
        #in phi from the central look direction. 
        #Also convert from radians to degrees. 
        self.dthetap = np.rad2deg(self.model['sxi_theta'] - self.thetap)
        self.dphip = np.rad2deg(self.model['sxi_phi'] - self.phip) 
        
        #Also get angle between r and p to determine how close to a tangent each point is. 
        self.costangent = (self.mpx_unit*self.px_unit)+(self.mpy_unit*self.py_unit)+(self.mpz_unit*self.pz_unit) 
        self.tangent = np.rad2deg(np.arccos(self.costangent)) 
        
    
        
    def __repr__(self):
        return f"analyse model object."
    
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def trapezium_rule(self, p_spacing, eta_LOS):
        '''This will integrate a function using the trapezium rule. '''

        return (p_spacing/2)*(eta_LOS[0] + eta_LOS[-1] + 2*sum(eta_LOS[1:-1]))
        
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
        
    def convert_shue_to_xyz_coords(self, r, theta, phi):
        '''This will convert the Shue coordinates back to xyz coordinates. 
        
        Parameters
        ----------
        r, theta (rad), phi (rad)
        
        Returns
        -------
        x,y,z
        '''

        x = r*np.cos(theta)
        y = r*np.sin(theta)*np.cos(phi)
        z = r*np.sin(theta)*np.sin(phi)

        return x,y,z 
    
    def get_model_magnetopause(self, theta, phi): 
        '''This will get the magnetopause function out of the optimised model. If it's 
        'jorg', it returns a Shue model. If it's 'cmem', it returns a Lin model. 
        
        It also returns the subsolar magnetopause radius. 
        
        Parameters
        ----------
        theta - 1D array for theta
        phi - 1D array for phi 
        
        Returns
        -------
        rmp - 2D array of radial points for the magnetopause. 
        rmp_sub - subsolar magnetopause distance.
        
        ''' 
        
        theta = np.deg2rad(theta)
        phi = np.deg2rad(phi) 
        
        #Create 2D arrays. 
        self.theta2d, self.phi2d = np.meshgrid(theta, phi)
        
        if self.current_model == 'jorg':
            
            #Get the Jorgensen magnetopause for a range of theta and phi 
            self.rmp = bef.shue_func(self.theta2d, self.phi2d, self.model['params best nm'][0], self.model['params best nm'][7], self.model['params best nm'][8])
            
        elif self.current_model == 'cmem':
            
            #Get Lin coefficients. 
            self.lin_coeffs = bef.get_lin_coeffs(self.model['dipole'], self.model['pdyn'], self.model['pmag'], self.model['bz'])
            
            #Get Lin magnetopause for a range of theta and phi 
            self.rmp = bef.lin_scaled_func(self.theta2d, self.phi2d, *self.lin_coeffs, p0=self.model['params best nm'][0], p1=self.model['params best nm'][7], p2=self.model['params best nm'][8], p3=self.model['params best nm'][9]) 
            
        #Get subsolar magnetopause (theta = 0 and phi = 0) 
        sub_idx = np.where((self.theta2d == 0) & (self.phi2d == 0))
        self.rmp_sub = self.rmp[sub_idx] 
    
    def smile_in_or_out(self):
        '''This will call the magnetopause model again but only in the direction of SMILE to determine whether SMILE is inside or outside the optimised magnetopause.'''
        
        #Convert SMILE vector to Shue Coordinates. 
        smile_loc = self.model['smile_loc']
        self.smile_r = (smile_loc[0]**2 + smile_loc[1]**2 + smile_loc[2]**2)**0.5
        self.smile_theta = np.arccos(smile_loc[0]/self.smile_r)
        self.smile_phi = np.arccos(smile_loc[1]/((smile_loc[1]**2 + smile_loc[2]**2)**0.5))
        
        if self.current_model == 'jorg':
            
            #Get the Jorgensen magnetopause for a range of theta and phi 
            self.rmp_smile = bef.shue_func(self.smile_theta, self.smile_phi, self.model['params best nm'][0], self.model['params best nm'][7], self.model['params best nm'][8])
            
        elif self.current_model == 'cmem':
            
            #Get Lin magnetopause for a range of theta and phi 
            self.rmp_smile = bef.lin_scaled_func(self.smile_theta, self.smile_phi, *self.lin_coeffs, p0=self.model['params best nm'][0], p1=self.model['params best nm'][7], p2=self.model['params best nm'][8], p3=self.model['params best nm'][9]) 
        
        if self.smile_r > self.rmp_smile:
            self.inout = 'out'
        else:
            self.inout = 'in'
    #THIS IS WHERE WE CALCULATE THE OPTIMUM EMISSIVITY, AS IT'S NOT SAVED. 
    
#    def get_eta_model(self, params):
#        '''This function calculates the eta model values for each iteration. This function is intended to be run from the cost function. 
#        Parameters
#        ----------
#        params - tuple of the model parameters for the chosen model. '''
        
        #Get the function to work out the emissivity. ('jorg' or 'cmem')  
#        self.current_func = bef.get_model_func(self.current_model)
        
#        if self.current_model == "jorg": 
#            eta_model = self.current_func(self.r, self.theta, self.phi, *params)
#        elif self.current_model == "cmem":
#            eta_model = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *params)
#        else: 
#            raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
        
#        return eta_model
        

    
    def plot_change_in_parameters(self, save=False, savetag="", fname=None):
        '''This will plot how the parameters changed over the course of the
        optimisation procedure. '''

        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(hspace=0.4, wspace=0.4, top=0.85)
        fig.text(0.5, 0.9, "Parameter variation with optimisation\nn = {:.2f} cm".format(self.model['density'])+r"$^{-3}$"+", SMILE = ({},{},{}), Aim Point = ({},{},{}), nxm = {}x{}".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2], self.model['n_pixels'], self.model['m_pixels'])+" \nOptimisation Time = {:.1f}s\nModel = {}".format(self.model['opt time'], self.image_tag.capitalize()), ha="center")
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
        ax2b = ax2.twinx()
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)
        
        if self.model['model'] == "jorg":
            ax1.text(0.05, 1.05, self.parameter_names[0], c="r", transform=ax1.transAxes, fontsize=12)
            ax1.text(0.25, 1.05, self.parameter_names[1], c="b", transform=ax1.transAxes, fontsize=12)
            ax2.text(0.05, 1.05, self.parameter_names[2], c="r", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.25, 1.05, self.parameter_names[3], c="b", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.45, 1.05, self.parameter_names[4], c="g", transform=ax2.transAxes, fontsize=12)
            ax3.text(0.05, 1.05, self.parameter_names[5], c="r", transform=ax3.transAxes, fontsize=12)
            ax3.text(0.25, 1.05, self.parameter_names[6], c="b", transform=ax3.transAxes, fontsize=12)
            ax4.text(0.05, 1.05, self.parameter_names[7], c="r", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.25, 1.05, self.parameter_names[8], c="b", transform=ax4.transAxes, fontsize=12)
            ax5.text(0.05, 1.05, self.parameter_names[9], c="r", transform=ax5.transAxes, fontsize=12)
            ax5.text(0.25, 1.05, self.parameter_names[10], c="b", transform=ax5.transAxes, fontsize=12)
           
            
        elif self.model['model'] == "cmem":
            ax1.text(0.05, 1.05, self.parameter_names[0]+r'$r_0$', c="r", transform=ax1.transAxes, fontsize=12)
            ax1.text(0.25, 1.05, self.parameter_names[1], c="b", transform=ax1.transAxes, fontsize=12)
            ax2.text(0.05, 1.05, self.parameter_names[2], c="r", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.25, 1.05, self.parameter_names[3], c="b", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.85, 1.05, self.parameter_names[4], c="g", transform=ax2.transAxes, fontsize=12)
            ax3.text(0.05, 1.05, self.parameter_names[5], c="r", transform=ax3.transAxes, fontsize=12)
            ax3.text(0.25, 1.05, self.parameter_names[6], c="b", transform=ax3.transAxes, fontsize=12)
            ax4.text(0.05, 1.05, self.parameter_names[7], c="r", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.25, 1.05, self.parameter_names[8], c="b", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.45, 1.05, self.parameter_names[9], c="g", transform=ax4.transAxes, fontsize=12)
            ax5.text(0.05, 1.05, self.parameter_names[10], c="r", transform=ax5.transAxes, fontsize=12)
            ax5.text(0.25, 1.05, self.parameter_names[11], c="b", transform=ax5.transAxes, fontsize=12)
        
        # Sort cost axis and values. 
        if self.model['cost func'] == "sum squares":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost'],3)), c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']
            ax6.set_ylabel("(keV cm"+r+"$^{-3}$ s"+r"$^{-1}$ sr"+r"$^{-1})^2$", fontsize=10)
        elif self.model['cost func'] == "absolute":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost'],3)), c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']
            ax6.set_ylabel("keV cm"+r"$^{-3}$ s"+r"$^{-1}$ sr"+r"$^{-1}$", fontsize=10)
        elif self.model['cost func'] == "normalised":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost'],3)), c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']
            ax6.set_ylabel("Cost", fontsize=10)

        ax1.set_xlabel("Iterations", fontsize=10)
        ax2.set_xlabel("Iterations", fontsize=10)
        ax3.set_xlabel("Iterations", fontsize=10)
        ax4.set_xlabel("Iterations", fontsize=10)
        ax5.set_xlabel("Iterations", fontsize=10)
        ax6.set_xlabel("Iterations", fontsize=10)

        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax2.get_xticklabels() + ax2.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax3.get_xticklabels() + ax3.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax4.get_xticklabels() + ax4.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax5.get_xticklabels() + ax5.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax6.get_xticklabels() + ax6.get_yticklabels()): 
            label.set_fontsize(8)

        iteration = np.arange(len(self.model['param list']))
        param_list_t = np.array(self.model['param list']).transpose()
       
        
        if self.model['model'] == "jorg":
            ax1.plot(iteration, param_list_t[0], "r")
            ax1.plot(iteration, param_list_t[1], "b")
            ax2.plot(iteration, param_list_t[2]*100000, "r")
            ax2.plot(iteration, param_list_t[3]*100000, "b")
            ax2.plot(iteration, param_list_t[4]*100000, "g")
            ax3.plot(iteration, param_list_t[5], "r")
            ax3.plot(iteration, param_list_t[6], "b")
            ax4.plot(iteration, param_list_t[7], "r")
            ax4.plot(iteration, param_list_t[8], "b")
            ax5.plot(iteration, param_list_t[9], "r")
            ax5.plot(iteration, param_list_t[10], "b")
            ax6.plot(iteration, cpi, "k", label="Cost")
        elif self.model['model'] == "cmem":
            ax1.plot(iteration, param_list_t[0]*self.model['r0lin'], "r")
            ax1.plot(iteration, param_list_t[1], "b")
            ax2.plot(iteration, param_list_t[2]*100000, "r")
            ax2.plot(iteration, param_list_t[3]*100000, "b")
            ax2b.plot(iteration, param_list_t[4], "g")
            ax3.plot(iteration, param_list_t[5], "r")
            ax3.plot(iteration, param_list_t[6], "b")
            ax4.plot(iteration, param_list_t[7], "r")
            ax4.plot(iteration, param_list_t[8], "b")
            ax4.plot(iteration, param_list_t[9], "g")
            ax5.plot(iteration, param_list_t[10], "r")
            ax5.plot(iteration, param_list_t[11], "b")
            ax6.plot(iteration, cpi, "k", label="Cost")

        # If boundaries were applied to parameters, plot them on. NOT ADAPTED FOR JORGENSEN LIN. 
#        if self.model['param bounds'] is not None: 
#            pbounds = self.model['param bounds'] 
#            ax1.plot([iteration[0], iteration[-1]], [pbounds[0][0], pbounds[0][0]], 'r--')
#            ax1.plot([iteration[0], iteration[-1]], [pbounds[0][1], pbounds[0][1]], 'r--')
#            ax1.plot([iteration[0], iteration[-1]], [pbounds[1][0], pbounds[1][0]], 'b--')
#            ax1.plot([iteration[0], iteration[-1]], [pbounds[1][1], pbounds[1][1]], 'b--')
#            ax4.plot([iteration[0], iteration[-1]], [pbounds[7][0], pbounds[7][0]], 'r--')
#            ax4.plot([iteration[0], iteration[-1]], [pbounds[7][1], pbounds[7][1]], 'r--')
#            ax4.plot([iteration[0], iteration[-1]], [pbounds[8][0], pbounds[8][0]], 'b--')
#            ax4.plot([iteration[0], iteration[-1]], [pbounds[8][1], pbounds[8][1]], 'b--')
#            ax5.plot([iteration[0], iteration[-1]], [pbounds[9][0], pbounds[9][0]], 'r--')
#            ax5.plot([iteration[0], iteration[-1]], [pbounds[9][1], pbounds[9][1]], 'r--')
#            ax5.plot([iteration[0], iteration[-1]], [pbounds[10][0], pbounds[10][0]], 'b--')
#            ax5.plot([iteration[0], iteration[-1]], [pbounds[10][1], pbounds[10][1]], 'b--')
           

        # Put unit labels on axes where necessary. 
        ax1.set_ylabel(r"$R_E$", fontsize=8)
        ax2.set_ylabel(self.parameter_units[2], fontsize=8)
       
        if save: 
            if fname is None:
                fig.savefig(self.plot_path+"fitted_images/{}_parameter_changes{}.png".format(self.filename.split("_.pkl")[0], savetag))
            else:
                fig.savefig(self.plot_path+fname)

        self.fig_param = fig 


    def plot_images(self, cmap='hot', vmin=-8, vmax=-4, levels=100, los_max=12, save=False, savetag="", add_mp_projection=True, fname=None, ellipse=None, elev=45, azim=45, add_fov_projection=False, colour_cap=0):
        '''This will plot the final model image alongside the PPMLR image.
        
        Parameters
        ----------
        cmap - def = 'hot'
        vmin - redundant
        vmax - redundant
        levels - redundant
        los_max - max LOS intensity to plot to.
        save - boolean to save. def = False
        savetag - option to add extra text to filename. 
        add_mp_projection - boolean to add a projection of the magnetopause on. 
        fname - filename to save the plot to. 
        ellipse - ellipse object to add elliptical orbit to plot.
        elev - elevation of 3d plot
        azim - azimuth of 3d plot 
        
        ''' 
        
        if add_fov_projection:
            fig = plt.figure(figsize=(10,5))
            fig.subplots_adjust(left=0.05, wspace=0.2, bottom=0.20) 
            ax1 = fig.add_subplot(131) 
            ax2 = fig.add_subplot(132)
            ax3 = fig.add_subplot(133, projection='3d')
        else:
            fig = plt.figure(figsize=(8,5))
            fig.subplots_adjust(left=0.10, wspace=0.5, bottom=0.20) 
            ax1 = fig.add_subplot(121) 
            ax2 = fig.add_subplot(122)
            
            
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.model['n_pixels'], self.model['n_pixels']+1)-0.5
        j_array = np.linspace(0,self.model['m_pixels'], self.model['m_pixels']+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.model['theta_fov']/2.) + (self.model['theta_fov']/self.model['n_pixels'])*(I+0.5)
        phi_pixels = -(self.model['phi_fov']/2.) + (self.model['phi_fov']/self.model['m_pixels'])*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        #if not hasattr(self.smile, 'radial_constraint'):
        if not self.limb:
            print ('Reverse angular arrays for plotting if not updated SMILE FOV limb object.') 
            theta_pixels = -np.rad2deg(theta_pixels)
            phi_pixels = -np.rad2deg(phi_pixels)
        
        # Get contour levels. 
        #levels = np.linspace(vmin, vmax, levels+1)
        
        mesh1 = ax1.pcolormesh(phi_pixels, theta_pixels, self.model['ppmlr los intensity'], cmap=cmap, vmin=0, vmax=los_max)
        if add_fov_projection:
            ax1.set_title("PPMLR Image from SMILE\nn = {} cm".format(self.model['density'])+r"$^{-3}$", fontsize=10)
        else:
            ax1.set_title("PPMLR Image from SMILE\nSMILE = ({:.2f},{:.2f},{:.2f}), Aim Point = ({},{},{})".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2]), fontsize=10)
        ax1.set_xlabel('deg')
        if not add_fov_projection: ax1.set_ylabel('deg')
        ax1.set_aspect('equal')
        #if not add_fov_projection:
        cbar1 = plt.colorbar(mesh1, ax=ax1, shrink=0.8)
        if not add_fov_projection:
            cbar1.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        #Now add the model image. 
        mesh2 = ax2.pcolormesh(phi_pixels, theta_pixels, self.model['model los intensity'], cmap=cmap, vmin=0, vmax=los_max)
        if self.inout == 'out':
            ax2.set_title("{} Image from SMILE\nSMILE is outside MP".format(self.image_tag), fontsize=10)
        else:
            ax2.set_title("{} Image from SMILE\nSMILE is inside MP".format(self.image_tag), fontsize=10)
        
        ax2.set_xlabel('deg')
        if not add_fov_projection: ax2.set_ylabel('deg')
        ax2.set_aspect('equal')
        cbar2 = plt.colorbar(mesh2, ax=ax2, shrink=0.8)
        cbar2.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        #Add a projection of the magnetopause here. 
        #ax2.scatter(self.dphip, self.dthetap, c='w', s=5)
        #if add_mp_projection:
        
            #Add lines going for each constant value of phi. 
        #    for p in range(len(self.dphip)):
        #        front = np.where(self.tangent[p] > 90)
        #        back = np.where(self.tangent[p] <= 90)
        #        ax2.plot(self.dphip[p][front], self.dthetap[p][front], c='w', lw=0.5)
        #        ax2.plot(self.dphip[p][back], self.dthetap[p][back], c='gray', lw=0.5)
        
            #Transpose arrays to plot cylindrical lines. 
        
        #    for t in range(len(self.dphip[0])):
        #        front = np.where(self.tangent[:,t] > 90)
        #        back = np.where(self.tangent[:,t] <= 90)
        #        ax2.plot(self.dphip[:,t][front], self.dthetap[:,t][front], c='w', lw=0.5)
        #        ax2.plot(self.dphip[:,t][back], self.dthetap[:,t][back], c='gray', lw=0.5)
        
            #Make sure only FOV is shown, even if MP projection goes outside it. 
        #    ax2.set_xlim(phi_pixels.min(), phi_pixels.max())
        #    ax2.set_ylim(theta_pixels.min(), theta_pixels.max())
        
        #Add a label for the subsolar magnetopause radial position. 
        
        
        #Add a label for the PPMLR subsolar magnetopause position. 
        if add_fov_projection:
            fig.text(0.30, 0.2, r"$r_0$"+" = {:.2f}".format(self.model['maxIx'])+r"$R_E$", ha='center')
            fig.text(0.60, 0.2, r"$r_0$"+" = {:.2f}".format(self.rmp_sub[0])+r"$R_E$", ha='center')
        else:
            fig.text(0.42, 0.2, r"$r_0$"+" = {:.2f}".format(self.model['maxIx'])+r"$R_E$", ha='center')
            fig.text(0.90, 0.2, r"$r_0$"+" = {:.2f}".format(self.rmp_sub[0])+r"$R_E$", ha='center')
        
        if add_fov_projection:
            #Add the 3rd axis to show the projection angles better. 
            #Add the emissivity along the LOS.
            eta_model = self.model['eta_model']
            
            los_log = np.zeros(eta_model.shape)+vmin
            i = np.where(eta_model > 0)
            los_log[i] = np.log10(eta_model[i])
        
            #Only plot values above a certain emissivity. 
            bright = np.where(los_log > vmin+colour_cap) 
        
            emit = ax3.scatter(self.model['xpos'][bright], self.model['ypos'][bright], self.model['zpos'][bright], c=los_log[bright], cmap="hot", s=0.001, alpha=0.2, vmin=vmin, vmax=vmax)
            #cbar3 = plt.colorbar(emit, ax=ax3, shrink=0.8)
        
            #cticks = np.arange(vmin, vmax)
            #cbar3.set_ticks(cticks)
            #cbar3.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])
            #cbar3.set_label('SWCX Emissivity (eV cm'+r'$^{-3}$ s'+r'$^{-1}$)') 
        
            #Add FOV boundaries. 
            self.add_fov_boundaries(ax3)
        
            #Add the Earth on. 
            self.add_earth(ax3) 
        
            #Add orbit ellipse if it is provided.
            if ellipse is not None: 
                ax3.plot(ellipse.coords_gsm.x/ellipse.RE, ellipse.coords_gsm.y/ellipse.RE, ellipse.coords_gsm.z/ellipse.RE, 'k')
        
            #Add on the Lin magnetopause from the optimised model. 
            if add_mp_projection: 
                
                #Add lines going for each constant value of phi. 
                for p in range(len(self.mpx)):
                    #out = np.where((self.dthetap[p] > theta_pixels.max()) | (self.dthetap[p] < theta_pixels.min()) | (self.dphip[p] > phi_pixels.max()) | (self.dphip[p] < phi_pixels.min()))
                    
                    ax3.plot(self.mpx[p], self.mpy[p], self.mpz[p], c='k', lw=0.5, zorder=0, alpha=0.2) 
                for p in range(len(self.mpx[0])):
                    ax3.plot(self.mpx[:,p], self.mpy[:,p], self.mpz[:,p], c='k', lw=0.5, zorder=0, alpha=0.5) 
            ax3.set_xlabel('x')
            ax3.set_ylabel('y')
            #ax3.set_zlabel('z')
            ax3.set_xlim(-10,30)
            ax3.set_ylim(-30,30)
            ax3.set_zlim(-30,30)
            ax3.set_title("SMILE = ({:.2f},{:.2f},{:.2f})\nAim Point = ({:.2f},{:.2f},{:.2f})".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2]), fontsize=10)
            ax3.set_aspect('equal')
            ax3.view_init(elev,azim) 
        
        #Add a line from the Earth to the spacecraft.
        #ax3.plot([0,self.model['smile_loc'][0]], [0, self.model['smile_loc'][1]], [0, self.model['smile_loc'][2]], 'cyan')
        
        # Add a label to show the model parameters. 
        label = ""
        for p,pval in enumerate(self.model['params best nm']):
                pv = pval 
                label += "{}={} {}, ".format(self.parameter_names[p], self.sig_figs(pv,3), self.parameter_units[p])
                if len(self.parameter_names)//2 == p+1:
                    label += "\n"
        
        fig.text(0.5, 0.02, label, ha='center')
        
        if save: 
        
            if fname is None:
                #Add tag if mp is added. 
                mp_tag = 'mp' if add_mp_projection else ''
                fov_tag = 'fov' if add_fov_projection else ''
                fig.savefig(self.plot_path+"fitted_images/{}_images_{}_{}{}.png".format(self.filename.split("_.pkl")[0], mp_tag, fov_tag, savetag))
            else:
                #You need the full path here. 
                fig.savefig(self.plot_path+fname)
        self.fig_param = fig 
    
    def plot_images_sequence(self, cmap='hot', los_max=60):
        '''This will plot the CMEM image alongside the PPMLR image, but one figure after another. The images can then be stitched together into a gif.'''
        
        import matplotlib.animation as animation 
        
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.model['n_pixels'], self.model['n_pixels']+1)-0.5
        j_array = np.linspace(0,self.model['m_pixels'], self.model['m_pixels']+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.model['theta_fov']/2.) + (self.model['theta_fov']/self.model['n_pixels'])*(I+0.5)
        phi_pixels = -(self.model['phi_fov']/2.) + (self.model['phi_fov']/self.model['m_pixels'])*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        #Make initial figure. 
        plt.close("all")
        fig = plt.figure(figsize=(8,5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        fig.subplots_adjust(left=0.10, wspace=0.5, bottom=0.20) 
    
        #Add the PPMLR image, which is always the same. 
        mesh1 = ax1.pcolormesh(phi_pixels, theta_pixels, self.model['ppmlr los intensity'], cmap=cmap, vmin=0, vmax=los_max)
        
        ax1.set_title("PPMLR Image from SMILE\nSMILE = ({:.2f},{:.2f},{:.2f}), Aim Point = ({},{},{})".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2]), fontsize=10)
        ax1.set_xlabel('deg')
        ax1.set_ylabel('deg')
        ax1.set_aspect('equal')
        
        cbar1 = plt.colorbar(mesh1, ax=ax1, shrink=0.8)
        cbar1.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        #Add on subsolar magnetopause position for PPMLR. 
        fig.text(0.42, 0.2, r"$r_0$"+" = {:.2f}".format(self.model['maxIx'])+r"$ R_E$", ha='center')
            
        #Now we need to add on the initial model image. 
        
        #LOS coords in cartesian. 
        x = self.model['xpos']
        y = self.model['ypos']
        z = self.model['zpos'] 
        
        #LOS coords in Shue. 
        self.r, self.theta, self.phi = self.convert_xyz_to_shue_coords(x, y, z) 
        
        #Get the Lin coefficients and initial emissivity. 
        self.lin_coeffs = bef.get_lin_coeffs(self.model['dipole'], self.model['pdyn'], self.model['pmag'], self.model['bz'])
        self.current_func = bef.get_model_func('cmem')
        eta_model_init = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *self.model['params0'])
        
        #Calculate the LOS intensity. 
        los_intensity_init = np.zeros((self.model['n_pixels'], self.model['m_pixels']))
        
        # For each pixel: 
        for i in range(self.model['n_pixels']):
            for j in range(self.model['m_pixels']):
            
                #Added unit conversion factor from ev.RE to kev.cm
                los_intensity_init[i][j] = ((1/(4*np.pi))*self.trapezium_rule(0.5, eta_model_init[i][j]))*637100
                
        #Now add the model image. 
        mesh2 = ax2.pcolormesh(phi_pixels, theta_pixels, los_intensity_init, cmap=cmap, vmin=0, vmax=los_max)
        if self.inout == 'out':
            ax2.set_title("{} Image from SMILE\nSMILE is outside MP".format(self.image_tag), fontsize=10)
        else:
            ax2.set_title("{} Image from SMILE\nSMILE is inside MP".format(self.image_tag), fontsize=10)
        
        ax2.set_xlabel('deg')
        ax2.set_ylabel('deg')
        ax2.set_aspect('equal')
        cbar2 = plt.colorbar(mesh2, ax=ax2, shrink=0.8)
        cbar2.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        #Get initial subsolar magnetopause position for CMEM. 
        rmp_cmem = bef.lin_scaled_func(0, 0, *self.lin_coeffs, p0=self.model['params0'][0], p1=self.model['params0'][7], p2=self.model['params0'][8], p3=self.model['params0'][9]) 
        rmp_text = fig.text(0.90, 0.2, r"$r_0$"+" = {:.2f}".format(rmp_cmem)+r" $R_E$", ha='center')
        
        
        # Add a label to show the model parameters. 
        label = ""
        info = gnau.get_parameter_info(model='cmem')
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.model['params0']):
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pval,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"
        
        labeltext = fig.text(0.5, 0.02, label, ha='center')
        
        self.mesh2 = mesh2
        
        def update(frame):
            '''This will update the CMEM image for each frame using the next set of parameters.'''
            
            #Recalculate emissivity. 
            new_params = tuple(self.model['param list'][frame])
            print ('Frame = ', frame)
            new_eta = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *new_params)
            
            #Recalculate LOS intensity. 
            new_los_intensity = np.zeros((self.model['n_pixels'], self.model['m_pixels']))
            
            # For each pixel: 
            for i in range(self.model['n_pixels']):
                for j in range(self.model['m_pixels']):
            
                    #Added unit conversion factor from ev.RE to kev.cm
                    new_los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(0.5, new_eta[i][j]))*637100
            
            #Reset data in mesh2 object. 
            mesh2.set(array = new_los_intensity)
            
            #Update subsolar magnetopause distance from CMEM. 
            rmp_cmem = bef.lin_scaled_func(0, 0, *self.lin_coeffs, p0=new_params[0], p1=new_params[7], p2=new_params[8], p3=new_params[9])
            rmp_label = r"$r_0$"+" = {:.2f}".format(rmp_cmem)+r" $R_E$"
            rmp_text.set_text(rmp_label)
            
            #Update the label too. 
            label = "" 
            for p,pval in enumerate(new_params):
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pval,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"
            labeltext.set_text(label)
            
            return mesh2, rmp_text, labeltext 
        
        #Now make the animation. 
        ani = animation.FuncAnimation(fig=fig, func=update,  frames=len(self.model['param list']), interval=20) 
        ani.save(self.plot_path+'fitting_animation.gif')
        #plt.show()    
            
        
        
        
        
        
        
        
        
                    
    def add_fov_boundaries(self, ax2):
        '''This will add the FOV boundaries in black. '''
        
        #For corner pixels only. 
        ax2.plot(self.model['xpos'][0][0], self.model['ypos'][0][0], self.model['zpos'][0][0], 'k', lw=0.5, zorder=3)
        ax2.plot(self.model['xpos'][0][-1], self.model['ypos'][0][-1], self.model['zpos'][0][-1], 'k', lw=0.5, zorder=3)
        ax2.plot(self.model['xpos'][-1][0], self.model['ypos'][-1][0], self.model['zpos'][-1][0], 'k', lw=0.5, zorder=3)
        ax2.plot(self.model['xpos'][-1][-1], self.model['ypos'][-1][-1], self.model['zpos'][-1][-1], 'k', lw=0.5, zorder=3)
        
        #Join corners together. 
        ax2.plot([self.model['xpos'][0][0][-1],self.model['xpos'][0][-1][-1]], [self.model['ypos'][0][0][-1],self.model['ypos'][0][-1][-1]], [self.model['zpos'][0][0][-1],self.model['zpos'][0][-1][-1]], 'k', lw=0.5, zorder=3)
        ax2.plot([self.model['xpos'][0][-1][-1],self.model['xpos'][-1][-1][-1]], [self.model['ypos'][0][-1][-1],self.model['ypos'][-1][-1][-1]], [self.model['zpos'][0][-1][-1],self.model['zpos'][-1][-1][-1]], 'k', lw=0.5, zorder=3)
        ax2.plot([self.model['xpos'][-1][-1][-1],self.model['xpos'][-1][0][-1]], [self.model['ypos'][-1][-1][-1],self.model['ypos'][-1][0][-1]], [self.model['zpos'][-1][-1][-1],self.model['zpos'][-1][0][-1]], 'k', lw=0.5, zorder=3)
        ax2.plot([self.model['xpos'][-1][0][-1],self.model['xpos'][0][0][-1]], [self.model['ypos'][-1][0][-1],self.model['ypos'][0][0][-1]], [self.model['zpos'][-1][0][-1],self.model['zpos'][0][0][-1]], 'k', lw=0.5, zorder=3)
        
    def add_earth(self, ax):
        '''This will add a sphere for the Earth. '''
        
        #Create a spherical surface. 
        radius = 1
        u = np.linspace(0, 2*np.pi, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color='k', lw=0, alpha=1)
            
        
    def sig_figs(self, x: float, precision: int):
        """
        Rounds a number to number of significant figures
        Parameters:
        - x - the number to be rounded
        - precision (integer) - the number of significant figures
        Returns:
        - float
        """

        x = float(x)
        precision = int(precision)

        return np.round(x, -int(np.floor(np.log10(abs(x)))) + (precision - 1))
