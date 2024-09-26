#This will take the model and smile fov object and simulate an integrated LOS image through the simulation. 

import numpy as np 
import matplotlib.pyplot as plt 
import os
from matplotlib.gridspec import GridSpec
import matplotlib.image as mpl_image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

from . import boundary_emissivity_functions as bef 
from . import get_names_and_units as gnau 
from . import set_initial_params as sip

class image_sim(): 
    '''This is the object to simulate the image.''' 
    
    def __init__(self, smile, model="jorg", init_method=2, params0=None, temp=200000, density=5, vx=400, vy=0, vz=0, bx=0, by=0, bz=5, dipole=0): 
        '''This takes in the smile object and initial parameters and adds them to self.'''
        
        #Save the image to a standard name. 
        self.plot_path = os.environ.get("PLOT_PATH") 
        
        self.smile = smile 
        
        # Extract any useful solar wind parameters
        self.temp = temp
        self.density = density
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.bx = bx
        self.by = by
        self.bz = bz
        self.pdyn = self.calc_dynamic_pressure()
        self.pmag = self.calc_magnetic_pressure()
        self.dipole = dipole 
        
        #Extract the x, y and z coords along the lines of sight and 
        #convert them to theta, phi and r.         
        self.r, self.theta, self.phi = self.convert_xyz_to_shue_coords(smile.xpos, smile.ypos, smile.zpos) 
        
        #Get the initial model parameters. 
        self.get_init_model_params(model=model, init_method=init_method, params0=params0)
        self.calc_model_image()
        
        #Get emissivity for all coodinates if doing the outreach plot. 
        #if outreach: 
        #    self.get_emissivity_all_coords()
        #    print ('Completed emissivity calculation...')
            
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
    

            
    def get_init_model_params(self, model="jorg", init_method=2, params0=None):
        '''This function will get the initial parameters for the chosen model. 
        
        Parameters
        ----------
        model - "jorg" or "cmem" 
        init_method - 1 or 2. Method to initialise parameters. 
        params0 - option to set the parameters manually. 
        
        '''
        
        self.current_model = model 
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
        

        
    #WE WILL WORK OUT THE ACTUAL IMAGE DOWN HERE. 
    #############################################
    
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
        
    def calc_model_image(self):
        '''This is the function that will actually work out the emission and LOS intensities for the given spacecraft viewing direction.''' 
        
        
        #Get the function to work out the emissivity. ('jorg' or 'cmem')  
        self.current_func = bef.get_model_func(self.current_model)
        
        #Calculate the emissivity for each LOS coordinate. 
        self.eta_model = self.get_eta_model(self.params0) 
        
        #Calculate the LOS intensity. 
        self.los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
        
        # For each pixel: 
        for i in range(self.smile.n_pixels):
            for j in range(self.smile.m_pixels):
            
                #Added unit conversion factor from ev.RE to kev.cm
                self.los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.eta_model[i][j]))*637100
    
    
    #PLOTTING FUNCTIONS 
    ###################
    
    def plot_image(self, elev=45, azim=45, cmap='hot', vmin=-8, vmax=-4, levels=100, colour_cap=0):
        '''This will plot the simulated image it has created. '''
        
        fig = plt.figure(figsize=(8,5))
        fig.subplots_adjust(left=0.05, wspace=0.2, bottom=0.20) 
        ax = fig.add_subplot(121) 
        
        #Get Image tags for model names. 
        if self.current_model == 'cmem':
            image_tag = 'CMEM'
        elif self.current_model == 'jorg':
            image_tag = 'Jorg'
        else:
            raise ValueError('Invalid model selected') 
        
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.smile.n_pixels, self.smile.n_pixels+1)-0.5
        j_array = np.linspace(0,self.smile.m_pixels, self.smile.m_pixels+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.smile.theta_fov/2.) + (self.smile.theta_fov/self.smile.n_pixels)*(I+0.5)
        phi_pixels = -(self.smile.phi_fov/2.) + (self.smile.phi_fov/self.smile.m_pixels)*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)
        
        mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0)
        ax.set_title("LOS integration through \n{} model from SMILE".format(image_tag))
        ax.set_xlabel('deg')
        ax.set_ylabel('deg')
        ax.set_aspect('equal')
        cbar = plt.colorbar(mesh, ax=ax, shrink=0.8)
        cbar.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        ax2 = fig.add_subplot(122, projection='3d')

        #Add the emissivity along the LOS.     
        los_log = np.zeros(self.eta_model.shape)+vmin
        i = np.where(self.eta_model > 0)
        los_log[i] = np.log10(self.eta_model[i])
        
        #Only plot values above a certain emissivity. 
        bright = np.where(los_log > vmin+colour_cap) 
        
        emit = ax2.scatter(self.smile.xpos[bright], self.smile.ypos[bright], self.smile.zpos[bright], c=los_log[bright], cmap="hot", s=0.05, vmin=vmin, vmax=vmax)
        cbar2 = plt.colorbar(emit, ax=ax2, shrink=0.8)
        
        cticks = np.arange(vmin, vmax)
        cbar2.set_ticks(cticks)
        cbar2.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])
        cbar2.set_label('SWCX Emissivity (eV cm'+r'$^{-3}$ s'+r'$^{-1}$)') 
        #Add FOV boundaries. 
        self.add_fov_boundaries(ax2)
        
        #Add the Earth on. 
        self.add_earth(ax2) 
        
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel('z')
        ax2.set_xlim(-10,30)
        ax2.set_title('n = {} cm'.format(self.density)+r'$^{-3}$'+'\nSMILE Coords: ({},{},{})\nAim Point: ({},{},{})'.format(self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2]))
        ax2.set_aspect('equal')
        ax2.view_init(elev,azim) 
        
        # Add a label to show the model parameters. 
        label = ""
        info = gnau.get_parameter_info(model=self.current_model)
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.params0):
                pv = pval 
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"
        
        fig.text(0.5, 0.02, label, ha='center')
                    
        
        fig.savefig(self.plot_path+"{}_image_sim_n_{}_SMILE_{}_{}_{}_Target_{}_{}_{}_nxm_{}_{}.png".format(self.current_model, self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2],self.smile.n_pixels, self.smile.m_pixels)) 
    
    
    
    
    
    #SPECIAL ADAPTATIONS FOR OUTREACH PLOT. 
    #######################################
    def get_emissivity_all_coords(self):
        '''This will calculate the emissivity using CMEM for all coordinates'''
        
        print ('Calculating emissivity for all points in space...')
        
        #Define your own set of coordinates. 
        x = np.linspace(-20,40,300)
        y = np.linspace(-40,40,400)
        z = np.linspace(-40,40,400)
        
        #Make 3D
        self.X_all, self.Y_all, self.Z_all = np.meshgrid(x,y,z)
        
        #Convert to Shue coordinates. 
        r_all, theta_all, phi_all = self.convert_xyz_to_shue_coords(self.X_all,self.Y_all,self.Z_all)
        
        #Calculate the emissivity for each LOS coordinate. 
        self.eta_all = self.current_func(r_all, theta_all, phi_all, *self.lin_coeffs, *self.params0)
         
    def outreach_plot(self, elev=45, azim=45, cmap='bone', vmin=-8, vmax=-4, levels=100, colour_cap=2, fname='test'):
        '''This will make a more funky looking outreach plot for JAC''' 
        
        #Sort out the plot configuration first. 
        plt.style.use('dark_background')
        plt.rcParams['grid.color'] = 'k'

        fig = plt.figure(figsize=(8,6))
        fig.patch.set_facecolor('k')
        gs = GridSpec(2,3, left=0.01, top=0.99, right=0.99, bottom=0.01)
        ax1 = fig.add_subplot(gs[:,:], projection='3d')
        ax1.set_facecolor('k')
        ax1.xaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.yaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.zaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.xaxis.line.set_color('k')
        ax1.yaxis.line.set_color('k')
        ax1.zaxis.line.set_color('k') 
        ax1.set_xlim(-20,30)
        ax1.set_ylim(-30,30)
        ax1.set_zlim(-30,30)

        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_zticks([])
        ax1.set_aspect('equal')

        ax2 = fig.add_subplot(gs[0,2])
        ax2.set_xticks([])
        ax2.set_yticks([])
    
        #Add the Earth on. 
        self.add_earth(ax1, color='w') 
        
        #Make pixel information for image plot. 
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.smile.n_pixels, self.smile.n_pixels+1)-0.5
        j_array = np.linspace(0,self.smile.m_pixels, self.smile.m_pixels+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.smile.theta_fov/2.) + (self.smile.theta_fov/self.smile.n_pixels)*(I+0.5)
        phi_pixels = -(self.smile.phi_fov/2.) + (self.smile.phi_fov/self.smile.m_pixels)*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        #Make image plot. 
        mesh = ax2.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0)
        ax2.set_aspect('equal')
        
        #Add the emissivity along the LOS.     
        los_log = np.zeros(self.eta_model.shape)+vmin
        i = np.where(self.eta_model > 0)
        los_log[i] = np.log10(self.eta_model[i])
        
        
        #Add the emissivity for the whole space.     
        #los_log = np.zeros(self.eta_all.shape)+vmin
        #i = np.where(self.eta_all > 0)
        #los_log[i] = np.log10(self.eta_all[i])
        
        
        #Only plot values above a certain emissivity. 
        bright = np.where(los_log > vmin+colour_cap) 
        
        emit = ax1.scatter(self.smile.xpos[bright], self.smile.ypos[bright], self.smile.zpos[bright], c=los_log[bright], cmap=cmap, s=0.001, alpha=0.2, vmin=vmin, vmax=vmax)
        
        #emit = ax1.scatter(self.X_all[bright], self.Y_all[bright], self.Z_all[bright], c=los_log[bright], cmap="hot", s=0.001, alpha=0.2, vmin=vmin, vmax=vmax)
        
        #Add FOV boundaries. 
        self.add_fov_boundaries(ax1, color='w', lw=1)
        
        #Save figure.
        if fname is None:
            fig.savefig(self.plot_path+'outreach_sim/'+fname)
        else:
            fig.savefig(fname)
    
    def outreach_plot2(self, elev=45, azim=45, cmap='bone', vmin=-8, vmax=-4, levels=100, colour_cap=2, fname=None, max_alpha=0.2):
        '''This will make a more funky looking outreach plot for JAC.
        This will plot the whole space, not just the FOV. ''' 
        
        #Sort out the plot configuration first. 
        plt.style.use('dark_background')
        plt.rcParams['grid.color'] = 'k'

        fig = plt.figure(figsize=(8,6))
        fig.patch.set_facecolor('k')
        gs = GridSpec(8,12, left=0.01, top=0.99, right=0.99, bottom=0.01)
        ax1 = fig.add_subplot(gs[:,:], projection='3d')
        ax1.set_facecolor('k')
        ax1.xaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.yaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.zaxis.set_pane_color((0.0,0.0,0.0,0.0))
        ax1.xaxis.line.set_color('k')
        ax1.yaxis.line.set_color('k')
        ax1.zaxis.line.set_color('k') 
        ax1.set_xlim(-10,20)
        ax1.set_ylim(-25,25)
        ax1.set_zlim(-25,25)
        #ax1.set_xlabel('x', color='w')

        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_zticks([])
        ax1.set_aspect('equal')

        ax2 = fig.add_subplot(gs[0:4,8:12])
        ax2.set_xticks([])
        ax2.set_yticks([])
    
        #Add the Earth on. 
        self.add_earth2(ax1) 
        
        #Make pixel information for image plot. 
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.smile.n_pixels, self.smile.n_pixels+1)-0.5
        j_array = np.linspace(0,self.smile.m_pixels, self.smile.m_pixels+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.smile.theta_fov/2.) + (self.smile.theta_fov/self.smile.n_pixels)*(I+0.5)
        phi_pixels = -(self.smile.phi_fov/2.) + (self.smile.phi_fov/self.smile.m_pixels)*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        #Make image plot. 
        mesh = ax2.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0)
        ax2.set_aspect('equal')
        
        #Add the emissivity along the LOS.     
        #los_log = np.zeros(self.eta_model.shape)+vmin
        #i = np.where(self.eta_model > 0)
        #los_log[i] = np.log10(self.eta_model[i])
        
        
        #Add the emissivity for the whole space.     
        los_log = np.zeros(self.eta_all.shape)+vmin
        i = np.where(self.eta_all > 0)
        los_log[i] = np.log10(self.eta_all[i])
        
        #Create normalised array for alpha values. 
        los_alpha = (los_log - vmin)/(vmax-vmin)
        
        #Set all alpha values to go between 0 and 0.2 
        los_alpha = los_alpha*max_alpha
        
        #Only plot values above a certain emissivity. 
        bright = np.where(los_log > vmin+colour_cap)
        #bright_inner = np.where(los_log > vmin+colour_cap+1) 
        #bright_outer = np.where((los_log > vmin+colour_cap) & (los_log <= vmin+colour_cap+1))
        #emit = ax1.scatter(self.smile.xpos[bright], self.smile.ypos[bright], self.smile.zpos[bright], c=los_log[bright], cmap=cmap, s=0.001, alpha=0.2, vmin=vmin, vmax=vmax)
        
        emit_inner = ax1.scatter(self.X_all[bright], self.Y_all[bright], self.Z_all[bright], c=los_log[bright], cmap=cmap, s=0.001, alpha=los_alpha[bright], vmin=vmin, vmax=vmax)
        #emit_outer = ax1.scatter(self.X_all[bright_outer], self.Y_all[bright_outer], self.Z_all[bright_outer], c=los_log[bright_outer], cmap=cmap, s=0.001, alpha=0.1, vmin=vmin, vmax=vmax)
        
        #Add FOV boundaries. 
        self.add_fov_boundaries(ax1, color='w', lw=1)
        
        #Now we are going to add the SMILE logo. 
        #https://towardsdatascience.com/how-to-add-an-image-to-a-matplotlib-plot-in-python-76098becaf53
        
        smile_logo_file = "/home/s/sw682/Pictures/SMILE_logo.png"
        smile_logo = mpl_image.imread(smile_logo_file) 
        
        #Create Outreach box for SMILE logo. 
        smile_imagebox = OffsetImage(smile_logo, zoom=0.3)
        smile_ab = AnnotationBbox(smile_imagebox, (0.5,0.5), frameon=False)
        
        #Create axis for the image. 
        ax3 = fig.add_subplot(gs[6:8,0:2])
        ax3.set_xlim(0,1)
        ax3.set_ylim(0,1)
        ax3.add_artist(smile_ab)
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.spines['bottom'].set_color('k')
        ax3.spines['top'].set_color('k')
        ax3.spines['left'].set_color('k')
        ax3.spines['right'].set_color('k')
        
        #Get the UOL logo.
        uol_logo_file = "/home/s/sw682/Pictures/uol_logo.png"
        uol_logo = mpl_image.imread(uol_logo_file) 
        
        #print ('Shape of UOL logo:', uol_logo.shape) 
        
        #Create Outreach box for UOL logo. 
        uol_imagebox = OffsetImage(uol_logo, zoom=0.15)
        uol_ab = AnnotationBbox(uol_imagebox, (0.5,0.5), frameon=False)
        
        #Create axis for the image. 
        ax4 = fig.add_subplot(gs[6:8,8:12])
        ax4.set_xlim(0,1)
        ax4.set_ylim(0,1)
        ax4.add_artist(uol_ab)
        ax4.set_xticks([])
        ax4.set_yticks([])
        ax4.spines['bottom'].set_color('k')
        ax4.spines['top'].set_color('k')
        ax4.spines['left'].set_color('k')
        ax4.spines['right'].set_color('k')
        
        
        #Save figure.
        if fname is None:
            fig.savefig(self.plot_path+'outreach_sim/test_v2.png')
        else:
            fig.savefig(fname)
            
            
                
    def add_fov_boundaries(self, ax2, color='k', lw=2):
        '''This will add the FOV boundaries in black/white. '''
        
        #For corner pixels only. 
        ax2.plot(self.smile.xpos[0][0], self.smile.ypos[0][0], self.smile.zpos[0][0], color, lw=lw)
        ax2.plot(self.smile.xpos[0][-1], self.smile.ypos[0][-1], self.smile.zpos[0][-1], color, lw=lw)
        ax2.plot(self.smile.xpos[-1][0], self.smile.ypos[-1][0], self.smile.zpos[-1][0], color, lw=lw)
        ax2.plot(self.smile.xpos[-1][-1], self.smile.ypos[-1][-1], self.smile.zpos[-1][-1], color, lw=lw)
        
        #Join corners together. 
        ax2.plot([self.smile.xpos[0][0][-1],self.smile.xpos[0][-1][-1]], [self.smile.ypos[0][0][-1],self.smile.ypos[0][-1][-1]], [self.smile.zpos[0][0][-1],self.smile.zpos[0][-1][-1]], color, lw=lw)
        ax2.plot([self.smile.xpos[0][-1][-1],self.smile.xpos[-1][-1][-1]], [self.smile.ypos[0][-1][-1],self.smile.ypos[-1][-1][-1]], [self.smile.zpos[0][-1][-1],self.smile.zpos[-1][-1][-1]], color, lw=lw)
        ax2.plot([self.smile.xpos[-1][-1][-1],self.smile.xpos[-1][0][-1]], [self.smile.ypos[-1][-1][-1],self.smile.ypos[-1][0][-1]], [self.smile.zpos[-1][-1][-1],self.smile.zpos[-1][0][-1]], color, lw=lw)
        ax2.plot([self.smile.xpos[-1][0][-1],self.smile.xpos[0][0][-1]], [self.smile.ypos[-1][0][-1],self.smile.ypos[0][0][-1]], [self.smile.zpos[-1][0][-1],self.smile.zpos[0][0][-1]], color, lw=lw)
        
    def add_earth(self, ax, color='k'):
        '''This will add a sphere for the Earth. '''
        
        #Create a spherical surface. 
        radius = 1
        u = np.linspace(0, 2*np.pi, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color=color, lw=0, alpha=1)
    
    def add_earth2(self, ax):
        '''This will add a sphere for the Earth. White on dayside, navy on nightside'''
        
        #Create a hemis-spherical surface that is white. 
        radius = 1
        u = np.linspace(-np.pi/2, np.pi/2, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color='w', lw=0, alpha=1, edgecolor='w', antialiased=False, shade=False)
        
        u = np.linspace(np.pi/2, 3*np.pi/2, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color='navy', lw=0, alpha=1, edgecolor='navy', antialiased=False, shade=False)
            
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
        
        
