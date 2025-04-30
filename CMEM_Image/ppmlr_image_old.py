#This will resample the PPMLR data into the SMILE FOV and work out an image. 

import numpy as np
import matplotlib.pyplot as plt
import os

from SXI_Core import calc_pressures
from SXI_Core import get_earth 
from SXI_Core import add_fov_boundaries

class ppmlr_image():
    '''This class takes in the ppmlr simulation object and the smile fov object and calculates an image through the simulation.'''
    
    def __init__(self, ppmlr, smile):
        '''This takes in the ppmlr and smile objects.'''
        
        self.ppmlr = ppmlr
        self.smile = smile 
        
        # Extract any useful solar wind parameters
        self.temp = ppmlr.temp
        self.density = ppmlr.density
        self.vx = ppmlr.vx
        self.vy = ppmlr.vy
        self.vz = ppmlr.vz
        self.bx = ppmlr.bx
        self.by = ppmlr.by
        self.bz = ppmlr.bz
        self.pdyn = calc_pressures.calc_dynamic_pressure(self.vx, self.vy, self.vz, self.density)
        self.pmag = calc_pressures.calc_magnetic_pressure(self.bx, self.by, self.bz)
        
        #Calculate the emissivity along the LOS. 
        self.get_weighted_eta_in_fov()
        
        #Calculate the LOS intensity that forms the image. 
        self.calc_model_image()
        
        self.plot_path = os.environ.get("PLOT_PATH") 
    

    #THIS CODE WILL WORK OUT THE EMISSIVITY ALONG THE LOS FROM THE PPMLR MODEL.
    ###########################################################################
    
    def get_weighted_eta_in_fov(self):
        '''This will get the weighted eta values for all points in the FOV. 
        '''
        
        print ("Calculating emissivity from weighted average method...") 
        #Loop through each pixel and each point along the LOS of each pixel. 
        self.peta = np.zeros((self.smile.xpos.shape))
        for i in range(self.smile.n_pixels):
            print ("Pixel i = ", i)
            for j in range(self.smile.m_pixels):
                for p in range(self.smile.xpos[0][0].size):
                    px = self.smile.xpos[i][j][p]
                    py = self.smile.ypos[i][j][p]
                    pz = self.smile.zpos[i][j][p]
                    self.peta[i][j][p] = self.get_weighted_eta_single_value(px, py, pz)
        
        
    def get_weighted_eta_single_value(self, px, py, pz):
        '''This will get the nearest x, y and z values for a given point. 
        
        Parameters
        ----------
        px - x coordinate of point in Smile Fov
        py - y coordinate of point in Smile Fov
        pz - z coordinate of point in Smile Fov 
        
        Returns
        -------
        peta - weighted emissivity value for the coords (px, py, pz) 
        
        '''
        
        #If px, py or pz are outside the boundaries of the simulation, eta is zero. 
        if (px < self.ppmlr.x[0]) or (px > self.ppmlr.x[-1]) or (py < self.ppmlr.y[0]) or (py > self.ppmlr.y[-1]) or (pz < self.ppmlr.z[0]) or (pz > self.ppmlr.z[-1]):
            peta = 0 
            return peta
            
        else: 
            #It must be inside the cube.
            
            #Get the indices of x0, y0 and z0.  
            ix = self.get_x0_x1(px) 
            iy = self.get_y0_y1(py)
            iz = self.get_z0_z1(pz) 
            
            #Get the x, y, z and eta values for the vertices. 
            self.vert_x = self.ppmlr.x_3d[iz:iz+2, iy:iy+2, ix:ix+2]
            self.vert_y = self.ppmlr.y_3d[iz:iz+2, iy:iy+2, ix:ix+2]
            self.vert_z = self.ppmlr.z_3d[iz:iz+2, iy:iy+2, ix:ix+2]
            self.etav = self.ppmlr.eta_3d[iz:iz+2, iy:iy+2, ix:ix+2]
            
            #Get the radial distances to each vertex. 
            r_vertices = self.get_r_to_vertex(px, py, pz, self.vert_x, self.vert_y, self.vert_z)
            
            #Need to account for possibility (px,py,pz) is an exact 
            #point in the ppmlr grid. 
            if 0 in r_vertices:
                peta = self.etav[r_vertices==0][0]
            else:
                #Get weights for each vertex.
                weights = 1/r_vertices 
                
                #Calculate the weighted eta value for the coordinates (px, py, pz) 
                peta = (weights*self.etav).sum()/weights.sum()
            
            return peta 
    
    def get_r_to_vertex(self, px, py, pz, vx, vy, vz):
        '''This calculates the radial distance to a vertex. '''
        
        #r = np.sqrt((vertex[0]-px)**2 + (vertex[1]-py)**2 + (vertex[2]-pz)**2)
        r = np.sqrt((vx-px)**2 + (vy-py)**2 + (vz-pz)**2)
        return r 
        
        
    def get_x0_x1(self, px):
        '''This will get the index of x0'''
        
        diff = self.ppmlr.x - px
        ix = diff[diff <= 0].argmax()    
        return ix
        
        
    def get_y0_y1(self, py):
        '''This will get the index of y0'''
            
        diff = self.ppmlr.y - py
        iy = diff[diff <= 0].argmax()    
        return iy
        
        
    def get_z0_z1(self, pz):
        '''This will get the index of z0'''
        
        diff = self.ppmlr.z - pz
        iz = diff[diff <= 0].argmax()    
        return iz
        
    #THIS CODE WILL WORK OUT THE INTENSITY IMAGE
    ######################################################
    
    def trapezium_rule(self, p_spacing, eta_LOS):
        '''This will integrate a function using the trapezium rule. '''

        return (p_spacing/2)*(eta_LOS[0] + eta_LOS[-1] + 2*sum(eta_LOS[1:-1]))
    
    def calc_model_image(self):
        '''This is the function that will actually work out the emission and LOS intensities for the given spacecraft viewing direction.''' 
        
        print ("Calculating LOS intensity...") 
        #Calculate the LOS intensity. 
        self.los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
        
        # For each pixel: 
        for i in range(self.smile.n_pixels):
            for j in range(self.smile.m_pixels):
            
                #Added unit conversion factor from ev.RE to kev.cm
                self.los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.peta[i][j]))*637100
    
    def calc_model_image_sections(self, r_bounds = [0,8,16,24,32,40]):
        '''This will create a set of images of the emission from different distances, like Andy's code.'''
        
        print ("Calculating LOS intensity for a range of distances...")
        
        self.r_bounds = r_bounds
        
        self.los_intensity_sections = np.zeros((len(self.r_bounds)-1, self.smile.n_pixels, self.smile.m_pixels)) 
        
        #Loop through each section. 
        for r in range(len(self.r_bounds)-1):
            
            # For each pixel: 
            for i in range(self.smile.n_pixels):
                for j in range(self.smile.m_pixels):
                    
                    #Filter to just use emissivities in this range. 
                    idx1 = int(self.r_bounds[r]/self.smile.p_spacing)
                    idx2 = int(self.r_bounds[r+1]/self.smile.p_spacing)
                    
                    #Added unit conversion factor from ev.RE to kev.cm
                    self.los_intensity_sections[r][i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.peta[i][j][idx1:idx2]))*637100
    
    
    
    #PLOTTING FUNCTIONS 
    ###################
    
    def plot_image(self, elev=45, azim=45, cmap='hot', vmin=-8, vmax=-4, levels=100, colour_cap=0, los_max=6, image_name=None, ellipse=None):
        '''This will plot the simulated image it has created. 
        
        Parameters
        ----------
        elev - viewing elevation in degrees for 3d viewing model
        azim - viewing azimuth in degrees for 3d viewing model 
        cmap - colourmap
        vmin - min emissivity on colour bar (logged)
        vmax - max emissivity on colour bar (logged)
        levels - number of levels on contour maps
        colour_cap - order of magnitude above vmin to start plotting emissivity 
        los_max - max los intensity on colourbar
        image_name - will default to standard name if not specified. must be full path
        ellipse - ellipse object if you wish to add orbit ellipse. 
        
        '''
        
        fig = plt.figure(figsize=(8,5))
        fig.subplots_adjust(left=0.05, wspace=0.2, bottom=0.20) 
        ax = fig.add_subplot(121) 
        
        
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
        
        mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0, vmax=los_max)
        ax.set_title("LOS integration through PPMLR\n simulation from SMILE")
        ax.set_xlabel('deg')
        ax.set_ylabel('deg')
        ax.set_aspect('equal')
        cbar = plt.colorbar(mesh, ax=ax, shrink=0.8)
        cbar.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        ax2 = fig.add_subplot(122, projection='3d')

        #Add the emissivity along the LOS.     
        los_log = np.zeros(self.peta.shape)+vmin
        i = np.where(self.peta > 0)
        los_log[i] = np.log10(self.peta[i])
        
        #Only plot values above a certain emissivity. 
        bright = np.where(los_log > vmin+colour_cap) 
        
        emit = ax2.scatter(self.smile.xpos[bright], self.smile.ypos[bright], self.smile.zpos[bright], c=los_log[bright], cmap="hot", s=0.05, vmin=vmin, vmax=vmax)
        cbar2 = plt.colorbar(emit, ax=ax2, shrink=0.8)
        
        cticks = np.arange(vmin, vmax)
        cbar2.set_ticks(cticks)
        cbar2.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])
        cbar2.set_label('SWCX Emissivity (eV cm'+r'$^{-3}$ s'+r'$^{-1}$)') 
        #Add FOV boundaries. 
        add_fov_boundaries.add_fov_boundaries(ax2, self.smile)
        
        #Add the Earth on. 
        get_earth.make_earth_3d(ax2) 
        
        #Add orbit ellipse if it is provided.
        if ellipse is not None: 
            ax2.plot(ellipse.x3/ellipse.RE, ellipse.y3/ellipse.RE, ellipse.z3/ellipse.RE, 'k')
        
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel('z')
        ax2.set_xlim(-10,30)
        ax2.set_ylim(-30,30)
        ax2.set_zlim(-30,30)
        ax2.set_title('n = {} cm'.format(self.density)+r'$^{-3}$'+'\nSMILE Coords: ({:.2f},{:.2f},{:.2f})\nAim Point: ({},{},{})'.format(self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2]))
        ax2.set_aspect('equal')
        ax2.view_init(elev,azim) 
        
        # Add a label to show the model parameters. 
        #label = ""
        #info = gnau.get_parameter_info(model=self.current_model)
        #parameter_names = [info[i][0] for i in info.keys()]
        #parameter_units = [info[i][1] for i in info.keys()]
        #for p,pval in enumerate(self.params0):
        #        pv = pval 
        #        label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
        #        if len(parameter_names)//2 == p+1:
        #            label += "\n"
        
        #fig.text(0.5, 0.02, label, ha='center')
                    
        #Save the image to a standard name. 
        if image_name is None: 
            
        
            fig.savefig(self.plot_path+"PPMLR_image_sim_n_{}_SMILE_{}_{}_{}_Target_{}_{}_{}_nxm_{}_{}.png".format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2],self.smile.n_pixels, self.smile.m_pixels)) 
        
        else:
            fig.savefig(image_name) 
    
    def plot_image_sections(self, elev=45, azim=45, cmap='hot', vmin=-8, vmax=-4, levels=100, colour_cap=0, los_max=None, image_name=None, ellipse=None, save=False, savetag=''):
        '''This will plot the series of images created from emission at different distances, plus the total 
        at the bottom. Akin to Andy Read's code. 
        
        Parameters
        ----------
        elev - viewing elevation in degrees for 3d viewing model
        azim - viewing azimuth in degrees for 3d viewing model 
        cmap - colourmap
        vmin - min emissivity on colour bar (logged)
        vmax - max emissivity on colour bar (logged)
        levels - number of levels on contour maps
        colour_cap - order of magnitude above vmin to start plotting emissivity 
        los_max - max los intensity on colourbar
        image_name - will default to standard name if not specified. must be full path
        ellipse - ellipse object if you wish to add orbit ellipse. 
        
        '''
        
        plt.close("all") 
        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(left=0.10, wspace=0.3, bottom=0.10) 
        
        n_axes = len(self.r_bounds) 
        rows = n_axes/2
        
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
    
        print (n_axes)
        for n in range(n_axes):
            #print (n, self.r_bounds[n], n%3, n//3)
            ax = fig.add_subplot(int(n_axes/rows), int(rows), n+1)
            if n == n_axes-1:
                #Do the total image. 
                mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0, vmax=los_max*3)
                ax.set_title("{} < r < {}".format(0, self.smile.p_max))
                
            else: 
                #Do the sections images. 
                mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity_sections[n], cmap=cmap, vmin=0, vmax=los_max)
                ax.set_title("{} < r < {}".format(self.r_bounds[n],self.r_bounds[n+1]))
            
            #Add customisation to the plot. 
            if n//3 == 1: ax.set_xlabel('deg')
            if n%3 == 0 : ax.set_ylabel('deg')
            ax.set_aspect('equal')
            cbar = plt.colorbar(mesh, ax=ax, shrink=0.8)
            if n%3 == int(n_axes/rows): cbar.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
    
        #Add title with information about FOV and simulation. 
        fig.text(0.5, 0.95, 'n = {} cm'.format(self.density)+r'$^{-3}$'+'\nSMILE Coords: ({:.2f},{:.2f},{:.2f}),   Aim Point: ({},{},{})'.format(self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2]), ha='center', va='top', fontsize=12)
    
        if save:
            print ('Saved to: ', self.plot_path+'PPMLR_Image_Sections_n_{:.1f}_SMILE_({:.2f},{:.2f},{:.2f})_Target_({:.2f},{:.2f},{:.2f})_{}.png'.format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2], savetag))  
            fig.savefig(self.plot_path+'PPMLR_Image_Sections_n_{:.1f}_SMILE_({:.2f},{:.2f},{:.2f})_Target_({:.2f},{:.2f},{:.2f})_{}.png'.format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2], savetag))
    


        
