#This will analyse the data in the pickle file for the run_fit_image_orbit experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from matplotlib.patches import Rectangle

class analyse_orbit():
    '''This will contain any plotting functions needed'''
    
    def __init__(self, fname='orbital_limb_output_corrected.pkl', orbit_folder_name='20260402-20260403/', ppmlr_folder='S05D20V400B0000-05rad.fits/'):
        '''Sorts out paths and reads in the file.'''
        
        self.fname = fname 
        self.plot_path = os.environ.get('PLOT_PATH')+'fitted_orbit_sim_limb/'+orbit_folder_name+ppmlr_folder 
        self.pickle_path = os.environ.get('PICKLE_PATH')
        
        self.target_dict = self.read_pickle(self.plot_path+self.fname)
        
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_orbit(self, add_cost=False, figname='Orbit_1_analysis.png', aim_ylim=[4,12]):
        '''This will plot a graph of how the subsolar magnetopause radii 
        determined by CMEM and extracted from the PPMLR simulation vary
        with the number of pixels. This will show us the systematic error.'''
        
        #Extract key values. 
        maxIx = self.target_dict['maxIx']
        maxdIx = self.target_dict['maxdIx']
        f25 = self.target_dict['f.25']
        n_pixels = np.array(self.target_dict['n_pixels'])
        m_pixels = np.array(self.target_dict['m_pixels']) 
        cmem_mp = np.array(self.target_dict['cmem_mp']) 
        inout = np.array(self.target_dict['inout'])
        smile_loc = np.array(self.target_dict['smile_loc']) 
        density = self.target_dict['density'] 
        time = np.array(self.target_dict['t_list'])
        xgsm = np.array(self.target_dict['x_gsm'])
        ygsm = np.array(self.target_dict['y_gsm'])
        zgsm = np.array(self.target_dict['z_gsm']) 
        target_loc = np.array(self.target_dict['target_loc']).T[0]
        min_cost = np.array(self.target_dict['min_cost'])
        tilt = np.array(self.target_dict['tilt']) 
        
        #Get time axis in hours. 
        delta_time = np.arange(time.size) 
        
        #Select data where s/c is inside MP and where it is outside MP. 
        inside = np.where(inout == 'in')
        outside = np.where(inout == 'out')
        
        #Get time resolution. 
        #time_res = time[1] - time[0] 
        time_res=1
        #Identify midpoint time for exit and entry to magnetosphere. 
        intime = delta_time[inside]
        outtime = delta_time[outside] 
        
        #Make the figure. 
        if add_cost:
            fig = plt.figure(figsize=(6,8))
            fig.subplots_adjust(bottom=0.10)
            ax = fig.add_subplot(311) 
            ax2 = fig.add_subplot(312)
            ax3 = fig.add_subplot(313)
            
            #Add labels. 
            ax3.set_xlabel('Time Since Periapsis (hours)') 
            ax.set_ylabel('Subsolar Magnetopause Position [RE]')
            ax2.set_ylabel('Spacecraft Position [RE]') 
            ax3.set_ylabel('Min. Cost [keV cm'+r'$^{-3}$'+' s'+r'$^{-1}$'+' sr'+r'$^{-1}$')
        
        else:
            fig = plt.figure(figsize=(6,6))
            fig.subplots_adjust(bottom=0.10)
            ax = fig.add_subplot(311) 
            ax2 = fig.add_subplot(312)
            ax3 = fig.add_subplot(313)
                        
            #Add labels. 
            ax3.set_xlabel('Time (hours)') 
            ax.set_ylabel('Subsolar Magnetopause\n Position [RE]')
            ax2.set_ylabel('Spacecraft \nPosition [RE]') 
            ax3.set_ylabel('Aim Point [RE]') 
            #ax3b = ax3.twinx()
            #ax3b.set_ylabel('Tilt Angle [deg]', color='b') 
            
        
        #Add on the PPMLR values. 
        #xlims = [self.total_pixels[0], self.total_pixels[-1]] 
        xlims = [delta_time[0], delta_time[-1]]
        ax.plot(xlims, [maxIx, maxIx], 'b-', label='maxIx')
        ax.plot(xlims, [f25, f25], 'g-', label='f.25')  
        ax.plot(xlims, [maxdIx, maxdIx], 'r-', label='maxdIx')
        
        #Add on the CMEM determinations. 
        ax.plot(delta_time, cmem_mp, 'k-', label='CMEM', marker='x') 
        

        
        #Add orbit information onto plot below. 
        r = np.sqrt(xgsm**2 + ygsm**2 + zgsm**2) 
        ax2.plot(delta_time, xgsm, 'b-', label='x', marker='x')
        ax2.plot(delta_time, ygsm, 'r-', label='y', marker='x')
        ax2.plot(delta_time, zgsm, 'g-', label='z', marker='x') 
        ax2.plot(delta_time, r, 'k-', label='r', marker='x')
        
        #Add on aim point and tilt angle. 
        ax3.plot(delta_time, target_loc, 'k-', marker='x') 
        #ax3b.plot(time, np.rad2deg(tilt), 'b-', marker='x') 
        
        
        if add_cost:
            #Add the cost values on. 
            ax3.plot(time, min_cost, 'k') 
        
        #Loop through the intime array to show when SMILE is thought to be inside the magnetopause. 
        for i, ival in enumerate(intime):
            
            #Add patch to top figure. 
            rect1 = Rectangle([ival,0], height=15, width=time_res, color='gray', alpha=0.4)
            ax.add_patch(rect1)
            
            #Add patch to second figure. 
            rect2 = Rectangle([ival, -20], height=50, width=time_res, color='gray', alpha=0.4)
            ax2.add_patch(rect2) 
            
            #Add patch to third figure. 
            rect3 = Rectangle([ival, 0], height=50, width=time_res, color='gray', alpha=0.4)
            ax3.add_patch(rect3) 
        
        
        #if add_cost:
        #    rect5 = Rectangle([time[0], 0], height=10, width=exit_time, color='gray', alpha=0.4)
        #    rect6 = Rectangle([entry_time, 0], height=10, width=time[-1]-entry_time, color='gray', alpha=0.4, linecolor=None)
        
       
        #if add_cost:
        #    ax3.add_patch(rect5)
        #    ax3.add_patch(rect6) 
        
        #Reset your limits. 
        ax.set_xlim(xlims) 
        ax.set_ylim(6.5,10.5)
        #ax.legend(loc='best', fontsize=8)
        ax.grid() 
        
        #Add labels to side. 
        ax.text(1.01, 0.99, 'max Ix', color='b', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.89, 'f.25', color='g', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.79, 'max dIx', color='r', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.69, 'CMEM', color='k', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        
        #Reset your limits. 
        ax2.set_xlim(xlims) 
        ax2.set_ylim(-20,30)
        #ax2.legend(loc='best', fontsize=8)
        ax2.grid()
        
        #Add labels to side. 
        ax2.text(1.01, 0.99, 'x', color='b', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        ax2.text(1.01, 0.89, 'y', color='r', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        ax2.text(1.01, 0.79, 'z', color='g', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        ax2.text(1.01, 0.69, 'r', color='k', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        
        ax3.set_xlim(xlims)
        ax3.set_ylim(aim_ylim)
        ax3.grid()
        
        #Add a title. 
        start = time[0]
        end = time[-1]
        fig.text(0.5, 0.9, start.strftime("%Y%m%d %H:%M")+' - '+start.strftime("%Y%m%d %H:%M"), ha='center')
        
        if add_cost:
            #Reset your limits. 
            ax3.set_xlim(xlims) 
            ax3.set_ylim(0,5)
            #ax3.legend(loc='best', fontsize=8)
            ax3.grid()
        
        if add_cost:
            fig.savefig(self.plot_path+'orbit_variation_analysis_with_cost.png')
        else:
            #Save the plot. 
            fig.savefig(self.plot_path+figname) 
        
