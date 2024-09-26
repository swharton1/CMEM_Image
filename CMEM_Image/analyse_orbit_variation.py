#This will analyse the data in the pickle file for the run_fit_image_orbit experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from matplotlib.patches import Rectangle

class analyse_orbit():
    '''This will contain any plotting functions needed'''
    
    def __init__(self, fname='parameter_output.pkl', orbit_folder_name='ra2_rp19_inc70_raan180_omega300/', ppmlr_folder='S05D20V400B0000-05rad.dat/'):
        '''Sorts out paths and reads in the file.'''
        
        self.fname = fname 
        self.plot_path = os.environ.get('PLOT_PATH')+'fitted_orbit_sim/'+orbit_folder_name+ppmlr_folder 
        self.pickle_path = os.environ.get('PICKLE_PATH')
        
        self.target_dict = self.read_pickle(self.plot_path+self.fname)
        
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_orbit(self, add_cost=False):
        '''This will plot a graph of how the subsolar magnetopause radii 
        determined by CMEM and extracted from the PPMLR simulation vary
        with the number of pixels. This will show us the systematic error.'''
        
        #Extract key values. 
        maxIx = self.target_dict['maxIx']
        maxdIx = self.target_dict['maxdIx']
        f25 = self.target_dict['f.25'][0] 
        n_pixels = np.array(self.target_dict['n_pixels'])
        m_pixels = np.array(self.target_dict['m_pixels']) 
        cmem_mp = np.array(self.target_dict['cmem_mp']) 
        inout = np.array(self.target_dict['inout'])
        smile_loc = np.array(self.target_dict['smile_loc']) 
        density = self.target_dict['density'] 
        time = self.target_dict['t_list'] 
        xgsm = self.target_dict['x_gsm']
        ygsm = self.target_dict['y_gsm']
        zgsm = self.target_dict['z_gsm'] 
        target_loc = self.target_dict['target_loc']
        min_cost = np.array(self.target_dict['min_cost'])
        
        
        #Select data where s/c is inside MP and where it is outside MP. 
        inside = np.where(inout == 'in')
        outside = np.where(inout == 'out')
        
        #Get time resolution. 
        time_res = time[1] - time[0] 
        
        #Identify midpoint time for exit and entry to magnetosphere. 
        intime = time[inside]
        outtime = time[outside] 
        
        exit_time = outtime[0] - time_res/2. 
        entry_time = outtime[-1] + time_res/2. 
        
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
            ax = fig.add_subplot(211) 
            ax2 = fig.add_subplot(212)
            
            #Add labels. 
            ax2.set_xlabel('Time Since Periapsis (hours)') 
            ax.set_ylabel('Subsolar Magnetopause Position [RE]')
            ax2.set_ylabel('Spacecraft Position [RE]') 
            
        
        
        #Get array for total number of pixels. 
        #self.total_pixels = n_pixels*m_pixels 
        
        #Add on the PPMLR values. 
        #xlims = [self.total_pixels[0], self.total_pixels[-1]] 
        xlims = [time[0], time[-1]]
        ax.plot(xlims, [maxIx, maxIx], 'b-', label='maxIx')
        ax.plot(xlims, [f25, f25], 'g-', label='f.25')  
        ax.plot(xlims, [maxdIx, maxdIx], 'r-', label='maxdIx')
        
        #Add on the CMEM determinations. 
        ax.plot(time, cmem_mp, 'k-', label='CMEM', marker='x') 
        
        #ax.set_xticks(m_pixels) 
        #xticks = ["{}\n{}".format(m_pixels[x], self.total_pixels[x]) for x in range(len(m_pixels))]
        #ax.set_xticklabels(xticks)
        #ax.set_xlim(0,55)
        
        
        ax.set_title('Target = ({:.2f},{:.2f},{:.2f}), n = {:.2f} cm'.format(target_loc[0], target_loc[1], target_loc[2], density)+r'$^{-3}$')
        
        #Add orbit information onto plot below. 
        r = np.sqrt(xgsm**2 + ygsm**2 + zgsm**2) 
        ax2.plot(time, xgsm, 'b-', label='x', marker='x')
        ax2.plot(time, ygsm, 'r-', label='y', marker='x')
        ax2.plot(time, zgsm, 'g-', label='z', marker='x') 
        ax2.plot(time, r, 'k-', label='r', marker='x')
        
        if add_cost:
            #Add the cost values on. 
            ax3.plot(time, min_cost, 'k') 
        
        
        #Add rectangles in grey to show areas inside the magnetopause. 
        rect1 = Rectangle([time[0], 0], height=15, width=exit_time, color='gray', alpha=0.4)
        rect2 = Rectangle([entry_time, 0], height=15, width=time[-1]-entry_time, color='gray', alpha=0.4)
        rect3 = Rectangle([time[0], -20], height=50, width=exit_time, color='gray', alpha=0.4)
        rect4 = Rectangle([entry_time, -20], height=50, width=time[-1]-entry_time, color='gray', alpha=0.4)
        if add_cost:
            rect5 = Rectangle([time[0], 0], height=10, width=exit_time, color='gray', alpha=0.4)
            rect6 = Rectangle([entry_time, 0], height=10, width=time[-1]-entry_time, color='gray', alpha=0.4)
        
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax2.add_patch(rect3)
        ax2.add_patch(rect4) 
        if add_cost:
            ax3.add_patch(rect5)
            ax3.add_patch(rect6) 
        
        #Reset your limits. 
        ax.set_xlim(xlims) 
        ax.set_ylim(6.5,10.5)
        ax.legend(loc='best', fontsize=8)
        ax.grid()
        
        #Reset your limits. 
        ax2.set_xlim(xlims) 
        ax2.set_ylim(-20,30)
        ax2.legend(loc='best', fontsize=8)
        ax2.grid()
        
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
            fig.savefig(self.plot_path+'orbit_variation_analysis.png') 
        
