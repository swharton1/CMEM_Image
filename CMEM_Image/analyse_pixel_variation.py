#This will analyse the data in the pickle file for the run_fit_image_pixel experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import string
from matplotlib.ticker import MultipleLocator

class analyse_pixel():
    '''This will contain any plotting functions needed'''
    
    def __init__(self, fname='pixel_variation_output.pkl'):
        '''Sorts out paths and reads in the file.'''
        
        self.fname = fname 
        self.plot_path = os.environ.get('PLOT_PATH')+'pixel_variation/' 
        self.pickle_path = os.environ.get('PICKLE_PATH')
        
        self.target_dict = self.read_pickle(self.plot_path+self.fname)
        
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_pixel(self):
        '''This will plot a graph of how the subsolar magnetopause radii 
        determined by CMEM and extracted from the PPMLR simulation vary
        with the number of pixels. This will show us the systematic error.'''
        
        #Extract key values. 
        maxIx = self.target_dict['maxIx']
        maxdIx = self.target_dict['maxdIx']
        f25 = self.target_dict['f.25'] 
        n_pixels = np.array(self.target_dict['n_pixels'])
        m_pixels = np.array(self.target_dict['m_pixels']) 
        cmem_mp = self.target_dict['cmem_mp'] 
        inout = self.target_dict['inout']
        smile_loc = self.target_dict['smile_loc'] 
        density = self.target_dict['density'] 
        
        fig = plt.figure(figsize=(6,4))
        ax = fig.add_subplot(111) 
        fig.subplots_adjust(bottom=0.20)
        
        #Add labels. 
        ax.set_xlabel('M pixel number\n Total Number of Pixels') 
        ax.set_ylabel('Subsolar Magnetopause Position [RE]')
        
        #Get array for total number of pixels. 
        self.total_pixels = n_pixels*m_pixels 
        
        #Add on the PPMLR values. 
        #xlims = [self.total_pixels[0], self.total_pixels[-1]] 
        xlims = [m_pixels[0]-5, m_pixels[-1]+5]
        ax.plot(xlims, [maxIx, maxIx], 'b-', label='maxIx')
        ax.plot(xlims, [f25, f25], 'g-', label='f.25')  
        ax.plot(xlims, [maxdIx, maxdIx], 'r-', label='maxdIx')
        
        #Add on the CMEM determinations. 
        ax.plot(m_pixels, cmem_mp, 'k-', label='CMEM', marker='x') 
        
        ax.set_xticks(m_pixels) 
        xticks = ["{}\n{}".format(m_pixels[x], self.total_pixels[x]) for x in range(len(m_pixels))]
        ax.set_xticklabels(xticks)
        ax.set_xlim(0,55)
        ax.set_yticks(np.linspace(7.8,8.8,11))
        ax.legend(loc='best')
        ax.grid()
        
        ax.set_title('SMILE = ({:.2f},{:.2f},{:.2f}), n = {:.2f} cm'.format(smile_loc[0], smile_loc[1], smile_loc[2], density)+r'$^{-3}$')
        
        #Save the plot. 
        fig.savefig(self.plot_path+'pixel_variation_analysis.png') 

class analyse_pixel_all():
    '''This will create a master plot for all solar wind densities.'''
    
    def __init__(self):
        '''Sorts out paths and reads in the file.'''
        
        self.filenames = ['pixel_variation_output_n_5.0.pkl', 
            'pixel_variation_output_n_7.5.pkl',
            'pixel_variation_output_n_12.3.pkl',
            'pixel_variation_output_n_20.0.pkl',
            'pixel_variation_output_n_25.0.pkl',
            'pixel_variation_output_n_35.0.pkl'] 
        
        
        self.plot_path = os.environ.get('PLOT_PATH')+'pixel_variation/' 
        
        self.data = []
        for f in self.filenames: 
            self.data.append(self.read_pickle(self.plot_path+f))
        
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict        
        
    def make_plot(self):
        '''This will make the plot that appears in the paper.'''
        plt.close("all")
        fig = plt.figure(figsize=(6,8))
        
        for d, data in enumerate(self.data):
        
            #Add an axis. 
            ax = fig.add_subplot(6,1,d+1)
            ax.set_ylabel(r'$r_{mp}$ '+r'$(R_E)$')
            
            #Add three definitions. 
            ax.plot([data['m_pixels'][0], data['m_pixels'][-1]], [data['maxIx'], data['maxIx']], 'b')
            ax.plot([data['m_pixels'][0], data['m_pixels'][-1]], [data['f.25'], data['f.25']], 'g')
            ax.plot([data['m_pixels'][0], data['m_pixels'][-1]], [data['maxdIx'], data['maxdIx']], 'r')
            
            #Add CMEM. 
            ax.plot(data['m_pixels'], data['cmem_mp'], 'k', marker='x')
            
            #Set xlims. 
            ax.set_xlim(0,55) 
            
            #Add density label. 
            ax.text(0.99, 0.95, '{:.2f} cm'.format(data['density'])+r'$^{-3}$', fontsize=8, transform=ax.transAxes, ha='right', va='top', rotation=0, bbox=dict(boxstyle='round', facecolor='w', alpha=1, edgecolor='none'))
            
            #Subplot labels 
            letters = string.ascii_lowercase
            ax.text(0.01, 0.95, '({})'.format(letters[d]), fontsize=8, transform=ax.transAxes, ha='left', va='top', rotation=0, bbox=dict(boxstyle='round', facecolor='w', alpha=1, edgecolor='none'))
            
            if d == 5:
                xticklabels = ['{}\n{}'.format(data['m_pixels'][i], data['m_pixels'][i]*data['n_pixels'][i]) for i in range(len(data['m_pixels']))]
                ax.set_xlabel('m pixel number\nTotal pixel number') 
            else:
                xticklabels = []
            ax.set_xticks(data['m_pixels'], labels=xticklabels, fontsize=10) 
            
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.set_ylim(data['f.25']-1, data['f.25']+1)
            ax.grid(which='both')
            
            #Add title to top plot. 
            if d == 0: 
                ax.set_title('Variation of Subsolar Magnetopause Position\nwith pixel number\nSMILE = ({:.2f},{:.2f},{:.2f}), Aim Point = ({:.2f},{:.2f},{:.2f})'.format(*data['smile_loc'], *data['target_loc']), fontsize=10)

        #Save the plot. 
        fig.savefig(self.plot_path+'pixel_analysis_all.png')           
