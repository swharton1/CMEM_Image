#This will analyse the data in the pickle file for the run_fit_image_target experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import string
from matplotlib.ticker import MultipleLocator


class analyse_target():
    '''This will contain any plotting functions needed'''
    
    def __init__(self, fname='target_variation_output.pkl'):
        '''Sorts out paths and reads in the file.'''
        
        self.fname = fname 
        self.plot_path = os.environ.get('PLOT_PATH')+'target_variation/' 
        self.pickle_path = os.environ.get('PICKLE_PATH')
        
        self.target_dict = self.read_pickle(self.plot_path+self.fname)
        
    def read_pickle(self, filename):
        '''This will read a single pickle file. '''

        with open(filename, 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_targetx(self):
        '''This will plot a graph of how the subsolar magnetopause radii 
        determined by CMEM and extracted from the PPMLR simulation vary
        with the target x position. This will show us the systematic error.'''
        
        #Extract key values. 
        maxIx = self.target_dict['maxIx']
        maxdIx = self.target_dict['maxdIx']
        f25 = self.target_dict['f.25'] 
        targetx = self.target_dict['target_x']
        cmem_mp = self.target_dict['cmem_mp'] 
        inout = self.target_dict['inout']
        smile_loc = self.target_dict['smile_loc'] 
        density = self.target_dict['density'] 
        
        fig = plt.figure(figsize=(6,4))
        ax = fig.add_subplot(111) 
        
        #Add labels. 
        ax.set_xlabel('Aim Point [RE]') 
        ax.set_ylabel('Subsolar Magnetopause Position [RE]')
        
        #Add on the PPMLR values. 
        print (targetx)
        xlims = [targetx[0], targetx[-1]] 
        
        ax.plot(xlims, [maxIx, maxIx], 'b-', label='maxIx')
        ax.plot(xlims, [f25, f25], 'g-', label='f.25')  
        ax.plot(xlims, [maxdIx, maxdIx], 'r-', label='maxdIx')
        
        #Add on the CMEM determinations. 
        ax.plot(targetx, cmem_mp, 'k-', label='CMEM', marker='x') 
        ax.set_xlim([xlims[0]-0.5,xlims[1]+0.5])
        ax.set_ylim(6.5,9)
        ax.set_xticks(np.linspace(6.0,13.0,8))
        ax.minorticks_on()
        #ax.set_yticks(np.linspace(7.8,8.8,11))
        ax.legend(loc='best')
        ax.grid(which='both')
        
        ax.set_title('SMILE = ({:.2f},{:.2f},{:.2f}), n = {:.2f} cm'.format(smile_loc[0], smile_loc[1], smile_loc[2], density)+r'$^{-3}$')
        
        #Save the plot. 
        fig.savefig(self.plot_path+'target_variation_analysis.png') 



class analyse_target_all():
    '''This will create a master plot for all solar wind densities.'''
    
    def __init__(self):
        '''Sorts out paths and reads in the file.'''
        
        self.filenames = ['target_variation_output_n_5.0.pkl', 
            'target_variation_output_n_7.5.pkl',
            'target_variation_output_n_12.3.pkl',
            'target_variation_output_n_20.0.pkl',
            'target_variation_output_n_25.0.pkl',
            'target_variation_output_n_35.0.pkl'] 
        
        
        self.plot_path = os.environ.get('PLOT_PATH')+'target_variation/' 
        
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
            ax.plot([data['target_x'][0], data['target_x'][-1]], [data['maxIx'], data['maxIx']], 'b')
            ax.plot([data['target_x'][0], data['target_x'][-1]], [data['f.25'], data['f.25']], 'g')
            ax.plot([data['target_x'][0], data['target_x'][-1]], [data['maxdIx'], data['maxdIx']], 'r')
            
            #Add CMEM. 
            ax.plot(data['target_x'], data['cmem_mp'], 'k', marker='x')
            
            #Set xlims. 
            ax.set_xlim(5.5,13.5)
            
            #Add density label.  
            ax.text(0.99, 0.95, '{:.2f} cm'.format(data['density'])+r'$^{-3}$', fontsize=8, transform=ax.transAxes, ha='right', va='top', rotation=0, bbox=dict(boxstyle='round', facecolor='w', alpha=1, edgecolor='none'))
            
            #Subplot labels 
            letters = string.ascii_lowercase
            ax.text(0.01, 0.95, '({})'.format(letters[d]), fontsize=8, transform=ax.transAxes, ha='left', va='top', rotation=0, bbox=dict(boxstyle='round', facecolor='w', alpha=1, edgecolor='none'))
            
            #if d == 5:
            #    xticklabels = ['{}\n{}'.format(data['m_pixels'][i], data['m_pixels'][i]*data['n_pixels'][i]) for i in range(len(data['m_pixels']))]
            #    ax.set_xlabel('m pixel number\nTotal pixel number') 
            #else:
            if d < 5: 
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Aim Point '+r'$(R_E)$')
            #ax.set_xticks(data['m_pixels'], labels=xticklabels, fontsize=10) 
            #ax.minorticks_on()
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.set_ylim(data['f.25']-1, data['f.25']+1)
            ax.grid(which='both')
            
            #Add vertical lines. 
            ax.plot([data['maxIx'],data['maxIx']], ax.get_ylim(), 'b', zorder=0)
            ax.plot([data['f.25'],data['f.25']], ax.get_ylim(), 'g', zorder=0)
            ax.plot([data['maxdIx'],data['maxdIx']], ax.get_ylim(), 'r', zorder=0)
            
            
            #Add title to top plot. 
            if d == 0: 
                ax.set_title('Variation of Subsolar Magnetopause Position\nwith Aim Point   SMILE = ({:.2f},{:.2f},{:.2f})'.format(*data['smile_loc'], ), fontsize=10)

        #Save the plot. 
        fig.savefig(self.plot_path+'target_analysis_all.png')     
