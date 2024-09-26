#This will analyse the data in the pickle file for the run_fit_image_target experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

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
        ax.set_xlabel('Target x Position [RE]') 
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
        
