#This will analyse the data in the pickle file for the run_fit_image_density experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

class analyse_density():
	'''This will contain any plotting functions needed'''
	
	def __init__(self, fname='density_variation_output.pkl'):
		'''Sorts out paths and reads in the file.'''
		
		self.fname = fname 
		self.plot_path = os.environ.get('PLOT_PATH')+'density_variation/' 
		self.pickle_path = os.environ.get('PICKLE_PATH')
		
		self.target_dict = self.read_pickle(self.plot_path+self.fname)
		
	def read_pickle(self, filename):
		'''This will read a single pickle file. '''

		with open(filename, 'rb') as f: 
			pickle_dict = pickle.load(f)
		return pickle_dict
	
	def plot_density(self):
		'''This will plot a graph of how the subsolar magnetopause radii 
		determined by CMEM and extracted from the PPMLR simulation vary
		with the density. This will show us the systematic error.'''
		
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
		ax.set_xlabel('Solar Wind Density (cm'+r'$^{-3}$)') 
		ax.set_ylabel('Subsolar Magnetopause Position [RE]')
		
		#Get array for total number of pixels. 
		#self.total_pixels = n_pixels*m_pixels 
		
		#Add on the PPMLR values. 
		#xlims = [self.total_pixels[0], self.total_pixels[-1]] 
		#xlims = [m_pixels[0]-5, m_pixels[-1]+5]
		ax.plot(density, maxIx, 'b-', label='maxIx')
		ax.plot(density, f25, 'g-', label='f.25')  
		ax.plot(density, maxdIx, 'r-', label='maxdIx')
		
		#Add on the CMEM determinations. 
		ax.plot(density, cmem_mp, 'k-', label='CMEM', marker='x') 
		
		#ax.set_xticks(m_pixels) 
		#xticks = ["{}\n{}".format(m_pixels[x], self.total_pixels[x]) for x in range(len(m_pixels))]
		#ax.set_xticklabels(xticks)
		ax.set_xlim(0,40)
		#ax.set_yticks(np.linspace(7.8,8.8,11))
		ax.legend(loc='best')
		ax.grid()
		
		ax.set_title('SMILE = ({:.2f},{:.2f},{:.2f})'.format(smile_loc[0], smile_loc[1], smile_loc[2]))
		
		#Save the plot. 
		fig.savefig(self.plot_path+'density_variation_analysis.png') 
		
