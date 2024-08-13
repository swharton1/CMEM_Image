#This class will visualise the output of the fit image class. 

import numpy as np
import matplotlib.pyplot as plt 
import pickle 
import os

from . import get_names_and_units as gnau

class analyse_fit():
	'''This class will contain the plotting functions that 
	the CMEM/visualise_models.py analysis class used.'''
	
	def __init__(self, filename='fit_image_n_5.0_SMILE_-10_-30_0_Target_10_0_0_nxm_100_50_cmem_normalised_im2_.pkl', model='cmem'): 
		'''This takes in the filename for the pickle file.''' 
		
		
		self.current_model = model 
		self.filename = filename
		
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
		
	def __repr__(self):
		return f"analyse model object."
    
	def read_pickle(self, filename):
		'''This will read a single pickle file. '''

		with open(filename, 'rb') as f: 
			pickle_dict = pickle.load(f)
		return pickle_dict
	
	
	def plot_change_in_parameters(self, save=False, savetag=""):
		'''This will plot how the parameters changed over the course of the
		optimisation procedure. '''

		fig = plt.figure(figsize=(8,8))
		fig.subplots_adjust(hspace=0.4, wspace=0.4, top=0.85)
		fig.text(0.5, 0.9, "Parameter variation with optimisation\nn = {:.2f} cm".format(self.model['density'])+r"$^{-3}$"+", SMILE = ({},{},{}), Target = ({},{},{}), nxm = {}x{}".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2], self.model['n_pixels'], self.model['m_pixels'])+" \nOptimisation Time = {:.1f}s\nModel = {}".format(self.model['opt time'], self.image_tag.capitalize()), ha="center")
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
			fig.savefig(self.plot_path+"fitted_images/{}_parameter_changes_{}.png".format(self.filename, savetag))

		self.fig_param = fig 


	def plot_images(self, cmap='hot', vmin=-8, vmax=-4, levels=100, los_max=12, save=False, savetag=""):
		'''This will plot the final model image alongside the PPMLR image.''' 
		
		fig = plt.figure(figsize=(8,5))
		fig.subplots_adjust(left=0.05, wspace=0.2, bottom=0.20) 
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
		theta_pixels = -np.rad2deg(theta_pixels)
		phi_pixels = -np.rad2deg(phi_pixels)
        
		# Get contour levels. 
		#levels = np.linspace(vmin, vmax, levels+1)
        
		mesh1 = ax1.pcolormesh(phi_pixels, theta_pixels, self.model['ppmlr los intensity'], cmap=cmap, vmin=0, vmax=los_max)
		ax1.set_title("PPMLR Image from SMILE\nSMILE = ({},{},{}), Target = ({},{},{})".format(self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2]), fontsize=10)
		ax1.set_xlabel('deg')
		ax1.set_ylabel('deg')
		ax1.set_aspect('equal')
		cbar1 = plt.colorbar(mesh1, ax=ax1, shrink=0.8)
		cbar1.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
		
		#Now add the model image. 
		mesh2 = ax2.pcolormesh(phi_pixels, theta_pixels, self.model['model los intensity'], cmap=cmap, vmin=0, vmax=los_max)
		ax2.set_title("{} Image from SMILE\nSMILE = ({},{},{}), Target = ({},{},{})".format(self.image_tag, self.model['smile_loc'][0], self.model['smile_loc'][1], self.model['smile_loc'][2], self.model['target_loc'][0], self.model['target_loc'][1], self.model['target_loc'][2]), fontsize=10)
		
		ax2.set_xlabel('deg')
		ax2.set_ylabel('deg')
		ax2.set_aspect('equal')
		cbar2 = plt.colorbar(mesh2, ax=ax2, shrink=0.8)
		cbar2.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
		
		
		# Add a label to show the model parameters. 
		label = ""
		for p,pval in enumerate(self.model['params best nm']):
				pv = pval 
				label += "{}={} {}, ".format(self.parameter_names[p], self.sig_figs(pv,3), self.parameter_units[p])
				if len(self.parameter_names)//2 == p+1:
					label += "\n"
		
		fig.text(0.5, 0.02, label, ha='center')
		
		if save: 
			fig.savefig(self.plot_path+"fitted_images/{}_images_{}.png".format(self.filename, savetag))

		self.fig_param = fig 
		
		
		
		
		
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