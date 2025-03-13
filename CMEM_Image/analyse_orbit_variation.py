#This will analyse the data in the pickle file for the run_fit_image_orbit experiment. 

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge, Polygon, Circle

class analyse_orbit():
    '''This will contain any plotting functions needed'''
    
    def __init__(self, fname='orbital_limb_output_newsmile.pkl', orbit_folder_name='orbit_10/', ppmlr_folder='S05D20V400B0000-05rad.fits/'):
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
    
    def plot_orbit(self, add_cost=False, figname='Orbit_04_analysis.png', aim_ylim=[4,12], smile_xlim=[-10,15], smile_ylim=[-10,10], smile_zlim=[-5,22]):
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
        xgse = np.array(self.target_dict['x_gse'])
        ygse = np.array(self.target_dict['y_gse'])
        zgse = np.array(self.target_dict['z_gse']) 
        
        target_loc = np.array(self.target_dict['target_loc']).T[0]
        min_cost = np.array(self.target_dict['min_cost'])
        #tilt = np.array(self.target_dict['tilt']) 
        
        #Get time axis in hours. 
        delta_time = np.arange(time.size) 
        
        #Select data where s/c is inside MP and where it is outside MP. 
        inside = np.where(inout == 'in')
        outside = np.where(inout == 'out')
        
        #Get time resolution. 
        #time_res = time[1] - time[0] 
        time_res=1
        #Identify times inside and outside magnetopause.  
        intime = delta_time[inside]
        outtime = delta_time[outside] 
        
        #Make the figure. 
        if add_cost:
            fig = plt.figure(figsize=(6,6))
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
            fig = plt.figure(figsize=(6,8))
            gs = GridSpec(15,3,figure=fig)
            fig.subplots_adjust(bottom=0.10, wspace=0.7, right=0.9, hspace=0.2, top=0.9)
            ax = fig.add_subplot(gs[0:3, 0:3]) 
            #ax2 = fig.add_subplot(gs[1,0:2])
            ax3 = fig.add_subplot(gs[4:7,0:3])
            
            #Add small plots to show SMILE position in GSE.  
            ax4 = fig.add_subplot(gs[8:11,0])
            ax5 = fig.add_subplot(gs[8:11,1])
            ax6 = fig.add_subplot(gs[8:11,2])
            
            #Add small plots to show SMILE position in GSM.  
            ax7 = fig.add_subplot(gs[12:,0])
            ax8 = fig.add_subplot(gs[12:,1])
            ax9 = fig.add_subplot(gs[12:,2])
                      
            #Add labels. 
            ax.set_xlabel('Observing Time (hours)', fontsize=8) 
            ax3.set_xlabel('Observing Time (hours)', fontsize=8) 
            ax.set_ylabel('Subsolar Magnetopause\n Position [RE]', fontsize=8)
            #ax2.set_ylabel('GSM Spacecraft \nPosition [RE]') 
            ax3.set_ylabel('Aim Point [RE]', fontsize=8) 
            
            #Sort axes for orbital plots. 
            ax4.set(xlim=smile_xlim, ylim=smile_zlim)
            ax5.set(xlim=smile_xlim, ylim=smile_ylim)
            ax6.set(xlim=smile_ylim, ylim=smile_zlim) 
            ax4.set_xlabel(r'$X_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            ax4.set_ylabel(r'$Z_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            ax5.set_xlabel(r'$X_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            ax5.set_ylabel(r'$Y_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            ax6.set_xlabel(r'$Y_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            ax6.set_ylabel(r'$Z_{GSE}$'+' '+r'$[R_E]$', fontsize=8)
            
            #Sort axes for orbital plots. 
            ax7.set(xlim=smile_xlim, ylim=smile_zlim)
            ax8.set(xlim=smile_xlim, ylim=smile_ylim)
            ax9.set(xlim=smile_ylim, ylim=smile_zlim) 
            ax7.set_xlabel(r'$X_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            ax7.set_ylabel(r'$Z_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            ax8.set_xlabel(r'$X_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            ax8.set_ylabel(r'$Y_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            ax9.set_xlabel(r'$Y_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            ax9.set_ylabel(r'$Z_{GSM}$'+' '+r'$[R_E]$', fontsize=8)
            
            
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
        #r = np.sqrt(xgsm**2 + ygsm**2 + zgsm**2) 
        #ax2.plot(delta_time, xgsm, 'b-', label='x', marker='x')
        #ax2.plot(delta_time, ygsm, 'r-', label='y', marker='x')
        #ax2.plot(delta_time, zgsm, 'g-', label='z', marker='x') 
        #ax2.plot(delta_time, r, 'k-', label='r', marker='x')
        
        #Add on aim point. 
        ax3.plot(xlims, [maxIx, maxIx], 'b-', label='maxIx')
        ax3.plot(xlims, [f25, f25], 'g-', label='f.25')  
        ax3.plot(xlims, [maxdIx, maxdIx], 'r-', label='maxdIx')
        
        ax3.plot(delta_time, target_loc, 'k-', marker='x') 
        
        
        
        
        #ax3b.plot(time, np.rad2deg(tilt), 'b-', marker='x') 
        
        #Add orbital positions to right hand plots. 
        ax4.grid()
        ax5.grid()
        ax6.grid()
        ax4.plot(xgse, zgse, 'k', lw=0.5)
        ax5.plot(xgse, ygse, 'k', lw=0.5)
        ax6.plot(ygse, zgse, 'k', lw=0.5) 
        self.make_earth(ax4, rotation=-90)
        self.make_earth(ax5, rotation=-90)
        circle = Circle((0,0), 1, facecolor='w', edgecolor='navy')
        ax6.add_patch(circle)
        ax4.set_aspect('equal')
        ax5.set_aspect('equal')
        ax6.set_aspect('equal')
        
        #Add orbital positions to right hand plots. 
        ax7.grid()
        ax8.grid()
        ax9.grid()
        ax7.plot(xgsm, zgsm, 'k', lw=0.5)
        ax8.plot(xgsm, ygsm, 'k', lw=0.5)
        ax9.plot(ygsm, zgsm, 'k', lw=0.5) 
        self.make_earth(ax7, rotation=-90)
        self.make_earth(ax8, rotation=-90)
        circle = Circle((0,0), 1, facecolor='w', edgecolor='navy')
        ax9.add_patch(circle)
        ax7.set_aspect('equal')
        ax8.set_aspect('equal')
        ax9.set_aspect('equal')
        
        #Add a mark every five hours. 
        ax4.scatter(xgse[::5], zgse[::5], c='k', marker='x', s=10, lw=0.5)
        ax5.scatter(xgse[::5], ygse[::5], c='k', marker='x', s=10, lw=0.5)
        ax6.scatter(ygse[::5], zgse[::5], c='k', marker='x', s=10, lw=0.5) 
        
        #Add a mark every five hours. 
        ax7.scatter(xgsm[::5], zgsm[::5], c='k', marker='x', s=10, lw=0.5)
        ax8.scatter(xgsm[::5], ygsm[::5], c='k', marker='x', s=10, lw=0.5)
        ax9.scatter(ygsm[::5], zgsm[::5], c='k', marker='x', s=10, lw=0.5) 
        
        
        #Add time in hours since start of observing time. 
        xmax4 = xgse[zgse==zgse.max()]
        xmax5 = xgse[ygse==ygse.max()]
        ymax6 = ygse[zgse==zgse.max()]
        
        for t in range(len(xgse[::5])):
             
            ha4='right' if xgse[::5][t] < xmax4 else 'left'
            ax4.text(xgse[::5][t], zgse[::5][t], delta_time[::5][t], fontsize=8, ha=ha4)
            
            ha5='right' if xgse[::5][t] < xmax5 else 'left'
            ax5.text(xgse[::5][t], ygse[::5][t], delta_time[::5][t], fontsize=8, ha=ha5)
            
            ha6='right' if ygse[::5][t] < ymax6 else 'left'
            ax6.text(ygse[::5][t], zgse[::5][t], delta_time[::5][t], fontsize=8, ha=ha6)
            
        #Add time in hours since start of observing time. 
        xmax7 = xgsm[zgsm==zgsm.max()]
        xmax8 = xgsm[ygsm==ygsm.max()]
        ymax9 = ygsm[zgsm==zgsm.max()]
        
        for t in range(len(xgsm[::5])):
             
            ha7='right' if xgsm[::5][t] < xmax7 else 'left'
            ax7.text(xgsm[::5][t], zgsm[::5][t], delta_time[::5][t], fontsize=8, ha=ha7)
            
            ha8='right' if xgsm[::5][t] < xmax8 else 'left'
            ax8.text(xgsm[::5][t], ygsm[::5][t], delta_time[::5][t], fontsize=8, ha=ha8)
            
            ha9='right' if ygsm[::5][t] < ymax9 else 'left'
            ax9.text(ygsm[::5][t], zgsm[::5][t], delta_time[::5][t], fontsize=8, ha=ha9)
            
        #Shrink fontsize for tickmarks. 
        fontsize=8
        ax.tick_params(labelsize=fontsize)
        ax3.tick_params(labelsize=fontsize)
        ax4.tick_params(labelsize=fontsize)
        ax5.tick_params(labelsize=fontsize)
        ax6.tick_params(labelsize=fontsize)
        ax7.tick_params(labelsize=fontsize)
        ax8.tick_params(labelsize=fontsize)
        ax9.tick_params(labelsize=fontsize)
        
        
        
        
        if add_cost:
            #Add the cost values on. 
            ax3.plot(time, min_cost, 'k') 
        
        #Loop through the intime array to show when SMILE is thought to be inside the magnetopause. 
        for i, ival in enumerate(intime):
            
            #Add patch to top figure. 
            rect1 = Rectangle([ival-time_res/2,0], height=15, width=time_res, color='lightgrey', alpha=0.4)
            ax.add_patch(rect1)
            
            #Add patch to second figure. 
            #rect2 = Rectangle([ival-time_res/2, -20], height=50, width=time_res, color='lightgrey', alpha=0.4)
            #ax2.add_patch(rect2) 
            
            #Add patch to third figure. 
            rect3 = Rectangle([ival-time_res/2, 0], height=50, width=time_res, color='lightgrey', alpha=0.4)
            ax3.add_patch(rect3) 
        
        
        #if add_cost:
        #    rect5 = Rectangle([time[0], 0], height=10, width=exit_time, color='gray', alpha=0.4)
        #    rect6 = Rectangle([entry_time, 0], height=10, width=time[-1]-entry_time, color='gray', alpha=0.4, linecolor=None)
        
       
        #if add_cost:
        #    ax3.add_patch(rect5)
        #    ax3.add_patch(rect6) 
        
        #Reset your limits. 
        ax.set_xlim(0,41) #HARDCODED HERE SO ORBITS MATCH. 
        ax.set_ylim(7,10)
        #ax.legend(loc='best', fontsize=8)
        ax.grid() 
        
        #Add labels to side. 
        ax.text(1.01, 0.99, 'max Ix', color='b', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.89, 'f.25', color='g', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.79, 'max dIx', color='r', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        ax.text(1.01, 0.69, 'CMEM', color='k', ha='left', va='top', transform=ax.transAxes, fontsize=8)
        
        #Reset your limits. 
        #ax2.set_xlim(xlims) 
        #ax2.set_ylim(-20,30)
        #ax2.legend(loc='best', fontsize=8)
        #ax2.grid()
        
        #Add labels to side. 
        #ax2.text(1.01, 0.99, 'x', color='b', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        #ax2.text(1.01, 0.89, 'y', color='r', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        #ax2.text(1.01, 0.79, 'z', color='g', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        #ax2.text(1.01, 0.69, 'r', color='k', ha='left', va='top', transform=ax2.transAxes, fontsize=8)
        
        ax3.set_xlim(0,41) #HARDCODED HERE SO ORBITS MATCH. 
        ax3.set_ylim(aim_ylim)
        ax3.grid()
        
        #Add a title. 
        start = time[0]
        end = time[-1]
        ax.set_title(start.strftime("%Y%m%d %H:%M")+' - '+end.strftime("%Y%m%d %H:%M"), ha='center')
        
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
            print ('Saved at: ', self.plot_path+figname)
            fig.savefig(self.plot_path+figname) 

    def make_earth(self, ax, rotation=0):
        '''This will add a little plot of the Earth on top for reference. '''

        # Add white circle first. 
        r=1
        circle = Circle((0,0), r, facecolor='w', edgecolor='navy')
        ax.add_patch(circle)

        # Add nightside. 
        theta2 = np.arange(181)-180+rotation
        xval2 = np.append(r*np.cos(theta2*(np.pi/180)),0)
        yval2 = np.append(r*np.sin(theta2*(np.pi/180)),0)
        verts2 = [[xval2[i],yval2[i]] for i in range(len(xval2))]
        
        polygon2 = Polygon(verts2, closed=True, edgecolor='navy', facecolor='navy', alpha=1) 
        ax.add_patch(polygon2)       
