#This will vary the orbital position and show what the pointing looks like throughout the orbit using my code. 

import numpy as np
import matplotlib.pyplot as plt
from . import smile_fov_limb as sfl 
from . import load_ephemeris_vlocal as lev
import datetime as dt
import matplotlib.animation as animation 
import matplotlib.gridspec as gridspec
import matplotlib.dates as dates
import os
from matplotlib.patches import Rectangle 

from SXI_Core import get_earth 

class orbit_variation():
    def __init__(self, stime=(2025,10,1), etime=(2025,10,3,3)): 
        '''This takes in any initial parameters, like orbital date.''' 
        
        self.stime = stime
        self.etime = etime 
        self.get_orbit_data()
        
        self.plot_path = os.environ.get('PLOT_PATH') 
        
    def get_orbit_data(self):
        '''This will read in the orbital positions.'''
        
        orbit = lev.orbit(stime=self.stime, etime=self.etime, calc_gsm=True)
        self.orbit = orbit.new_data
    
    
    def get_orbit_parameters(self, p_max=10):
        '''This will loop through the smile objects and investigate their calculations.'''
        
        #Get SMILE locations throughout the orbit. 
        x = self.orbit['x_gse']
        y = self.orbit['y_gse']
        z = self.orbit['z_gse']
        r = np.sqrt(x**2 + y**2 + z**2) 
        t = self.orbit['dtime']
        
        #For each time along the orbit. 
        #print ("[x,y,z],             Aim,   Lx,   bx,   n,              cos(t),   Tilt,  bz/by, sz/sy")
        #for i, ival in enumerate(t):
        #    smile = sfl.smile_limb(smile_loc=[x[i], y[i], z[i]], p_max=p_max) 
        #    print ("[{:.2f},{:.2f},{:.2f}], {:.2f}, {:.2f}, {:.2f},  ({:.2f},{:.2f},{:.2f}),  {:.2f}, {:.2f},    {:.2f}, {:.2f}".format(x[i], y[i], z[i], smile.target_loc[0], smile.L[0], smile.b[0], smile.n[0], smile.n[1], smile.n[2], np.cos(smile.sxi_tilt), np.rad2deg(smile.sxi_tilt), smile.b[2]/smile.b[1], z[i]/y[i]))
    
        
    def animate_orbit(self, p_max=10, xlim=[-8,12], ylim=[-10,10], zlim=[-5,20], aimlim=[0,20]):
        '''This will create an initial FOV plot then animate it to change throughout the full orbit.'''
        
        #Get SMILE locations throughout the orbit. 
        x = self.orbit['x_gse']
        y = self.orbit['y_gse']
        z = self.orbit['z_gse']
        r = np.sqrt(x**2 + y**2 + z**2) 
        t = self.orbit['dtime']
        
        #Get initial SMILE object. 
        smile = sfl.smile_limb(n_pixels=8, m_pixels=4, smile_loc=[x[0], y[0], z[0]], p_max=p_max) 
        
        #Make initial figure. 
        plt.close("all")
        fig = plt.figure(figsize=(10,8))
        fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.05)
        spec = gridspec.GridSpec(nrows=4, ncols=3, figure=fig)
        ax = fig.add_subplot(spec[0:4,0:2], projection='3d') 
        ax2 = fig.add_subplot(spec[0,2]) 
        ax3 = fig.add_subplot(spec[1,2])
        ax4 = fig.add_subplot(spec[2,2])
        #ax5 = fig.add_subplot(spec[3,2])
        ax6 = fig.add_subplot(spec[3,2])
        
        
        #Add the Earth. 
        get_earth.make_earth_3d_2(ax) 
        
        #Add GSE Axes. 
        ax.plot(xlim,[0,0],[0,0],'k', lw=0.5)
        ax.plot([0,0],ylim,[0,0],'k', lw=0.5)
        ax.plot([0,0],[0,0],zlim,'k', lw=0.5)
        
        #Plot SMILE location in blue. 
        smile_pos = ax.plot(x[0], y[0], z[0], c='b', marker='o')[0]
        self.smile_pos = smile_pos
        
        #Add SMILE vector. 
        smile_vect = ax.plot([0, x[0]], [0, y[0]], [0, z[0]], 'b-', label='SMILE')[0] 
        self.smile_vect = smile_vect
        
        #Add SMILE orbit. 
        smile_orbit = ax.plot(x[0], y[0], z[0], c='k', lw=0.5)[0]
        
        #Add projections of SMILE onto the panels. 
        smile_x = ax.plot(xlim[0], y[0], z[0], c='b', marker='o', alpha=0.5)[0]
        smile_y = ax.plot(x[0], ylim[0], z[0], c='b', marker='o', alpha=0.5)[0]
        smile_z = ax.plot(x[0], y[0], zlim[0], c='b', marker='o', alpha=0.5)[0]
        
        smile_x_orbit = ax.plot(xlim[0], y[0], z[0], c='k', alpha=0.5, lw=0.5)[0]
        smile_y_orbit = ax.plot(x[0], ylim[0], z[0], c='k', alpha=0.5, lw=0.5)[0]
        smile_z_orbit = ax.plot(x[0], y[0], zlim[0], c='k', alpha=0.5, lw=0.5)[0]
        
        
        #Add Look vector. 
        look_vect = ax.plot([x[0], x[0]+smile.L[0]], [y[0], y[0]+smile.L[1]], [z[0], z[0]+smile.L[2]], 'r-', label='Look')[0]
        
        #Add Target vector (aim). 
        aim_vect = ax.plot([0, smile.target_loc[0]], [0, smile.target_loc[1]], [0, smile.target_loc[2]], 'g-', label='Target')[0] 
        
        #Add b vector. 
        b_vect = ax.plot([0, smile.b[0]], [0, smile.b[1]], [0, smile.b[2]], 'c-', label='b')[0] 
        
        #Add the FOV boundaries. 
        #For corner pixels only. 
        c1 = ax.plot(smile.xpos[0][0], smile.ypos[0][0], smile.zpos[0][0], 'k', lw=1)[0]
        c2 = ax.plot(smile.xpos[0][-1], smile.ypos[0][-1], smile.zpos[0][-1], 'k', lw=1)[0]
        c3 = ax.plot(smile.xpos[-1][0], smile.ypos[-1][0], smile.zpos[-1][0], 'k', lw=1)[0]
        c4 = ax.plot(smile.xpos[-1][-1], smile.ypos[-1][-1], smile.zpos[-1][-1], 'k', lw=1)[0]
        
        c5 = ax.plot([smile.xpos[0][0][-1],smile.xpos[0][-1][-1]], [smile.ypos[0][0][-1],smile.ypos[0][-1][-1]], [smile.zpos[0][0][-1],smile.zpos[0][-1][-1]], 'k', lw=1)[0]
        c6 = ax.plot([smile.xpos[0][-1][-1],smile.xpos[-1][-1][-1]], [smile.ypos[0][-1][-1],smile.ypos[-1][-1][-1]], [smile.zpos[0][-1][-1],smile.zpos[-1][-1][-1]], 'k', lw=1)[0]
        c7 = ax.plot([smile.xpos[-1][-1][-1],smile.xpos[-1][0][-1]], [smile.ypos[-1][-1][-1],smile.ypos[-1][0][-1]], [smile.zpos[-1][-1][-1],smile.zpos[-1][0][-1]], 'k', lw=1)[0]
        c8 = ax.plot([smile.xpos[-1][0][-1],smile.xpos[0][0][-1]], [smile.ypos[-1][0][-1],smile.ypos[0][0][-1]], [smile.zpos[-1][0][-1],smile.zpos[0][0][-1]], 'k', lw=1)[0]
        #c1, c2, c3, c4, c5, c6, c7, c8 = self.add_fov_boundaries(ax, smile, lw=1) 
        
        #Add the FOV Rectangle. 
        #rect_obj = self.add_fov_rectangle(ax, smile, color='gray') 
        
        # For each pixel: 
        #for i in range(smile.n_pixels):
        #    for j in range(smile.m_pixels):
        #        ax.plot(smile.xpos[i][j], smile.ypos[i][j], smile.zpos[i][j], 'k', lw=0.2)
        
        #Add title to show smile location. 
        title = ax.set_title('SMILE: ({:.2f},{:.2f},{:.2f})\nAim: ({:.2f},{:.2f},{:.2f})\n{} - {} '.format(smile.smile_loc[0], smile.smile_loc[1], smile.smile_loc[2], smile.target_loc[0], smile.target_loc[1], smile.target_loc[2], t[0].strftime('%Y%m%d %H:%M'), t[-1].strftime('%Y%m%d %H:%M')))
        
        ax.set_xlabel(r'$x_{GSE}$')
        ax.set_ylabel(r'$y_{GSE}$')
        ax.set_zlabel(r'$z_{GSE}$')
        ax.set(xlim=xlim, ylim=ylim, zlim=zlim)
               
        ax.set_aspect('equal')
        ax.view_init(elev=30, azim=45)
        
        #Sort out aim point axis below. 
        aim_data = ax2.plot(t[0], smile.target_loc[0], 'k')[0]
        ax2.set_ylabel('Aim Point [RE]')
        #ax2.set_xlabel('UT') 
        ax2.set_xlim(t[0], t[-1]) 
        ax2.set_ylim(aimlim)
        
        #Format time axis. 
        t_locate_major = dates.HourLocator(interval = 12)
        t_locate_minor = dates.HourLocator(interval = 2)
        t_form = dates.DateFormatter('%H:%M')
        ax2.xaxis.set_major_locator(t_locate_major)
        ax2.xaxis.set_minor_locator(t_locate_minor)
        ax2.xaxis.set_major_formatter(t_form)
        ax2.grid(which='both')
        
        #Sort out the alpha angle axis. 
        alpha_data = ax3.plot(t[0], np.rad2deg(smile.alpha_angle), 'k')[0]
        ax3.set_ylabel('Alpha [deg]')
        #ax2.set_xlabel('UT') 
        ax3.set_xlim(t[0], t[-1]) 
        ax3.set_ylim(-90,90)
        
        #Format time axis. 
        ax3.xaxis.set_major_locator(t_locate_major)
        ax3.xaxis.set_minor_locator(t_locate_minor)
        ax3.xaxis.set_major_formatter(t_form)
        ax3.grid(which='both')
        
        #Sort out the limbc-alpha axis. 
        limb_data = ax4.plot(t[0], np.rad2deg(smile.limb_c - smile.alpha_angle), 'k')[0]
        ax4.set_ylabel("l+r-alpha [deg]")
        #ax2.set_xlabel('UT') 
        ax4.set_xlim(t[0], t[-1]) 
        ax4.set_ylim(-90,90)
        
        #Format time axis. 
        ax4.xaxis.set_major_locator(t_locate_major)
        ax4.xaxis.set_minor_locator(t_locate_minor)
        ax4.xaxis.set_major_formatter(t_form)
        ax4.grid(which='both')
        
        #Sort out the tilt axis. 
        #tilt_data = ax5.plot(t[0], np.rad2deg(smile.sxi_tilt), 'k')[0]
        #ax5.set_ylabel("Tilt [deg]")
        #ax2.set_xlabel('UT') 
        #ax5.set_xlim(t[0], t[-1]) 
        #ax5.set_ylim(-180,190)
        
        #Format time axis. 
        #ax5.xaxis.set_major_locator(t_locate_major)
        #ax5.xaxis.set_minor_locator(t_locate_minor)
        #ax5.xaxis.set_major_formatter(t_form)
        #ax5.grid(which='both')
        
        #Sort out the radius axis. 
        r_data = ax6.plot(t[0], r[0], 'k')[0]
        ax6.set_ylabel("r [RE]")
        #ax2.set_xlabel('UT') 
        ax6.set_xlim(t[0], t[-1]) 
        ax6.set_ylim(0,22)
        
        #Format time axis. 
        ax6.xaxis.set_major_locator(t_locate_major)
        ax6.xaxis.set_minor_locator(t_locate_minor)
        ax6.xaxis.set_major_formatter(t_form)
        ax6.grid(which='both')
        
        #Add line to show where SMILE is too close to Earth to operate. 
        rect1 = Rectangle((t[0],0), t[-1]-t[0], 7.84, color='red', alpha=0.5)
        ax6.add_patch(rect1)
        #ax6.plot([t[0],t[-1]], [7.84,7.84], 'k')
        
        #Add line to show sun pointing constraint. 
        rect2 = Rectangle((t[0],90-35.83), t[-1]-t[0], 35.83, color='red', alpha=0.5)
        ax4.add_patch(rect2)
        #ax4.plot([t[0],t[-1]], [90-35.83,90-35.83], 'k')
        
        
        #Make a list to store the aim data. 
        aim_list = []
        alpha_list = []
        limb_list = [] 
        #tilt_list = []
        
        #This is the function that will update the image. 
        def update(frame):
            print (frame)
            #print (x[frame], y[frame], z[frame])
            #Recalculate SMILE. 
            smile = sfl.smile_limb(n_pixels=8, m_pixels=4, smile_loc=[x[frame], y[frame], z[frame]], p_max=p_max) 
            aim_list.append(smile.target_loc[0]) 
            alpha_list.append(np.rad2deg(smile.alpha_angle))
            limb_list.append(np.rad2deg(smile.limb_c-smile.alpha_angle))
            #tilt_list.append(np.rad2deg(smile.sxi_tilt)) 
            
            #Plot SMILE location in blue. 
            smile_pos.set_data_3d([x[frame]], [y[frame]], [z[frame]])
             
            #Add SMILE vector. 
            smile_vect.set_data_3d([0, x[frame]], [0, y[frame]], [0, z[frame]]) 
            
            #Add SMILE orbit. 
            smile_orbit.set_data_3d(x[:frame], y[:frame], z[:frame])
            
            #Add projection onto planes. 
            smile_x.set_data_3d([xlim[0]], [y[frame]], [z[frame]])
            smile_y.set_data_3d([x[frame]], [ylim[0]], [z[frame]])
            smile_z.set_data_3d([x[frame]], [y[frame]], [zlim[0]])
            
            smile_x_orbit.set_data_3d(np.full(len(x[:frame]),xlim[0]), y[:frame], z[:frame])
            smile_y_orbit.set_data_3d(x[:frame], np.full(len(y[:frame]),ylim[0]), z[:frame])
            smile_z_orbit.set_data_3d(x[:frame], y[:frame], np.full(len(z[:frame]),zlim[0]))
            
            
            #Add Look vector. 
            look_vect.set_data_3d([x[frame], x[frame]+smile.L[0]], [y[frame], y[frame]+smile.L[1]], [z[frame], z[frame]+smile.L[2]])
            
            #Add Aim vector. 
            aim_vect.set_data_3d([0, smile.target_loc[0]], [0, smile.target_loc[1]], [0, smile.target_loc[2]])
            
            #Add b vector. 
            b_vect.set_data_3d([0, smile.b[0]], [0, smile.b[1]], [0, smile.b[2]])
            
            #Add FOV info. 
            c1.set_data_3d(smile.xpos[0][0], smile.ypos[0][0], smile.zpos[0][0])
            c2.set_data_3d(smile.xpos[0][-1], smile.ypos[0][-1], smile.zpos[0][-1])
            c3.set_data_3d(smile.xpos[-1][0], smile.ypos[-1][0], smile.zpos[-1][0])
            c4.set_data_3d(smile.xpos[-1][-1], smile.ypos[-1][-1], smile.zpos[-1][-1])
            
            c5.set_data_3d([smile.xpos[0][0][-1],smile.xpos[0][-1][-1]], [smile.ypos[0][0][-1],smile.ypos[0][-1][-1]], [smile.zpos[0][0][-1],smile.zpos[0][-1][-1]])
            c6.set_data_3d([smile.xpos[0][-1][-1],smile.xpos[-1][-1][-1]], [smile.ypos[0][-1][-1],smile.ypos[-1][-1][-1]], [smile.zpos[0][-1][-1],smile.zpos[-1][-1][-1]])
            c7.set_data_3d([smile.xpos[-1][-1][-1],smile.xpos[-1][0][-1]], [smile.ypos[-1][-1][-1],smile.ypos[-1][0][-1]], [smile.zpos[-1][-1][-1],smile.zpos[-1][0][-1]])
            c8.set_data_3d([smile.xpos[-1][0][-1],smile.xpos[0][0][-1]], [smile.ypos[-1][0][-1],smile.ypos[0][0][-1]], [smile.zpos[-1][0][-1],smile.zpos[0][0][-1]])
            
            #Update the title. 
            title.set_text('SMILE: ({:.2f},{:.2f},{:.2f})\nAim: ({:.2f},{:.2f},{:.2f})\n{} - {} '.format(smile.smile_loc[0], smile.smile_loc[1], smile.smile_loc[2], smile.target_loc[0], smile.target_loc[1], smile.target_loc[2], t[0].strftime('%Y%m%d %H:%M'), t[-1].strftime('%Y%m%d %H:%M')))
            
            #Update the aim data in the bottom panel. 
            aim_data.set_data(t[:frame], aim_list[:frame])
            
            #Update the alpha angle data in the next axis. 
            alpha_data.set_data(t[:frame], alpha_list[:frame])
            
            #Update limb-alpha angle. 
            limb_data.set_data(t[:frame], limb_list[:frame]) 
            
            #Update the tilt angle. 
            #tilt_data.set_data(t[:frame], tilt_list[:frame])
            
            #Update the radial distance of SMILE. 
            r_data.set_data(t[:frame], r[:frame])
            
            return (smile_pos, smile_vect, smile_orbit, smile_x, smile_y, smile_z, smile_x_orbit, smile_y_orbit, smile_z_orbit, look_vect, aim_vect, b_vect, c1, c2, c3, c4, c5, c6, c7, c8, title, aim_data, alpha_data, limb_data, r_data)
            
        #Now make the animation. 
        ani = animation.FuncAnimation(fig=fig, func=update,  frames=len(x), interval=20) 
        ani.save(self.plot_path+'orbital_animation_{}_{}.gif'.format(t[0].strftime('%Y%m%d_%H:%M'), t[-1].strftime('%Y%m%d_%H:%M')))
        print ('Saved: ', self.plot_path+'orbital_animation_{}_{}.gif'.format(t[0].strftime('%Y%m%d_%H:%M'), t[-1].strftime('%Y%m%d_%H:%M')))
        
    def add_fov_boundaries(self, ax2, smile, color='k', lw=2):
        '''This will add the FOV boundaries in black/white. '''
        
        #For corner pixels only. 
        c1 = ax2.plot(smile.xpos[0][0], smile.ypos[0][0], smile.zpos[0][0], color, lw=lw)
        c2 = ax2.plot(smile.xpos[0][-1], smile.ypos[0][-1], smile.zpos[0][-1], color, lw=lw)
        c3 = ax2.plot(smile.xpos[-1][0], smile.ypos[-1][0], smile.zpos[-1][0], color, lw=lw)
        c4 = ax2.plot(smile.xpos[-1][-1], smile.ypos[-1][-1], smile.zpos[-1][-1], color, lw=lw)
        
        #Join corners together. 
        c5 = ax2.plot([smile.xpos[0][0][-1],smile.xpos[0][-1][-1]], [smile.ypos[0][0][-1],smile.ypos[0][-1][-1]], [smile.zpos[0][0][-1],smile.zpos[0][-1][-1]], color, lw=lw)
        c6 = ax2.plot([smile.xpos[0][-1][-1],smile.xpos[-1][-1][-1]], [smile.ypos[0][-1][-1],smile.ypos[-1][-1][-1]], [smile.zpos[0][-1][-1],smile.zpos[-1][-1][-1]], color, lw=lw)
        c7 = ax2.plot([smile.xpos[-1][-1][-1],smile.xpos[-1][0][-1]], [smile.ypos[-1][-1][-1],smile.ypos[-1][0][-1]], [smile.zpos[-1][-1][-1],smile.zpos[-1][0][-1]], color, lw=lw)
        c8 = ax2.plot([smile.xpos[-1][0][-1],smile.xpos[0][0][-1]], [smile.ypos[-1][0][-1],smile.ypos[0][0][-1]], [smile.zpos[-1][0][-1],smile.zpos[0][0][-1]], color, lw=lw)
        
        return c1, c2, c3, c4, c5, c6, c7, c8

#    def add_fov_rectangle(self, ax, smile, color='gray'):
#        '''This will hopefully add a rectangle to the end of the FOV to make its shape clearer.'''
        
#        v1 = [smile.xpos[0][0][-1], smile.ypos[0][0][-1], smile.zpos[0][0][-1]] 
#        v2 = [smile.xpos[0][-1][-1], smile.ypos[0][-1][-1], smile.zpos[0][-1][-1]] 
#        v3 = [smile.xpos[-1][-1][-1], smile.ypos[-1][-1][-1], smile.zpos[-1][-1][-1]] 
#        v4 = [smile.xpos[-1][0][-1], smile.ypos[-1][0][-1], smile.zpos[-1][0][-1]] 
        
#        rects = [[v1, v2, v3, v4, v1]]
#        rect_obj = ax.add_collection3d(Poly3DCollection(rects, color=color, alpha=0.5, edgecolor=None))
#        return rect_obj
                   
#    def add_earth(self, ax):
#        '''This will add a sphere for the Earth. '''
        
        #Create a spherical surface. 
#        radius = 1
#        u = np.linspace(np.pi/2, 1.5*np.pi, 100) 
#        v = np.linspace(0, np.pi, 100) 
#        x = radius* np.outer(np.cos(u), np.sin(v)) 
#        y = radius* np.outer(np.sin(u), np.sin(v))
#        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))
        
#        ax.plot_surface(x, y, z, color='k', lw=0, alpha=1)      
        
#        u = np.linspace(-np.pi/2, np.pi/2, 100) 
#        v = np.linspace(0, np.pi, 100) 
#        x = radius* np.outer(np.cos(u), np.sin(v)) 
#        y = radius* np.outer(np.sin(u), np.sin(v))
#        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))
        
#        ax.plot_surface(x, y, z, color='cyan', lw=0, alpha=1)  
        
    
        

