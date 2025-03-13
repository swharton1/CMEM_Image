#This script will load in one of the new orbits created by Andy and Steve. 

import numpy as np 
import datetime as dt
import matplotlib.pyplot as plt 
from spacepy import coordinates as coord
from spacepy.time import Ticktock 
from matplotlib.patches import Wedge, Polygon, Circle
import os

class orbit():
    '''This will contain functions to read in the orbit data.'''
    
    def __init__(self, orbit_num=1):
        
        self.orbit_num = orbit_num 
        
        self.orbit_path = '/data/smile/shared/sxi_skybgd/' 
        self.filename = 'time_doy_xorb_yorb_zorb_XaimY0Z0_2026_{:0>2}.txt'.format(orbit_num)
        
        self.read_orbit_file()
        self.calc_gsm() 
        
        self.plot_path = os.environ.get("PLOT_PATH")
        
    def read_orbit_file(self):
        '''This actually reads the file and returns the information.'''
        
        print ('Read orbit file...') 
        fullname = self.orbit_path+self.filename
        with open(fullname, 'r') as f:
            lines = f.readlines() 
        
        time = []
        doy = []
        x = []
        y = []
        z = []
        aim = []
        
            
        for l, line in enumerate(lines):
            if l > 0: 
                lsplit = line.split(' ') 
                time.append(float(lsplit[1]))
                doy.append(float(lsplit[3]))
                x.append(float(lsplit[5]))
                y.append(float(lsplit[7]))
                z.append(float(lsplit[9]))
                aim.append(float(lsplit[11])) 
        
        self.data = {}         
        self.data['time'] = np.array(time)
        self.data['doy'] = np.array(doy)
        self.data['x_gse'] = np.array(x)
        self.data['y_gse'] = np.array(y)
        self.data['z_gse'] = np.array(z)
        self.data['aim'] = np.array(aim) 
        
        #Need to get datetime object.
        start = dt.datetime(2026,1,1) 
        self.data['dtime'] = np.array([start + dt.timedelta(days=d) for d in self.data['doy']])


    def calc_gsm(self):
        '''This will work out the GSM positions of the spacecraft.'''
        
        #Need coordinates in appropriate form for conversion.
        print ('Convert coordinates to GSM with spacepy...')
        coords_gse = np.zeros((self.data['x_gse'].size,3))
        coords_gse[:,0] = self.data['x_gse']
        coords_gse[:,1] = self.data['y_gse']
        coords_gse[:,2] = self.data['z_gse']
        
        #Needs to take data in a a list of xyz points. 
        coord_obj = coord.Coords(coords_gse, 'GSE', 'car')
        
        #Add time information. 
        coord_obj.ticks = Ticktock(self.data['dtime'], 'UTC') 
        
        #To convert to GSM. 
        coords_gsm = coord_obj.convert('GSM', 'car') 
    
        self.data['x_gsm'] = coords_gsm.x
        self.data['y_gsm'] = coords_gsm.y
        self.data['z_gsm'] = coords_gsm.z 
        

    def plot_ellipse_2d(self, scatter=False, gse_col='darkblue', gsm_col='r'):
        '''This will plot the xy and xz planes of the orbit in GSE and GSM. 
        
        Parameters
        ----------
        scatter - Boolean to plot a scatter plot or a line plot. 
        '''
        
        self.dtime = self.data['dtime']
        self.x_gse = self.data['x_gse']
        self.y_gse = self.data['y_gse']
        self.z_gse = self.data['z_gse']
        self.x_gsm = self.data['x_gsm']
        self.y_gsm = self.data['y_gsm']
        self.z_gsm = self.data['z_gsm']
        
        
        #Set conditions for inbound and outbound. 
        #outward = np.where(self.phi > 0)
        #inward = np.where(self.phi < 0)
        
        fig = plt.figure(figsize=(6,6))
        fig.subplots_adjust(wspace=0.3)
        
        #GSE
        ax1 = fig.add_subplot(221)
        if scatter:
            ax1.scatter(self.x_gse, self.z_gse, c=gse_col, marker='x', s=5)
            #ax1.scatter(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, c='darkblue', marker='x', s=5)
        else:
            ax1.plot(self.x_gse, self.z_gse, gse_col)
            #ax1.plot(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, 'darkblue')
        
        ax1.set_xlabel(r'$x_{GSE}$')
        ax1.set_ylabel(r'$z_{GSE}$') 
        ax1.set_title('Orbit in GSE')
        self.make_earth_2d(ax1, rotation=-90)
        ax1.set_aspect('equal')
        
        ax2 = fig.add_subplot(223)
        if scatter:
            ax2.scatter(self.x_gse, self.y_gse, c=gse_col, marker='x', s=5)
            #ax2.scatter(self.x_gse, self.y_gse, c='darkblue', marker='x', s=5)
        else:
            ax2.plot(self.x_gse, self.y_gse, gse_col)
            #ax2.plot(self.x_gse, self.y_gse, 'darkblue')
        ax2.set_xlabel(r'$x_{GSE}$')
        ax2.set_ylabel(r'$y_{GSE}$') 
        #ax2.set_title('Orbit in GSE')
        self.make_earth_2d(ax2, rotation=-90)
        ax2.set_aspect('equal')
        
        ax3 = fig.add_subplot(222)
        if scatter:
            ax3.scatter(self.x_gsm, self.y_gsm, c=gsm_col, marker='x', s=5)
            #ax3.scatter(self.x_gsm, self.z_gsm, c='darkblue', marker='x', s=5)
        else:
            ax3.plot(self.x_gsm, self.z_gsm, gsm_col)
            #ax3.plot(self.x_gsm, self.z_gsm, 'darkblue')
        ax3.set_xlabel(r'$x_{GSM}$')
        ax3.set_ylabel(r'$z_{GSM}$') 
        ax3.set_title('Orbit in GSM')
        self.make_earth_2d(ax3, rotation=-90)
        ax3.set_aspect('equal')
        
        ax4 = fig.add_subplot(224)
        if scatter:
            ax4.scatter(self.x_gsm, self.y_gsm, c=gsm_col, marker='x', s=5)
            #ax4.scatter(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, c='darkblue', marker='x', s=5)
        else:
            ax4.plot(self.x_gsm, self.y_gsm, gsm_col)
            #ax4.plot(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, 'darkblue')
        ax4.set_xlabel(r'$x_{GSM}$')
        ax4.set_ylabel(r'$y_{GSM}$') 
        #ax4.set_title('Orbit in GSM')
        self.make_earth_2d(ax4, rotation=-90)
        ax4.set_aspect('equal')
        
        #fig.text(0.95, 0.05, 'Outbound', color='cyan', fontsize=8, ha='right')
        #fig.text(0.95, 0.03, 'Inbound', color='darkblue', fontsize=8, ha='right')
        #fig.text(0.05, 0.06, 'Inclination = {:.1f} deg'.format(np.rad2deg(self.inc)), fontsize=8, ha='left')
        #fig.text(0.05, 0.04, 'Arg. Periapsis = {:.1f} deg'.format(np.rad2deg(self.omega)), fontsize=8, ha='left')
        #fig.text(0.05, 0.02, 'RAAN* = {:.1f} deg'.format(np.rad2deg(self.raan)), fontsize=8, ha='left') 
        
        fig.text(0.5,0.95,'{} - {}'.format(self.dtime[0].strftime("%Y%m%d %H:%M"),self.dtime[-1].strftime("%Y%m%d %H:%M")), ha='center', fontsize=12)
        
        #Save the figure. 
        print('Saved: ',self.plot_path+'new_orbits/orbit_gse_gsm_{:0>2}.png'.format(self.orbit_num))
        fig.savefig(self.plot_path+'new_orbits/orbit_gse_gsm_{:0>2}.png'.format(self.orbit_num))
        
        
    def make_earth_2d(self, ax, rotation=0):
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
        
        
        
def compare_aims(orbit_num=1):
    '''This will compare the aim from my calculation with that in Andy's file.'''
    
    from . import smile_fov_limb 
    
    #Get information in orbit file. 
    orbit_info = orbit(orbit_num=orbit_num)  
    
    #Loop through each orbital position for just a few positions. 
    for i in range(len(orbit_info.data['x_gse'])):
    
        #Get position. 
        smile_loc = (orbit_info.data['x_gse'][i], orbit_info.data['y_gse'][i], orbit_info.data['z_gse'][i])
        print ('SMILE Loc: ', smile_loc)
        
        #Get aim in file. 
        aim_file = orbit_info.data['aim'][i]
        
        #Make SMILE object.
        smile = smile_fov_limb.smile_limb(n_pixels=20, m_pixels=10, smile_loc=smile_loc)
        
        print ('Aim File: ', aim_file, 'Aim Sam: ', smile.target_loc[0], 'Diff: ', aim_file-smile.target_loc[0])
       
